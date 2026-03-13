/* generate_data.cpp — Numerical experiments for nonstandard FD scheme analysis.
 *
 * Links against the QuantLib build via CMake find_package(QuantLib).
 * Outputs CSV files to ../data/ (relative to cwd) for consumption by
 * plot_figures.py.
 *
 * Compile with C++17.  Run from the build/ directory:
 *     ./generate_data            # writes to ../data/
 *     ./generate_data /tmp/out   # writes to /tmp/out/
 */

#include <ql/qldefines.hpp>
#include <ql/exercise.hpp>
#include <ql/instruments/payoffs.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/matrixutilities/sparsematrix.hpp>
#include <ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmeshercomposite.hpp>
#include <ql/methods/finitedifferences/meshers/uniform1dmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesop.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/operators/fdmhyperboliccot.hpp>
#include <ql/methods/finitedifferences/operators/firstderivativeop.hpp>
#include <ql/methods/finitedifferences/operators/secondderivativeop.hpp>
#include <ql/methods/finitedifferences/operators/triplebandlinearop.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp>
#include <ql/methods/finitedifferences/solvers/fdmsolverdesc.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmdiscretebarrierstepcondition.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp>
#include <ql/methods/finitedifferences/utilities/fdminnervaluecalculator.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/settings.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace QuantLib;

// =====================================================================
// CSV writer with metadata
// =====================================================================

class CsvWriter {
  public:
    CsvWriter(const std::string& path) : ofs_(path) {
        if (!ofs_)
            throw std::runtime_error("Cannot open " + path);
        ofs_ << std::setprecision(15);
    }

    void meta(const std::string& key, const std::string& value) {
        ofs_ << "# " << key << ": " << value << "\n";
    }

    void meta(const std::string& key, Real value) {
        ofs_ << "# " << key << ": " << value << "\n";
    }

    void meta(const std::string& key, Size value) {
        ofs_ << "# " << key << ": " << value << "\n";
    }

    template <typename... Args>
    void header(Args&&... cols) {
        writeRow(std::forward<Args>(cols)...);
    }

    template <typename... Args>
    void row(Args&&... vals) {
        writeRow(std::forward<Args>(vals)...);
    }

  private:
    std::ofstream ofs_;

    template <typename T>
    void writeOne(const T& v) { ofs_ << v; }

    void writeOne(Real v) { ofs_ << v; }

    template <typename T, typename... Rest>
    void writeRow(const T& first, Rest&&... rest) {
        writeOne(first);
        ((ofs_ << "," , writeOne(rest)), ...);
        ofs_ << "\n";
    }

    template <typename T>
    void writeRow(const T& only) {
        writeOne(only);
        ofs_ << "\n";
    }
};

// =====================================================================
// Helpers
// =====================================================================

static ext::shared_ptr<GeneralizedBlackScholesProcess>
makeProcess(Real spot, Rate r, Rate q, Volatility vol) {
    DayCounter dc = Actual365Fixed();
    Date today = Settings::instance().evaluationDate();
    return ext::make_shared<GeneralizedBlackScholesProcess>(
        Handle<Quote>(ext::make_shared<SimpleQuote>(spot)),
        Handle<YieldTermStructure>(
            ext::make_shared<FlatForward>(today, q, dc)),
        Handle<YieldTermStructure>(
            ext::make_shared<FlatForward>(today, r, dc)),
        Handle<BlackVolTermStructure>(
            ext::make_shared<BlackConstantVol>(today, NullCalendar(), vol, dc)));
}

static std::string schemeName(FdmBlackScholesSpatialDesc::Scheme s) {
    switch (s) {
      case FdmBlackScholesSpatialDesc::Scheme::StandardCentral:
        return "StandardCentral";
      case FdmBlackScholesSpatialDesc::Scheme::ExponentialFitting:
        return "ExponentialFitting";
      case FdmBlackScholesSpatialDesc::Scheme::MilevTaglianiCNEffectiveDiffusion:
        return "MilevTaglianiCN";
      default:
        return "Unknown";
    }
}

static std::string policyName(FdmBlackScholesSpatialDesc::MMatrixPolicy p) {
    switch (p) {
      case FdmBlackScholesSpatialDesc::MMatrixPolicy::None:
        return "None";
      case FdmBlackScholesSpatialDesc::MMatrixPolicy::DiagnosticsOnly:
        return "DiagnosticsOnly";
      case FdmBlackScholesSpatialDesc::MMatrixPolicy::FailFast:
        return "FailFast";
      case FdmBlackScholesSpatialDesc::MMatrixPolicy::FallbackToExponentialFitting:
        return "FallbackToExponentialFitting";
      default:
        return "Unknown";
    }
}

static FdmBlackScholesSpatialDesc descForScheme(
        FdmBlackScholesSpatialDesc::Scheme s,
        FdmBlackScholesSpatialDesc::MMatrixPolicy policy =
            FdmBlackScholesSpatialDesc::MMatrixPolicy::None) {
    FdmBlackScholesSpatialDesc d;
    d.scheme = s;
    d.mMatrixPolicy = policy;
    return d;
}

// Write the standard metadata block to a CSV for solver-based experiments.
static void writeStandardMeta(CsvWriter& csv,
                              const std::string& schemeName_,
                              const std::string& effectiveScheme,
                              Size xGrid, Size tGrid,
                              Rate r, Rate q, Volatility sigma,
                              Real strike, Time maturity,
                              const std::string& mMatrixPol,
                              const std::string& mesh) {
    csv.meta("scheme", schemeName_);
    csv.meta("effective_scheme", effectiveScheme);
    csv.meta("xGrid", xGrid);
    csv.meta("tGrid", tGrid);
    csv.meta("r", r);
    csv.meta("q", q);
    csv.meta("sigma", sigma);
    csv.meta("strike", strike);
    csv.meta("maturity", maturity);
    csv.meta("mMatrixPolicy", mMatrixPol);
    csv.meta("mesh", mesh);
}

// Black-Scholes analytical call price
static Real bsCall(Real S, Real K, Rate r, Rate q,
                    Volatility vol, Time T) {
    if (T <= 0.0) return std::max(S - K, 0.0);
    Real d1 = (std::log(S / K) + (r - q + 0.5 * vol * vol) * T)
              / (vol * std::sqrt(T));
    Real d2 = d1 - vol * std::sqrt(T);
    CumulativeNormalDistribution N;
    return S * std::exp(-q * T) * N(d1)
         - K * std::exp(-r * T) * N(d2);
}

// Look up the value of a solution array at a given spot level
// by finding the nearest grid node.
static Real valueAtSpot(const Array& v,
                        const ext::shared_ptr<FdmMesher>& mesher,
                        Real spotLevel, Size direction = 0) {
    const Real xTarget = std::log(spotLevel);
    Size bestIdx = 0;
    Real bestDist = 1e30;
    for (const auto& iter : *mesher->layout()) {
        const Real x = mesher->location(iter, direction);
        const Real d = std::fabs(x - xTarget);
        if (d < bestDist) {
            bestDist = d;
            bestIdx = iter.index();
        }
    }
    return v[bestIdx];
}

// Build a payoff vector on the grid.
static Array buildPayoff(const ext::shared_ptr<FdmMesher>& mesher,
                         const ext::shared_ptr<Payoff>& payoff,
                         Size direction = 0) {
    Array rhs(mesher->layout()->size());
    for (const auto& iter : *mesher->layout()) {
        const Real S = std::exp(mesher->location(iter, direction));
        rhs[iter.index()] = (*payoff)(S);
    }
    return rhs;
}

// Truncated call payoff: max(S-K, 0) * 1_{[K,U]}(S)
class TruncatedCallPayoff : public Payoff {
  public:
    TruncatedCallPayoff(Real strike, Real upper)
    : strike_(strike), upper_(upper) {}
    std::string name() const override { return "TruncatedCallPayoff"; }
    std::string description() const override { return "TruncatedCallPayoff"; }
    Real operator()(Real s) const override {
        return (s >= strike_ && s <= upper_) ? s - strike_ : 0.0;
    }
  private:
    Real strike_, upper_;
};

// =====================================================================
// Experiment 1-2: Truncated call — one CSV per scheme (Figs 1-2)
// =====================================================================

void runTruncatedCall(const std::string& dataDir) {
    std::cout << "  Running truncated call experiments..." << std::flush;

    const Rate r = 0.05;
    const Rate q = 0.0;
    const Volatility vol = 0.001;
    const Real K = 50.0;
    const Real U = 70.0;
    const Time T = 5.0 / 12.0;
    const Real spot = 60.0;

    const Real xMin = std::log(1.0);
    const Real xMax = std::log(140.0);
    const Size xGrid = 2801;
    const Size tGrid = 2801;

    auto process = makeProcess(spot, r, q, vol);
    auto payoff = ext::make_shared<TruncatedCallPayoff>(K, U);

    using Scheme = FdmBlackScholesSpatialDesc::Scheme;
    Scheme schemes[] = {Scheme::StandardCentral,
                        Scheme::ExponentialFitting,
                        Scheme::MilevTaglianiCNEffectiveDiffusion};

    auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    const Array locs = mesher->locations(0);

    for (Size s = 0; s < 3; ++s) {
        auto desc = descForScheme(schemes[s]);
        auto calc = ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);
        auto conditions = ext::make_shared<FdmStepConditionComposite>(
            std::list<std::vector<Time>>(),
            FdmStepConditionComposite::Conditions());

        FdmSolverDesc solverDesc = {
            mesher, FdmBoundaryConditionSet(), conditions,
            calc, T, tGrid, 0};

        FdmBlackScholesSolver solver(
            Handle<GeneralizedBlackScholesProcess>(process),
            K, solverDesc,
            FdmSchemeDesc::CrankNicolson(),
            false, -Null<Real>(),
            Handle<FdmQuantoHelper>(), desc);

        std::string schName = schemeName(schemes[s]);
        CsvWriter csv(dataDir + "/truncated_call_" + schName + ".csv");
        writeStandardMeta(csv, schName, schName,
                          xGrid, tGrid, r, q, vol, K, T,
                          "None", "Uniform1dMesher");
        csv.meta("upper_barrier", U);
        csv.header("S", "price");

        for (Size i = 0; i < xGrid; ++i) {
            Real S = std::exp(locs[i]);
            csv.row(S, solver.valueAt(S));
        }
    }

    // Fine-grid reference using ExponentialFitting at 8x refinement
    {
        const Size fineX = xGrid * 8;
        const Size fineT = tGrid * 8;

        auto fineMesher = ext::make_shared<FdmMesherComposite>(
            ext::make_shared<Uniform1dMesher>(xMin, xMax, fineX));
        auto fineCalc = ext::make_shared<FdmLogInnerValue>(payoff, fineMesher, 0);
        auto fineCond = ext::make_shared<FdmStepConditionComposite>(
            std::list<std::vector<Time>>(),
            FdmStepConditionComposite::Conditions());
        auto fineDesc = descForScheme(
            Scheme::ExponentialFitting,
            FdmBlackScholesSpatialDesc::MMatrixPolicy::FallbackToExponentialFitting);

        FdmSolverDesc fineSolverDesc = {
            fineMesher, FdmBoundaryConditionSet(), fineCond,
            fineCalc, T, fineT, 0};

        FdmBlackScholesSolver fineSolver(
            Handle<GeneralizedBlackScholesProcess>(process),
            K, fineSolverDesc,
            FdmSchemeDesc::CrankNicolson(),
            false, -Null<Real>(),
            Handle<FdmQuantoHelper>(), fineDesc);

        CsvWriter csv(dataDir + "/truncated_call_reference.csv");
        writeStandardMeta(csv, "ExponentialFitting", "ExponentialFitting",
                          fineX, fineT, r, q, vol, K, T,
                          "FallbackToExponentialFitting", "Uniform1dMesher");
        csv.meta("upper_barrier", U);
        csv.meta("refinement_factor", Size(8));
        csv.header("S", "price");

        // Output on the original coarse grid spots for comparison
        for (Size i = 0; i < xGrid; ++i) {
            Real S = std::exp(locs[i]);
            csv.row(S, fineSolver.valueAt(S));
        }
    }

    std::cout << " done.\n";
}

// =====================================================================
// Experiment 5: Grid convergence (Fig 5)
// =====================================================================

void runGridConvergence(const std::string& dataDir) {
    std::cout << "  Running grid convergence study..." << std::flush;

    const Real spot = 100.0;
    const Rate r = 0.05;
    const Rate q = 0.02;
    const Volatility vol = 0.20;
    const Real K = 100.0;
    const Time T = 1.0;

    const Real refPrice = bsCall(spot, K, r, q, vol, T);
    auto process = makeProcess(spot, r, q, vol);

    using Scheme = FdmBlackScholesSpatialDesc::Scheme;
    Scheme schemes[] = {Scheme::StandardCentral,
                        Scheme::ExponentialFitting,
                        Scheme::MilevTaglianiCNEffectiveDiffusion};

    Size grids[] = {25, 50, 100, 200, 400, 800, 1600};

    ext::shared_ptr<StrikedTypePayoff> payoff(
        new PlainVanillaPayoff(Option::Call, K));

    for (Size si = 0; si < 3; ++si) {
        std::string schName = schemeName(schemes[si]);
        CsvWriter csv(dataDir + "/grid_convergence_" + schName + ".csv");
        writeStandardMeta(csv, schName, schName,
                          Size(0), Size(0), r, q, vol, K, T,
                          "None", "FdmBlackScholesMesher");
        csv.meta("spot", spot);
        csv.meta("reference_price", refPrice);
        csv.header("xGrid", "tGrid", "price", "error");

        for (Size xGrid : grids) {
            Size tGrid = 4 * xGrid;

            auto desc = descForScheme(schemes[si]);
            auto equityMesher = ext::make_shared<FdmBlackScholesMesher>(
                xGrid, process, T, K);
            auto mesher = ext::make_shared<FdmMesherComposite>(equityMesher);
            auto calc = ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);
            auto conditions = ext::make_shared<FdmStepConditionComposite>(
                std::list<std::vector<Time>>(),
                FdmStepConditionComposite::Conditions());

            FdmSolverDesc solverDesc = {
                mesher, FdmBoundaryConditionSet(), conditions,
                calc, T, tGrid, 0};

            FdmBlackScholesSolver solver(
                Handle<GeneralizedBlackScholesProcess>(process),
                K, solverDesc,
                FdmSchemeDesc::CrankNicolson(),
                false, -Null<Real>(),
                Handle<FdmQuantoHelper>(), desc);

            Real price = solver.valueAt(spot);

            csv.row(xGrid, tGrid,
                    price, std::fabs(price - refPrice));
        }
    }

    std::cout << " done.\n";
}

// =====================================================================
// Experiment 6-7: Operator diagnostics (Figs 6-7)
// =====================================================================

void runOperatorDiagnostics(const std::string& dataDir) {
    std::cout << "  Running operator diagnostics..." << std::flush;

    const Real spot = 100.0;
    const Rate r = 0.05;
    const Rate q = 0.0;
    const Volatility vol = 0.001;
    const Real K = 100.0;
    const Time T = 1.0;

    const Real xMin = std::log(50.0);
    const Real xMax = std::log(200.0);
    const Size xGrid = 200;
    const Real h = (xMax - xMin) / (xGrid - 1);
    const Time dt = T / 50.0;

    auto process = makeProcess(spot, r, q, vol);
    auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    using Scheme = FdmBlackScholesSpatialDesc::Scheme;
    Scheme schemes[] = {Scheme::StandardCentral,
                        Scheme::ExponentialFitting,
                        Scheme::MilevTaglianiCNEffectiveDiffusion};

    // Off-diagonal data (Fig 7) — one CSV per scheme
    for (Size si = 0; si < 3; ++si) {
        std::string schName = schemeName(schemes[si]);
        CsvWriter csv(dataDir + "/mmatrix_offdiag_" + schName + ".csv");
        writeStandardMeta(csv, schName, schName,
                          xGrid, Size(0), r, q, vol, K, T,
                          "None", "Uniform1dMesher");
        csv.header("node", "S", "lower", "upper");

        auto desc = descForScheme(schemes[si]);
        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), desc);
        op->setTime(0.0, dt);
        SparseMatrix mat = op->toMatrix();

        for (Size i = 1; i + 1 < xGrid; ++i) {
            Real S = std::exp(xMin + i * h);
            csv.row(i, S,
                    Real(mat(i, i-1)), Real(mat(i, i+1)));
        }
    }

    // Effective diffusion data (Fig 6) — one CSV per scheme
    for (Size si = 0; si < 3; ++si) {
        std::string schName = schemeName(schemes[si]);
        CsvWriter csv(dataDir + "/effective_diffusion_" + schName + ".csv");
        writeStandardMeta(csv, schName, schName,
                          xGrid, Size(0), r, q, vol, K, T,
                          "None", "Uniform1dMesher");
        csv.meta("h", h);
        csv.header("node", "S", "aUsed");

        auto desc = descForScheme(schemes[si]);
        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), desc);
        op->setTime(0.0, dt);
        SparseMatrix mat = op->toMatrix();

        for (Size i = 1; i + 1 < xGrid; ++i) {
            Real S = std::exp(xMin + i * h);
            Real lower = Real(mat(i, i-1));
            Real upper = Real(mat(i, i+1));
            Real aUsed = (lower + upper) * h * h / 2.0;
            csv.row(i, S, aUsed);
        }
    }

    std::cout << " done.\n";
}

// =====================================================================
// Experiment 9: xCothx / Peclet number (Fig 9)
// =====================================================================

void runXCothx(const std::string& dataDir) {
    std::cout << "  Running xCothx data generation..." << std::flush;

    CsvWriter csv(dataDir + "/xcothx.csv");
    csv.meta("scheme", "N/A");
    csv.meta("effective_scheme", "N/A");
    csv.meta("xGrid", "N/A");
    csv.meta("tGrid", "N/A");
    csv.meta("r", "N/A");
    csv.meta("q", "N/A");
    csv.meta("sigma", "N/A");
    csv.meta("strike", "N/A");
    csv.meta("maturity", "N/A");
    csv.meta("mMatrixPolicy", "N/A");
    csv.meta("mesh", "N/A");
    csv.meta("experiment", "xcothx");
    csv.meta("xSmall", Real(1e-6));
    csv.meta("xLarge", Real(50.0));
    csv.header("Pe", "xCothx", "abs_Pe");

    // Dense sampling across regimes
    std::vector<Real> peValues;
    // Fine near zero
    for (Real pe = -1.0; pe <= 1.0; pe += 0.001)
        peValues.push_back(pe);
    // Moderate range
    for (Real pe = -50.0; pe < -1.0; pe += 0.1)
        peValues.push_back(pe);
    for (Real pe = 1.0; pe <= 50.0; pe += 0.1)
        peValues.push_back(pe);
    // Extreme values
    for (Real pe = -100.0; pe < -50.0; pe += 1.0)
        peValues.push_back(pe);
    for (Real pe = 50.0; pe <= 100.0; pe += 1.0)
        peValues.push_back(pe);

    std::sort(peValues.begin(), peValues.end());

    for (Real pe : peValues) {
        Real val = detail::xCothx(pe, 1e-6, 50.0);
        csv.row(pe, val, std::fabs(pe));
    }

    std::cout << " done.\n";
}

// =====================================================================
// Experiment 8: Performance benchmark (Fig 8)
// =====================================================================

void runBenchmark(const std::string& dataDir) {
    std::cout << "  Running performance benchmark..." << std::flush;

    const Real spot = 100.0;
    const Rate r = 0.05;
    const Rate q = 0.02;
    const Volatility vol = 0.20;
    const Real K = 100.0;
    const Time T = 1.0;

    const Real refPrice = bsCall(spot, K, r, q, vol, T);
    auto process = makeProcess(spot, r, q, vol);

    using Scheme = FdmBlackScholesSpatialDesc::Scheme;
    Scheme schemes[] = {Scheme::StandardCentral,
                        Scheme::ExponentialFitting,
                        Scheme::MilevTaglianiCNEffectiveDiffusion};

    Size grids[] = {50, 100, 200, 400, 800, 1600};
    // No wall-clock timing — use deterministic cost proxy (xGrid*tGrid)

    ext::shared_ptr<StrikedTypePayoff> payoff(
        new PlainVanillaPayoff(Option::Call, K));

    // One CSV per scheme.  Use mMatrixPolicy=None for benchmark so
    // scheme == effective_scheme (no ambiguous fallback provenance).
    for (Size si = 0; si < 3; ++si) {
        std::string schName = schemeName(schemes[si]);

        CsvWriter csv(dataDir + "/benchmark_" + schName + ".csv");
        writeStandardMeta(csv, schName, schName,
                          Size(0), Size(0), r, q, vol, K, T,
                          "None", "FdmBlackScholesMesher");
        csv.meta("spot", spot);
        csv.meta("reference_price", refPrice);
        csv.meta("cost_note", "relative_cost = xGrid * tGrid (deterministic proxy for runtime)");
        csv.header("xGrid", "tGrid", "relative_cost", "price", "error");

        for (Size xGrid : grids) {
            Size tGrid = 4 * xGrid;

            auto desc = descForScheme(schemes[si]);

            auto equityMesher = ext::make_shared<FdmBlackScholesMesher>(
                xGrid, process, T, K);
            auto mesher = ext::make_shared<FdmMesherComposite>(equityMesher);
            auto calc = ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);
            auto conditions = ext::make_shared<FdmStepConditionComposite>(
                std::list<std::vector<Time>>(),
                FdmStepConditionComposite::Conditions());

            FdmSolverDesc solverDesc = {
                mesher, FdmBoundaryConditionSet(), conditions,
                calc, T, tGrid, 0};

            FdmBlackScholesSolver solver(
                Handle<GeneralizedBlackScholesProcess>(process),
                K, solverDesc,
                FdmSchemeDesc::CrankNicolson(),
                false, -Null<Real>(),
                Handle<FdmQuantoHelper>(), desc);

            Real price = solver.valueAt(spot);

            // Use xGrid*tGrid as a deterministic cost proxy
            csv.row(xGrid, tGrid, Size(xGrid) * Size(tGrid),
                    price, std::fabs(price - refPrice));
        }
    }

    std::cout << " done.\n";
}

// =====================================================================
// Experiment 3-4: Discrete double barrier knock-out (Figs 3-4)
//
// Uses FdmDiscreteBarrierStepCondition + FdmStepConditionComposite
// + FdmBackwardSolver with a single rollback call.
// =====================================================================

void runDiscreteBarrier(const std::string& dataDir,
                        const std::string& label,
                        Rate r, Rate q, Volatility vol,
                        Real K, Real L, Real U, Time T,
                        Size nMonitoring, Size xGrid, Size tGrid,
                        bool useUniformMesh = false) {
    std::cout << "  Running discrete barrier (" << label << ")..." << std::flush;

    auto process = makeProcess(100.0, r, q, vol);

    using Scheme = FdmBlackScholesSpatialDesc::Scheme;
    Scheme schemes[] = {Scheme::StandardCentral,
                        Scheme::ExponentialFitting,
                        Scheme::MilevTaglianiCNEffectiveDiffusion};

    // Monitoring dates (equally spaced, including maturity)
    std::vector<Time> monitoringTimes;
    for (Size i = 1; i <= nMonitoring; ++i)
        monitoringTimes.push_back(T * Real(i) / Real(nMonitoring));

    // Spot values to evaluate
    std::vector<Real> spots;
    for (Real S = L - 5.0; S <= U + 5.0; S += 0.1)
        spots.push_back(S);
    // Also include paper Table 1 spot values if moderate-vol
    if (vol > 0.1) {
        for (Real S : {95.0, 95.0001, 95.5, 99.5, 100.0, 100.5,
                       109.5, 109.9999, 110.0})
            spots.push_back(S);
    }
    std::sort(spots.begin(), spots.end());
    spots.erase(std::unique(spots.begin(), spots.end(),
        [](Real a, Real b){ return std::fabs(a-b) < 1e-8; }),
        spots.end());

    ext::shared_ptr<StrikedTypePayoff> payoff(
        new PlainVanillaPayoff(Option::Call, K));

    // At low vol, FdmBlackScholesMesher auto-domain is too narrow to include
    // the barriers.  Use Uniform1dMesher with explicit bounds instead.
    std::string meshType;
    ext::shared_ptr<FdmMesherComposite> mesher;
    if (useUniformMesh) {
        const Real xMin = std::log(L - 15.0);  // well below L
        const Real xMax = std::log(U + 20.0);  // well above U
        mesher = ext::make_shared<FdmMesherComposite>(
            ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));
        meshType = "Uniform1dMesher";
    } else {
        std::vector<std::tuple<Real, Real, bool>> cPoints = {
            {K, 0.1, true}, {L, 0.1, true}, {U, 0.1, true}
        };
        mesher = ext::make_shared<FdmMesherComposite>(
            ext::make_shared<FdmBlackScholesMesher>(
                xGrid, process, T, K,
                Null<Real>(), Null<Real>(), 0.0001, 1.5, cPoints));
        meshType = "FdmBlackScholesMesher";
    }

    // Step condition: discrete double barrier knock-out
    auto barrierCondition =
        ext::make_shared<FdmDiscreteBarrierStepCondition>(
            mesher, monitoringTimes, L, U);

    auto conditions = ext::make_shared<FdmStepConditionComposite>(
        std::list<std::vector<Time>>(1, monitoringTimes),
        FdmStepConditionComposite::Conditions(1, barrierCondition));

    const FdmBoundaryConditionSet bcSet;

    // One CSV per scheme
    for (Size si = 0; si < 3; ++si) {
        auto desc = descForScheme(schemes[si]);
        std::string schName = schemeName(schemes[si]);

        CsvWriter csv(dataDir + "/barrier_" + label + "_" + schName + ".csv");
        writeStandardMeta(csv, schName, schName,
                          xGrid, tGrid, r, q, vol, K, T,
                          "None", meshType);
        csv.meta("lower_barrier", L);
        csv.meta("upper_barrier", U);
        csv.meta("n_monitoring", nMonitoring);
        csv.header("S", "price");

        // Create operator
        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), Size(0),
            ext::shared_ptr<FdmQuantoHelper>(), desc);

        // Build payoff vector
        Array rhs = buildPayoff(mesher, payoff);

        // Single rollback call — the step condition handles barrier
        // application automatically at monitoring times.
        FdmBackwardSolver solver(op, bcSet, conditions,
                                 FdmSchemeDesc::CrankNicolson());
        solver.rollback(rhs, T, 0.0, tGrid, 0);

        // Look up values at requested spots via valueAtSpot
        for (Size j = 0; j < spots.size(); ++j) {
            csv.row(spots[j], valueAtSpot(rhs, mesher, spots[j]));
        }
    }

    std::cout << " done.\n";
}

// =====================================================================
// Experiment MC: Monte Carlo reference for barriers
// =====================================================================

void runBarrierMC(const std::string& dataDir,
                  const std::string& label,
                  Rate r, Rate q, Volatility vol,
                  Real K, Real L, Real U, Time T,
                  Size nMonitoring, Size nPaths) {
    std::cout << "  Running MC barrier reference (" << label
              << ", " << nPaths << " paths)..." << std::flush;

    std::vector<Real> spots;
    if (vol > 0.1) {
        spots = {95.0, 95.0001, 95.5, 99.5, 100.0, 100.5,
                 109.5, 109.9999, 110.0};
    } else {
        for (Real S = L; S <= U; S += 1.0)
            spots.push_back(S);
    }

    CsvWriter csv(dataDir + "/mc_barrier_" + label + ".csv");
    csv.meta("scheme", "N/A");
    csv.meta("effective_scheme", "N/A");
    csv.meta("xGrid", "N/A");
    csv.meta("tGrid", "N/A");
    csv.meta("r", r);
    csv.meta("q", q);
    csv.meta("sigma", vol);
    csv.meta("strike", K);
    csv.meta("maturity", T);
    csv.meta("mMatrixPolicy", "N/A");
    csv.meta("mesh", "N/A");
    csv.meta("lower_barrier", L);
    csv.meta("upper_barrier", U);
    csv.meta("n_monitoring", nMonitoring);
    csv.meta("num_paths", nPaths);
    csv.meta("min_price_threshold", Real(0.01));
    csv.header("S0", "price", "standard_error", "num_paths");

    Real dt = T / nMonitoring;
    Real drift = (r - q - 0.5 * vol * vol) * dt;
    Real diffusion = vol * std::sqrt(dt);
    Real discount = std::exp(-r * T);

    std::mt19937_64 rng(42);
    std::normal_distribution<Real> normal(0.0, 1.0);

    for (Real S0 : spots) {
        Real sumPayoff = 0.0;
        Real sumPayoff2 = 0.0;

        for (Size path = 0; path < nPaths; ++path) {
            Real S = S0;
            bool knocked = false;

            for (Size step = 0; step < nMonitoring; ++step) {
                S *= std::exp(drift + diffusion * normal(rng));
                if (S <= L || S >= U) {
                    knocked = true;
                    break;
                }
            }

            Real payoff = knocked ? 0.0 : std::max(S - K, 0.0);
            payoff *= discount;
            sumPayoff += payoff;
            sumPayoff2 += payoff * payoff;
        }

        Real mean = sumPayoff / nPaths;
        Real var = sumPayoff2 / nPaths - mean * mean;
        Real se = std::sqrt(var / nPaths);

        // Only report spots where price > 0.01 (skip near-zero tails
        // where relative SE is ill-defined)
        if (mean > 0.01) {
            csv.row(S0, mean, se, nPaths);
        }
    }

    std::cout << " done.\n";
}

// =====================================================================
// Main
// =====================================================================

int main(int argc, char* argv[]) {
    try {
        Date today(28, March, 2004);
        Settings::instance().evaluationDate() = today;

        // Output directory: command-line arg or default to ../data
        std::string dataDir = "../data";
        if (argc > 1) {
            dataDir = argv[1];
        }

        // Create output directory
        std::string mkdirCmd = "mkdir -p " + dataDir;
        std::system(mkdirCmd.c_str());

        std::cout << "Generating results data (output: " << dataDir << ")...\n";

        // Fig 9: xCothx/Peclet (fastest, no solver)
        runXCothx(dataDir);

        // Figs 6-7: Operator diagnostics
        runOperatorDiagnostics(dataDir);

        // Fig 5: Grid convergence
        runGridConvergence(dataDir);

        // Figs 1-2: Truncated call (with fine-grid reference)
        runTruncatedCall(dataDir);

        // Fig 8: Performance benchmark
        runBenchmark(dataDir);

        // Fig 3: Discrete barrier moderate-vol (FdmBlackScholesMesher, 4000 nodes)
        runDiscreteBarrier(dataDir, "moderate_vol",
                           0.05, 0.0, 0.25,
                           100.0, 95.0, 110.0, 0.5,
                           5, 4000, 2000,
                           /*useUniformMesh=*/false);

        // Fine-grid barrier reference for moderate-vol (16000 nodes, ExpFit)
        runDiscreteBarrier(dataDir, "moderate_vol_reference",
                           0.05, 0.0, 0.25,
                           100.0, 95.0, 110.0, 0.5,
                           5, 16000, 8000,
                           /*useUniformMesh=*/false);

        // Fig 4: Discrete barrier low-vol (Uniform1dMesher — auto domain too narrow)
        runDiscreteBarrier(dataDir, "low_vol",
                           0.05, 0.0, 0.001,
                           100.0, 95.0, 110.0, 1.0,
                           5, 800, 200,
                           /*useUniformMesh=*/true);

        // MC reference: moderate-vol barrier (5M paths)
        runBarrierMC(dataDir, "moderate_vol",
                     0.05, 0.0, 0.25,
                     100.0, 95.0, 110.0, 0.5,
                     5, 5000000);

        // MC reference: low-vol barrier (5M paths)
        runBarrierMC(dataDir, "low_vol",
                     0.05, 0.0, 0.001,
                     100.0, 95.0, 110.0, 1.0,
                     5, 5000000);

        std::cout << "All experiments complete. CSV files in " << dataDir << "/\n";
        return 0;
    }
    catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
