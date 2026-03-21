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
#include <ql/instruments/payoffs.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/matrixutilities/sparsematrix.hpp>
#include <ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmeshercomposite.hpp>
#include <ql/methods/finitedifferences/meshers/uniform1dmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesop.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/operators/fdmhyperboliccot.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp>
#include <ql/methods/finitedifferences/solvers/fdmsolverdesc.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmdiscretebarrierstepcondition.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp>
#include <ql/methods/finitedifferences/utilities/fdminnervaluecalculator.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/settings.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <random>
#include <string>
#include <vector>

using namespace QuantLib;

using Scheme = FdmBlackScholesSpatialDesc::Scheme;
using Policy = FdmBlackScholesSpatialDesc::MMatrixPolicy;

static const Scheme ALL_SCHEMES[] = {
    Scheme::StandardCentral,
    Scheme::ExponentialFitting,
    Scheme::MilevTaglianiCNEffectiveDiffusion
};

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

    template <typename T, typename... Rest>
    void row(const T& first, Rest&&... rest) {
        ofs_ << first;
        ((ofs_ << "," << rest), ...);
        ofs_ << "\n";
    }

    template <typename... Args>
    void header(Args&&... cols) { row(std::forward<Args>(cols)...); }

  private:
    std::ofstream ofs_;
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

static std::string schemeName(Scheme s) {
    switch (s) {
      case Scheme::StandardCentral:           return "StandardCentral";
      case Scheme::ExponentialFitting:        return "ExponentialFitting";
      case Scheme::MilevTaglianiCNEffectiveDiffusion: return "MilevTaglianiCN";
      default:                                return "Unknown";
    }
}

static FdmBlackScholesSpatialDesc descForScheme(
        Scheme s, Policy policy = Policy::None) {
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

    auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    const Array locs = mesher->locations(0);

    for (auto scheme : ALL_SCHEMES) {
        auto desc = descForScheme(scheme);
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

        std::string schName = schemeName(scheme);
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
            Policy::FallbackToExponentialFitting);

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

// Helper: build a BS FDM solver on a FdmBlackScholesMesher and return
// the price at the given spot.
static Real solveEuropean(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        const ext::shared_ptr<StrikedTypePayoff>& payoff,
        const FdmBlackScholesSpatialDesc& desc,
        Real K, Time T, Size xGrid, Size tGrid, Real spot) {
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

    return solver.valueAt(spot);
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

    Size grids[] = {25, 50, 100, 200, 400, 800, 1600};

    ext::shared_ptr<StrikedTypePayoff> payoff(
        new PlainVanillaPayoff(Option::Call, K));

    for (auto scheme : ALL_SCHEMES) {
        std::string schName = schemeName(scheme);
        CsvWriter csv(dataDir + "/grid_convergence_" + schName + ".csv");
        writeStandardMeta(csv, schName, schName,
                          Size(0), Size(0), r, q, vol, K, T,
                          "None", "FdmBlackScholesMesher");
        csv.meta("spot", spot);
        csv.meta("reference_price", refPrice);
        csv.header("xGrid", "tGrid", "price", "error");

        auto desc = descForScheme(scheme);
        for (Size xGrid : grids) {
            Size tGrid = 4 * xGrid;
            Real price = solveEuropean(process, payoff, desc,
                                       K, T, xGrid, tGrid, spot);
            csv.row(xGrid, tGrid,
                    price, std::fabs(price - refPrice));
        }
    }

    std::cout << " done.\n";
}

// Write sweep metadata block for volatility-sweep experiments.
static void writeSweepMeta(CsvWriter& csv,
                           const std::string& schemeName_,
                           const std::string& effectiveScheme,
                           Size xGrid, Size tGrid,
                           Rate r, Rate q,
                           Volatility sigmaMin, Volatility sigmaMax,
                           Size nSigma,
                           Real strike, Time maturity,
                           const std::string& mMatrixPol,
                           const std::string& mesh) {
    csv.meta("scheme", schemeName_);
    csv.meta("effective_scheme", effectiveScheme);
    csv.meta("xGrid", xGrid);
    csv.meta("tGrid", tGrid);
    csv.meta("r", r);
    csv.meta("q", q);
    csv.meta("sigma", "sweep");
    csv.meta("sigmaMin", sigmaMin);
    csv.meta("sigmaMax", sigmaMax);
    csv.meta("nSigma", nSigma);
    csv.meta("sweepSpacing", "log");
    csv.meta("strike", strike);
    csv.meta("maturity", maturity);
    csv.meta("mMatrixPolicy", mMatrixPol);
    csv.meta("mesh", mesh);
}

// =====================================================================
// Shared σ-sweep infrastructure for Experiments 6 & 7
// =====================================================================

static std::vector<Volatility> logSpacedSigmas(Volatility sigmaMin,
                                                Volatility sigmaMax,
                                                Size n) {
    std::vector<Volatility> v(n);
    for (Size i = 0; i < n; ++i) {
        Real t = Real(i) / Real(n - 1);
        v[i] = sigmaMin * std::pow(sigmaMax / sigmaMin, t);
    }
    return v;
}

using SweepRowWriter = std::function<void(CsvWriter&, Volatility,
                                          const SparseMatrix&, Size, Real)>;

static void runSigmaSweep(const std::string& dataDir,
                           const std::string& filePrefix,
                           const std::string& label,
                           const std::function<void(CsvWriter&, Real)>& writeMeta,
                           const SweepRowWriter& writeRow) {
    std::cout << "  Running " << label << "..." << std::flush;

    const Real spot = 100.0;
    const Rate r = 0.05;
    const Rate q = 0.0;
    const Real K = 100.0;
    const Time T = 1.0;

    const Real xMin = std::log(50.0);
    const Real xMax = std::log(200.0);
    const Size xGrid = 200;
    const Real h = (xMax - xMin) / (xGrid - 1);
    const Time dt = T / 50.0;
    const Size midNode = xGrid / 2;

    const Volatility sigmaMin = 0.001;
    const Volatility sigmaMax = 0.5;
    const Size nSigma = 50;

    auto sigmas = logSpacedSigmas(sigmaMin, sigmaMax, nSigma);

    for (auto scheme : ALL_SCHEMES) {
        std::string schName = schemeName(scheme);
        CsvWriter csv(dataDir + "/" + filePrefix + schName + ".csv");
        writeSweepMeta(csv, schName, schName,
                       xGrid, Size(0), r, q,
                       sigmaMin, sigmaMax, nSigma,
                       K, T, "None", "Uniform1dMesher");
        writeMeta(csv, h);

        for (Size j = 0; j < nSigma; ++j) {
            auto proc = makeProcess(spot, r, q, sigmas[j]);
            auto mesher = ext::make_shared<FdmMesherComposite>(
                ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));
            auto desc = descForScheme(scheme, Policy::None);
            auto op = ext::make_shared<FdmBlackScholesOp>(
                mesher, proc, K,
                false, -Null<Real>(), 0,
                ext::shared_ptr<FdmQuantoHelper>(), desc);
            op->setTime(0.0, dt);
            SparseMatrix mat = op->toMatrix();

            writeRow(csv, sigmas[j], mat, midNode, h);
        }
    }

    std::cout << " done.\n";
}

// =====================================================================
// Experiment 6: Effective diffusion σ-sweep (Fig 6)
// =====================================================================

void runEffectiveDiffusionSweep(const std::string& dataDir) {
    runSigmaSweep(
        dataDir, "effective_diffusion_sweep_",
        "effective diffusion sweep",
        [](CsvWriter& csv, Real h) {
            csv.meta("h", h);
            csv.header("sigma", "a_eff", "a_base", "ratio");
        },
        [](CsvWriter& csv, Volatility sigma,
           const SparseMatrix& mat, Size midNode, Real h) {
            Real lower = Real(mat(midNode, midNode - 1));
            Real upper = Real(mat(midNode, midNode + 1));
            Real aEff = (lower + upper) * h * h / 2.0;
            Real aBase = sigma * sigma / 2.0;
            csv.row(sigma, aEff, aBase,
                    aBase > 0.0 ? aEff / aBase : 0.0);
        });
}

// =====================================================================
// Experiment 7: M-matrix off-diagonal σ-sweep (Fig 7)
// =====================================================================

void runMMatrixSweep(const std::string& dataDir) {
    runSigmaSweep(
        dataDir, "mmatrix_sweep_",
        "M-matrix sweep",
        [](CsvWriter& csv, Real) {
            csv.header("sigma", "lower", "upper");
        },
        [](CsvWriter& csv, Volatility sigma,
           const SparseMatrix& mat, Size midNode, Real) {
            csv.row(sigma,
                    Real(mat(midNode, midNode - 1)),
                    Real(mat(midNode, midNode + 1)));
        });
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

    Size grids[] = {50, 100, 200, 400, 800, 1600};

    ext::shared_ptr<StrikedTypePayoff> payoff(
        new PlainVanillaPayoff(Option::Call, K));

    for (auto scheme : ALL_SCHEMES) {
        std::string schName = schemeName(scheme);

        CsvWriter csv(dataDir + "/benchmark_" + schName + ".csv");
        writeStandardMeta(csv, schName, schName,
                          Size(0), Size(0), r, q, vol, K, T,
                          "None", "FdmBlackScholesMesher");
        csv.meta("spot", spot);
        csv.meta("reference_price", refPrice);
        csv.meta("cost_note", "relative_cost = xGrid * tGrid (deterministic proxy for runtime)");
        csv.header("xGrid", "tGrid", "relative_cost", "price", "error");

        auto desc = descForScheme(scheme);
        for (Size xGrid : grids) {
            Size tGrid = 4 * xGrid;
            Real price = solveEuropean(process, payoff, desc,
                                       K, T, xGrid, tGrid, spot);
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

    for (auto scheme : ALL_SCHEMES) {
        auto desc = descForScheme(scheme);
        std::string schName = schemeName(scheme);

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

        // Fig 6: Effective diffusion σ-sweep
        runEffectiveDiffusionSweep(dataDir);

        // Fig 7: M-matrix off-diagonal σ-sweep
        runMMatrixSweep(dataDir);

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
