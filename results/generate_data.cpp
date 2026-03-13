/* generate_data.cpp — Numerical experiments for nonstandard FD scheme analysis.
 *
 * Links against the QuantLib build via CMake find_package(QuantLib).
 * Outputs CSV files to results/data/ for consumption by plot_figures.py.
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
#include <fstream>
#include <iomanip>
#include <iostream>
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

static FdmBlackScholesSpatialDesc descForScheme(
        FdmBlackScholesSpatialDesc::Scheme s,
        FdmBlackScholesSpatialDesc::MMatrixPolicy policy =
            FdmBlackScholesSpatialDesc::MMatrixPolicy::None) {
    FdmBlackScholesSpatialDesc d;
    d.scheme = s;
    d.mMatrixPolicy = policy;
    return d;
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
// Experiment 1-2: Truncated call (Figs 1-2)
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

    CsvWriter csv(dataDir + "/truncated_call.csv");
    csv.meta("experiment", "truncated_call");
    csv.meta("r", r);
    csv.meta("q", q);
    csv.meta("sigma", vol);
    csv.meta("strike", K);
    csv.meta("upper_barrier", U);
    csv.meta("maturity", T);
    csv.meta("xGrid", xGrid);
    csv.meta("tGrid", tGrid);
    csv.meta("mesh", "Uniform1dMesher");
    csv.meta("time_scheme", "CrankNicolson");
    csv.meta("mMatrixPolicy", "None");
    csv.header("S", "StandardCentral", "ExponentialFitting", "MilevTaglianiCN");

    // Solve for each scheme
    std::vector<Array> solutions(3);
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

        solutions[s].resize(xGrid);
        for (Size i = 0; i < xGrid; ++i) {
            solutions[s][i] = solver.valueAt(std::exp(locs[i]));
        }
    }

    // Write rows
    for (Size i = 0; i < xGrid; ++i) {
        Real S = std::exp(locs[i]);
        csv.row(S, solutions[0][i], solutions[1][i], solutions[2][i]);
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

    CsvWriter csv(dataDir + "/grid_convergence.csv");
    csv.meta("experiment", "grid_convergence");
    csv.meta("spot", spot);
    csv.meta("r", r);
    csv.meta("q", q);
    csv.meta("sigma", vol);
    csv.meta("strike", K);
    csv.meta("maturity", T);
    csv.meta("reference_price", refPrice);
    csv.meta("time_scheme", "CrankNicolson");
    csv.meta("mMatrixPolicy", "None");
    csv.header("xGrid", "tGrid",
               "StandardCentral", "ExponentialFitting", "MilevTaglianiCN",
               "err_StandardCentral", "err_ExponentialFitting",
               "err_MilevTaglianiCN");

    ext::shared_ptr<StrikedTypePayoff> payoff(
        new PlainVanillaPayoff(Option::Call, K));

    for (Size xGrid : grids) {
        Size tGrid = 4 * xGrid;
        Real prices[3];

        for (Size s = 0; s < 3; ++s) {
            auto desc = descForScheme(schemes[s]);
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

            prices[s] = solver.valueAt(spot);
        }

        csv.row(xGrid, tGrid,
                prices[0], prices[1], prices[2],
                std::fabs(prices[0] - refPrice),
                std::fabs(prices[1] - refPrice),
                std::fabs(prices[2] - refPrice));
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

    // Off-diagonal data (Fig 7)
    {
        CsvWriter csv(dataDir + "/mmatrix_offdiag.csv");
        csv.meta("experiment", "mmatrix_offdiag");
        csv.meta("r", r);
        csv.meta("q", q);
        csv.meta("sigma", vol);
        csv.meta("xGrid", xGrid);
        csv.meta("mesh", "Uniform1dMesher");
        csv.meta("mMatrixPolicy", "None");
        csv.header("node", "S",
                   "lower_SC", "upper_SC",
                   "lower_EF", "upper_EF",
                   "lower_MT", "upper_MT");

        SparseMatrix mats[3];
        for (Size s = 0; s < 3; ++s) {
            auto desc = descForScheme(schemes[s]);
            auto op = ext::make_shared<FdmBlackScholesOp>(
                mesher, process, K,
                false, -Null<Real>(), 0,
                ext::shared_ptr<FdmQuantoHelper>(), desc);
            op->setTime(0.0, dt);
            mats[s] = op->toMatrix();
        }

        for (Size i = 1; i + 1 < xGrid; ++i) {
            Real S = std::exp(xMin + i * h);
            csv.row(i, S,
                    Real(mats[0](i, i-1)), Real(mats[0](i, i+1)),
                    Real(mats[1](i, i-1)), Real(mats[1](i, i+1)),
                    Real(mats[2](i, i-1)), Real(mats[2](i, i+1)));
        }
    }

    // Effective diffusion data (Fig 6)
    {
        CsvWriter csv(dataDir + "/effective_diffusion.csv");
        csv.meta("experiment", "effective_diffusion");
        csv.meta("r", r);
        csv.meta("q", q);
        csv.meta("sigma", vol);
        csv.meta("xGrid", xGrid);
        csv.meta("h", h);
        csv.meta("mesh", "Uniform1dMesher");
        csv.meta("mMatrixPolicy", "None");
        csv.header("node", "S",
                   "aUsed_SC", "aUsed_EF", "aUsed_MT");

        SparseMatrix mats[3];
        for (Size s = 0; s < 3; ++s) {
            auto desc = descForScheme(schemes[s]);
            auto op = ext::make_shared<FdmBlackScholesOp>(
                mesher, process, K,
                false, -Null<Real>(), 0,
                ext::shared_ptr<FdmQuantoHelper>(), desc);
            op->setTime(0.0, dt);
            mats[s] = op->toMatrix();
        }

        for (Size i = 1; i + 1 < xGrid; ++i) {
            Real S = std::exp(xMin + i * h);
            Real aUsed[3];
            for (Size s = 0; s < 3; ++s) {
                Real lower = Real(mats[s](i, i-1));
                Real upper = Real(mats[s](i, i+1));
                aUsed[s] = (lower + upper) * h * h / 2.0;
            }
            csv.row(i, S, aUsed[0], aUsed[1], aUsed[2]);
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
    csv.meta("experiment", "xcothx");
    csv.meta("xSmall", 1e-6);
    csv.meta("xLarge", 50.0);
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
    const Size nRuns = 5;

    CsvWriter csv(dataDir + "/benchmark.csv");
    csv.meta("experiment", "benchmark");
    csv.meta("spot", spot);
    csv.meta("r", r);
    csv.meta("q", q);
    csv.meta("sigma", vol);
    csv.meta("strike", K);
    csv.meta("maturity", T);
    csv.meta("reference_price", refPrice);
    csv.meta("time_scheme", "CrankNicolson");
    csv.meta("mMatrixPolicy", "FallbackToExponentialFitting");
    csv.meta("num_runs", nRuns);
    csv.header("xGrid", "scheme", "median_time_ms", "price", "error");

    ext::shared_ptr<StrikedTypePayoff> payoff(
        new PlainVanillaPayoff(Option::Call, K));

    for (Size xGrid : grids) {
        Size tGrid = 4 * xGrid;

        for (Size s = 0; s < 3; ++s) {
            auto desc = descForScheme(schemes[s],
                FdmBlackScholesSpatialDesc::MMatrixPolicy::FallbackToExponentialFitting);

            std::vector<double> timings;
            Real price = 0.0;

            for (Size run = 0; run < nRuns; ++run) {
                auto t0 = std::chrono::high_resolution_clock::now();

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

                price = solver.valueAt(spot);

                auto t1 = std::chrono::high_resolution_clock::now();
                double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
                timings.push_back(ms);
            }

            std::sort(timings.begin(), timings.end());
            double medianTime = timings[nRuns / 2];

            csv.row(xGrid, schemeName(schemes[s]),
                    medianTime, price,
                    std::fabs(price - refPrice));
        }
    }

    std::cout << " done.\n";
}

// =====================================================================
// Experiment 3-4: Discrete double barrier knock-out (Figs 3-4)
// =====================================================================

void runDiscreteBarrier(const std::string& dataDir,
                        const std::string& label,
                        Rate r, Rate q, Volatility vol,
                        Real K, Real L, Real U, Time T,
                        Size nMonitoring, Size xGrid, Size tGrid) {
    std::cout << "  Running discrete barrier (" << label << ")..." << std::flush;

    auto process = makeProcess(100.0, r, q, vol);

    using Scheme = FdmBlackScholesSpatialDesc::Scheme;
    Scheme schemes[] = {Scheme::StandardCentral,
                        Scheme::ExponentialFitting,
                        Scheme::MilevTaglianiCNEffectiveDiffusion};

    // Monitoring dates (equally spaced)
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

    CsvWriter csv(dataDir + "/barrier_" + label + ".csv");
    csv.meta("experiment", "discrete_barrier_" + label);
    csv.meta("r", r);
    csv.meta("q", q);
    csv.meta("sigma", vol);
    csv.meta("strike", K);
    csv.meta("lower_barrier", L);
    csv.meta("upper_barrier", U);
    csv.meta("maturity", T);
    csv.meta("n_monitoring", nMonitoring);
    csv.meta("xGrid", xGrid);
    csv.meta("tGrid", tGrid);
    csv.meta("mesh", "FdmBlackScholesMesher");
    csv.meta("time_scheme", "CrankNicolson");
    csv.meta("mMatrixPolicy", "None");
    csv.header("S", "StandardCentral", "ExponentialFitting", "MilevTaglianiCN");

    // For each scheme, solve with barrier step condition
    // Using manual rollback with FdmBackwardSolver and step conditions
    ext::shared_ptr<StrikedTypePayoff> payoff(
        new PlainVanillaPayoff(Option::Call, K));

    std::vector<std::vector<Real>> results(3);
    for (Size s = 0; s < 3; ++s) {
        auto desc = descForScheme(schemes[s]);

        std::vector<std::tuple<Real, Real, bool>> cPoints = {
            {K, 0.1, true}, {L, 0.1, true}, {U, 0.1, true}
        };
        auto equityMesher = ext::make_shared<FdmBlackScholesMesher>(
            xGrid, process, T, K,
            Null<Real>(), Null<Real>(), 0.0001, 1.5,
            cPoints);
        auto mesher = ext::make_shared<FdmMesherComposite>(equityMesher);
        auto calc = ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);

        // Set up the terminal payoff
        const Size n = mesher->layout()->size();
        Array rhs(n);
        for (const auto& iter : *mesher->layout()) {
            Real S = std::exp(mesher->location(iter, 0));
            rhs[iter.index()] = (*payoff)(S);
            // Apply initial barrier condition
            if (S < L || S > U)
                rhs[iter.index()] = 0.0;
        }

        // Create operator
        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), desc);

        // Roll back with barrier application at monitoring dates
        ext::shared_ptr<FdmStepConditionComposite> noCondition;
        FdmBackwardSolver solver(op, FdmBoundaryConditionSet(),
                                 noCondition,
                                 FdmSchemeDesc::CrankNicolson());

        // Roll back from T to 0, applying barrier at monitoring dates
        Real tPrev = T;
        for (auto it = monitoringTimes.rbegin();
             it != monitoringTimes.rend(); ++it) {
            Real tMon = *it;
            Size nSteps = std::max(Size(1),
                Size(std::round(tGrid * (tPrev - tMon) / T)));
            solver.rollback(rhs, tPrev, tMon, nSteps, 0);

            // Apply barrier condition: zero outside [L, U]
            for (const auto& iter : *mesher->layout()) {
                Real S = std::exp(mesher->location(iter, 0));
                if (S < L || S > U)
                    rhs[iter.index()] = 0.0;
            }
            tPrev = tMon;
        }
        // Final rollback to t=0
        if (tPrev > 0.0) {
            Size nSteps = std::max(Size(1),
                Size(std::round(tGrid * tPrev / T)));
            solver.rollback(rhs, tPrev, 0.0, nSteps, 0);
        }

        // Interpolate at requested spot values
        const Array barrierLocs = mesher->locations(0);
        results[s].resize(spots.size());
        for (Size j = 0; j < spots.size(); ++j) {
            Real x = std::log(spots[j]);
            if (x <= barrierLocs[0]) {
                results[s][j] = rhs[0];
            } else if (x >= barrierLocs[n-1]) {
                results[s][j] = rhs[n-1];
            } else {
                // Binary search for bracket
                Size lo = 0, hi = n - 1;
                while (hi - lo > 1) {
                    Size mid = (lo + hi) / 2;
                    if (barrierLocs[mid] <= x)
                        lo = mid;
                    else
                        hi = mid;
                }
                Real w = (x - barrierLocs[lo]) / (barrierLocs[hi] - barrierLocs[lo]);
                results[s][j] = (1.0 - w) * rhs[lo] + w * rhs[hi];
            }
        }
    }

    for (Size j = 0; j < spots.size(); ++j) {
        csv.row(spots[j], results[0][j], results[1][j], results[2][j]);
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
    csv.meta("experiment", "mc_barrier_" + label);
    csv.meta("r", r);
    csv.meta("q", q);
    csv.meta("sigma", vol);
    csv.meta("strike", K);
    csv.meta("lower_barrier", L);
    csv.meta("upper_barrier", U);
    csv.meta("maturity", T);
    csv.meta("n_monitoring", nMonitoring);
    csv.meta("num_paths", nPaths);
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

        csv.row(S0, mean, se, nPaths);
    }

    std::cout << " done.\n";
}

// =====================================================================
// Main
// =====================================================================

int main() {
    try {
        Date today(28, March, 2004);
        Settings::instance().evaluationDate() = today;

        std::string dataDir = "data";
        // Create data directory if running from results/build
        // Adjust path relative to the source results/ directory
        {
            std::ifstream test(dataDir + "/.gitkeep");
            if (!test) {
                dataDir = "../data";
                std::ifstream test2(dataDir + "/.gitkeep");
                if (!test2) {
                    // Try creating from current dir
                    dataDir = "data";
                    system("mkdir -p data");
                }
            }
        }

        std::cout << "Generating results data...\n";

        // Fig 9: xCothx/Peclet (fastest, no solver)
        runXCothx(dataDir);

        // Figs 6-7: Operator diagnostics
        runOperatorDiagnostics(dataDir);

        // Fig 5: Grid convergence
        runGridConvergence(dataDir);

        // Figs 1-2: Truncated call
        runTruncatedCall(dataDir);

        // Fig 8: Performance benchmark
        runBenchmark(dataDir);

        // Fig 3: Discrete barrier moderate-vol
        runDiscreteBarrier(dataDir, "moderate_vol",
                           0.05, 0.0, 0.25,
                           100.0, 95.0, 110.0, 0.5,
                           5, 2000, 50000);

        // Fig 4: Discrete barrier low-vol
        runDiscreteBarrier(dataDir, "low_vol",
                           0.05, 0.0, 0.001,
                           100.0, 95.0, 110.0, 1.0,
                           5, 2000, 50000);

        // MC reference: moderate-vol barrier
        runBarrierMC(dataDir, "moderate_vol",
                     0.05, 0.0, 0.25,
                     100.0, 95.0, 110.0, 0.5,
                     5, 1000000);

        // MC reference: low-vol barrier
        runBarrierMC(dataDir, "low_vol",
                     0.05, 0.0, 0.001,
                     100.0, 95.0, 110.0, 1.0,
                     5, 1000000);

        std::cout << "All experiments complete. CSV files in " << dataDir << "/\n";
        return 0;
    }
    catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
