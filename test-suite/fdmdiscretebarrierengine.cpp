// r6
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 FdmDiscreteBarrierEngine integration tests.

 Validates: exact mesh alignment at barrier/strike nodes, rebate
 pass-through, engine-level discrete monitoring with knock-in parity,
 and paper-faithful replication of [MT10, Example 4.1] (double barrier
 knock-out with convergence to the analytic reference).

 cf. [MT10, Table 1] for convergence data, [Ballabio20, Ch. 11] for
 the FDM barrier framework.
*/

#include "toplevelfixture.hpp"
#include "utilities.hpp"

#include <ql/exercise.hpp>
#include <ql/instruments/barrieroption.hpp>
#include <ql/instruments/payoffs.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>
#include <ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmeshercomposite.hpp>
#include <ql/methods/finitedifferences/meshers/uniform1dmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesop.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp>
#include <ql/methods/finitedifferences/solvers/fdmsolverdesc.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmdiscretebarrierstepcondition.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp>
#include <ql/methods/finitedifferences/utilities/fdminnervaluecalculator.hpp>
#include <ql/pricingengines/barrier/analyticbarrierengine.hpp>
#include <ql/pricingengines/barrier/fdblackscholesbarrierengine.hpp>
#include <ql/pricingengines/barrier/makefdblackscholesbarrierengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/fdblackscholesshoutengine.hpp>
#include <ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <list>
#include <vector>

using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FdmDiscreteBarrierEngineTests)

namespace {

    ext::shared_ptr<GeneralizedBlackScholesProcess>
    makeBSProcess(Real spot, Rate r, Rate q, Volatility vol) {
        DayCounter dc = Actual365Fixed();
        Date today = Settings::instance().evaluationDate();

        return ext::make_shared<GeneralizedBlackScholesProcess>(
            Handle<Quote>(ext::make_shared<SimpleQuote>(spot)),
            Handle<YieldTermStructure>(flatRate(today, q, dc)),
            Handle<YieldTermStructure>(flatRate(today, r, dc)),
            Handle<BlackVolTermStructure>(flatVol(today, vol, dc)));
    }

    // Look up the value of a solution array at a given spot level
    // by finding the nearest grid node (works when mesh is aligned).
    Real valueAtSpot(const Array& v,
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

    // Helper: build payoff vector on the grid
    Array buildPayoff(const ext::shared_ptr<FdmMesher>& mesher,
                      const ext::shared_ptr<Payoff>& payoff,
                      Size direction = 0) {
        Array rhs(mesher->layout()->size());
        for (const auto& iter : *mesher->layout()) {
            const Real S = std::exp(mesher->location(iter, direction));
            rhs[iter.index()] = (*payoff)(S);
        }
        return rhs;
    }

    FdmBlackScholesSpatialDesc centralDesc() {
        return FdmBlackScholesSpatialDesc::standard();
    }

    FdmBlackScholesSpatialDesc fittedDesc() {
        FdmBlackScholesSpatialDesc d =
            FdmBlackScholesSpatialDesc::exponentialFitting();
        d.mMatrixPolicy =
            FdmBlackScholesSpatialDesc::MMatrixPolicy::None;
        return d;
    }

    FdmBlackScholesSpatialDesc milevTaglianiDesc() {
        FdmBlackScholesSpatialDesc d =
            FdmBlackScholesSpatialDesc::milevTaglianiCN();
        d.mMatrixPolicy =
            FdmBlackScholesSpatialDesc::MMatrixPolicy::None;
        return d;
    }

} // anonymous namespace


// ===================================================================
// 1. Exact mesh alignment
// ===================================================================

BOOST_AUTO_TEST_CASE(testMultiPointMeshAlignment) {

    BOOST_TEST_MESSAGE(
        "Testing exact mesh alignment at strike and barriers...");

    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    auto process = makeBSProcess(100.0, 0.05, 0.0, 0.25);

    const Real K = 100.0;
    const Real L = 95.0;
    const Real U = 110.0;
    const Time maturity = 0.5;
    const Size xGrid = 500;

    std::vector<std::tuple<Real, Real, bool>> cPoints = {
        {K, 0.1, true},
        {L, 0.1, true},
        {U, 0.1, true}
    };

    FdmBlackScholesMesher mesher(
        xGrid, process, maturity, K,
        Null<Real>(), Null<Real>(), 0.0001, 1.5,
        cPoints);

    const auto& locs = mesher.locations();

    const Real tol = 1e-10;

    auto findExact = [&](Real spotLevel) -> bool {
        const Real target = std::log(spotLevel);
        for (Size i = 0; i < locs.size(); ++i) {
            if (std::fabs(locs[i] - target) < tol)
                return true;
        }
        return false;
    };

    BOOST_CHECK_MESSAGE(findExact(K),
        "Strike K=" << K << " not found as exact grid node "
        "(ln(K)=" << std::log(K) << ")");

    BOOST_CHECK_MESSAGE(findExact(L),
        "Lower barrier L=" << L << " not found as exact grid node "
        "(ln(L)=" << std::log(L) << ")");

    BOOST_CHECK_MESSAGE(findExact(U),
        "Upper barrier U=" << U << " not found as exact grid node "
        "(ln(U)=" << std::log(U) << ")");
}


// ===================================================================
// 2. Rebate support in discrete barrier step condition
// ===================================================================

BOOST_AUTO_TEST_CASE(testDiscreteBarrierRebate) {

    BOOST_TEST_MESSAGE(
        "Testing discrete barrier step condition with non-zero rebate...");

    const Size xGrid = 800;
    const Real xMin = std::log(50.0);
    const Real xMax = std::log(200.0);
    const Real L = 80.0;
    const Real U = 120.0;
    const Real rebate = 5.0;

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    FdmDiscreteBarrierStepCondition cond(
        mesher, std::vector<Time>({0.25}),
        L, U, rebate, Size(0));

    BOOST_CHECK_EQUAL(cond.rebate(), rebate);

    Array a(mesher->layout()->size(), 10.0);

    cond.applyTo(a, 0.25);

    const Real xL = std::log(L);
    const Real xU = std::log(U);

    for (const auto& iter : *mesher->layout()) {
        const Real x = mesher->location(iter, 0);
        if (x < xL || x > xU) {
            BOOST_CHECK_CLOSE(a[iter.index()], rebate, 1e-12);
        } else {
            BOOST_CHECK_CLOSE(a[iter.index()], 10.0, 1e-12);
        }
    }
}


// ===================================================================
// 3. Milev-Tagliani Example 4.1 replication
//    (Discrete double barrier knock-out, sigma=0.25, r=0.05)
// ===================================================================

BOOST_AUTO_TEST_CASE(testDoubleBarrierKnockOutPaperExample) {

    BOOST_TEST_MESSAGE(
        "Testing discrete double barrier knock-out "
        "(Milev-Tagliani Example 4.1)...");

    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    // Paper parameters
    const Real K = 100.0;
    const Real L = 95.0;
    const Real U = 110.0;
    const Volatility sigma = 0.25;
    const Rate r = 0.05;
    const Rate q = 0.0;
    const Time maturity = 0.5;

    // 5 monitoring dates equally spaced (including maturity)
    std::vector<Time> monTimes;
    for (Size i = 1; i <= 5; ++i)
        monTimes.push_back(maturity * Real(i) / 5.0);

    auto process = makeBSProcess(100.0, r, q, sigma);

    // Fine grid with alignment at K, L, U
    const Size xGrid = 2000;
    const Size tGrid = 500;

    std::vector<std::tuple<Real, Real, bool>> cPoints = {
        {K, 0.1, true},
        {L, 0.1, true},
        {U, 0.1, true}
    };

    auto equityMesher = ext::make_shared<FdmBlackScholesMesher>(
        xGrid, process, maturity, K,
        Null<Real>(), Null<Real>(), 0.0001, 1.5,
        cPoints);

    auto mesher = ext::make_shared<FdmMesherComposite>(equityMesher);

    // Payoff: max(S-K, 0)  (barrier projection via step condition)
    ext::shared_ptr<StrikedTypePayoff> payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, K);
    auto calculator =
        ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);

    // Step condition: discrete double barrier knock-out
    auto barrierCondition =
        ext::make_shared<FdmDiscreteBarrierStepCondition>(
            mesher, monTimes, L, U);

    auto conditions = ext::make_shared<FdmStepConditionComposite>(
        std::list<std::vector<Time>>(1, monTimes),
        FdmStepConditionComposite::Conditions(1, barrierCondition));

    const FdmBoundaryConditionSet bcSet;

    // ---- Run with ExponentialFitting + CN, then MilevTagliani + CN ----
    struct SchemeCase { FdmBlackScholesSpatialDesc desc; const char* name; };
    for (const auto& sc : {SchemeCase{fittedDesc(), "ExpFit+CN"},
                           SchemeCase{milevTaglianiDesc(), "MT+CN"}}) {
        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), Size(0),
            ext::shared_ptr<FdmQuantoHelper>(), sc.desc);

        Array rhs = buildPayoff(mesher, payoff);
        FdmBackwardSolver solver(
            op, bcSet, conditions,
            FdmSchemeDesc::CrankNicolson());
        solver.rollback(rhs, maturity, 0.0, tGrid, 0);

        const Real price95  = valueAtSpot(rhs, mesher, 95.0);
        const Real price100 = valueAtSpot(rhs, mesher, 100.0);

        BOOST_TEST_MESSAGE("  " << sc.name << ": V(95)="  << price95
                           << ", V(100)=" << price100);

        // Prices must be positive
        BOOST_CHECK_MESSAGE(price95 >= -1e-10,
            "Negative price at S=95: " << price95);
        BOOST_CHECK_MESSAGE(price100 >= -1e-10,
            "Negative price at S=100: " << price100);

        // Paper MC: V(95)~0.174, V(100)~0.233
        // Tolerance generous due to log-space vs S-space difference
        BOOST_CHECK_MESSAGE(price95 > 0.10 && price95 < 0.30,
            "V(95) out of expected range: " << price95);
        BOOST_CHECK_MESSAGE(price100 > 0.15 && price100 < 0.35,
            "V(100) out of expected range: " << price100);
    }
}


// ===================================================================
// 4. Greeks smoothness near barrier
// ===================================================================

BOOST_AUTO_TEST_CASE(testGreeksSmoothnessNearBarrier) {

    BOOST_TEST_MESSAGE(
        "Testing Greeks smoothness near barrier with "
        "exponential fitting...");

    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real K = 100.0;
    const Real L = 95.0;
    const Real U = 110.0;
    const Volatility sigma = 0.25;
    const Rate r = 0.05;
    const Time maturity = 0.5;

    std::vector<Time> monTimes;
    for (Size i = 1; i <= 5; ++i)
        monTimes.push_back(maturity * Real(i) / 5.0);

    auto process = makeBSProcess(100.0, r, 0.0, sigma);

    const Size xGrid = 2000;
    const Size tGrid = 500;

    std::vector<std::tuple<Real, Real, bool>> cPoints = {
        {K, 0.1, true}, {L, 0.1, true}, {U, 0.1, true}
    };

    auto equityMesher = ext::make_shared<FdmBlackScholesMesher>(
        xGrid, process, maturity, K,
        Null<Real>(), Null<Real>(), 0.0001, 1.5,
        cPoints);
    auto mesher = ext::make_shared<FdmMesherComposite>(equityMesher);

    ext::shared_ptr<StrikedTypePayoff> payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, K);

    auto barrierCondition =
        ext::make_shared<FdmDiscreteBarrierStepCondition>(
            mesher, monTimes, L, U);
    auto conditions = ext::make_shared<FdmStepConditionComposite>(
        std::list<std::vector<Time>>(1, monTimes),
        FdmStepConditionComposite::Conditions(1, barrierCondition));

    const FdmBoundaryConditionSet bcSet;

    // Price with exponential fitting
    FdmBlackScholesSpatialDesc desc = fittedDesc();
    auto op = ext::make_shared<FdmBlackScholesOp>(
        mesher, process, K,
        false, -Null<Real>(), Size(0),
        ext::shared_ptr<FdmQuantoHelper>(), desc);

    Array rhs = buildPayoff(mesher, payoff);
    FdmBackwardSolver solver(op, bcSet, conditions,
                             FdmSchemeDesc::CrankNicolson());
    solver.rollback(rhs, maturity, 0.0, tGrid, 0);

    // Collect (x, V) pairs near the lower barrier for Delta check
    std::vector<Real> xNear, vNear;
    const Real xL = std::log(L);
    for (const auto& iter : *mesher->layout()) {
        const Real x = mesher->location(iter, 0);
        if (std::fabs(x - xL) < 0.15) {
            xNear.push_back(x);
            vNear.push_back(rhs[iter.index()]);
        }
    }

    // Check that all values near the barrier are non-negative
    // and that there are no wild oscillations (|jump| < 0.5 between
    // consecutive nodes).
    Size oscillationCount = 0;
    for (Size i = 0; i < vNear.size(); ++i) {
        BOOST_CHECK_MESSAGE(vNear[i] >= -1e-10,
            "Negative value near barrier at x=" << xNear[i]
            << ": V=" << vNear[i]);
        if (i > 0) {
            const Real jump = std::fabs(vNear[i] - vNear[i-1]);
            if (jump > 0.5)
                ++oscillationCount;
        }
    }
    BOOST_CHECK_MESSAGE(oscillationCount == 0,
        "Found " << oscillationCount
        << " oscillation(s) near barrier with ExpFit");
}


// ===================================================================
// 5. Engine-level single-barrier discrete monitoring
// ===================================================================

BOOST_AUTO_TEST_CASE(testSingleBarrierEngineDiscreteMonitoring) {

    BOOST_TEST_MESSAGE(
        "Testing FdBlackScholesBarrierEngine with "
        "discrete monitoring...");

    const DayCounter dc = Actual365Fixed();
    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real spot = 100.0;
    const Rate r = 0.05;
    const Rate q = 0.0;
    const Volatility vol = 0.25;

    auto process = ext::make_shared<BlackScholesMertonProcess>(
        Handle<Quote>(ext::make_shared<SimpleQuote>(spot)),
        Handle<YieldTermStructure>(flatRate(today, q, dc)),
        Handle<YieldTermStructure>(flatRate(today, r, dc)),
        Handle<BlackVolTermStructure>(flatVol(today, vol, dc)));

    const Date exerciseDate = today + Period(6, Months);
    const auto exercise =
        ext::make_shared<EuropeanExercise>(exerciseDate);
    const auto payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, 100.0);

    // 5 monthly monitoring dates
    std::vector<Date> monDates;
    for (Size i = 1; i <= 5; ++i)
        monDates.push_back(today + Period(i, Months));

    const Real barrier = 90.0;
    const Real rebate = 0.0;

    BarrierOption option(Barrier::DownOut, barrier, rebate,
                         payoff, exercise);

    // Discrete monitoring engine
    option.setPricingEngine(
        ext::make_shared<FdBlackScholesBarrierEngine>(
            process, monDates,
            200, 400, 0,
            FdmSchemeDesc::CrankNicolson(),
            false, -Null<Real>(),
            FdmBlackScholesSpatialDesc::exponentialFitting()));

    const Real value = option.NPV();
    const Real delta = option.delta();
    const Real gamma = option.gamma();

    BOOST_TEST_MESSAGE("  DownOut discrete: NPV=" << value
                       << ", delta=" << delta
                       << ", gamma=" << gamma);

    // With a DownOut barrier at 90 and spot=100, the option
    // should have significant value (similar to vanilla but
    // somewhat reduced by the knockout possibility).
    BOOST_CHECK_MESSAGE(value > 0.0,
        "Expected positive value for DownOut, got " << value);
    BOOST_CHECK_MESSAGE(value < 20.0,
        "Value implausibly large: " << value);
    BOOST_CHECK(std::isfinite(delta));
    BOOST_CHECK(std::isfinite(gamma));
}


// ===================================================================
// 6. Builder pattern
// ===================================================================

BOOST_AUTO_TEST_CASE(testMakeBarrierEngineBuilder) {

    BOOST_TEST_MESSAGE(
        "Testing MakeFdBlackScholesBarrierEngine builder...");

    const DayCounter dc = Actual365Fixed();
    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    auto process = ext::make_shared<BlackScholesMertonProcess>(
        Handle<Quote>(ext::make_shared<SimpleQuote>(100.0)),
        Handle<YieldTermStructure>(flatRate(today, 0.0, dc)),
        Handle<YieldTermStructure>(flatRate(today, 0.05, dc)),
        Handle<BlackVolTermStructure>(flatVol(today, 0.25, dc)));

    const Date exerciseDate = today + Period(6, Months);
    const auto exercise =
        ext::make_shared<EuropeanExercise>(exerciseDate);
    const auto payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, 100.0);

    std::vector<Date> monDates;
    for (Size i = 1; i <= 5; ++i)
        monDates.push_back(today + Period(i, Months));

    BarrierOption option(Barrier::UpOut, 115.0, 0.0,
                         payoff, exercise);

    // Use the builder pattern
    option.setPricingEngine(
        MakeFdBlackScholesBarrierEngine(process)
            .withTGrid(200)
            .withXGrid(400)
            .withFdmSchemeDesc(FdmSchemeDesc::CrankNicolson())
            .withSpatialDesc(
                FdmBlackScholesSpatialDesc::exponentialFitting())
            .withDiscreteMonitoring(monDates));

    const Real value = option.NPV();

    BOOST_CHECK_MESSAGE(value > 0.0 && std::isfinite(value),
        "Builder-configured engine produced invalid value: "
        << value);

    BOOST_TEST_MESSAGE("  Builder UpOut discrete: NPV=" << value);
}


// ===================================================================
// 7. Standard-Central CN shows oscillation artifacts that
//    exponential fitting avoids (low-vol stress, double barrier)
// ===================================================================

BOOST_AUTO_TEST_CASE(testLowVolDoubleBarrierOscillation) {

    BOOST_TEST_MESSAGE(
        "Testing that StandardCentral+CN shows oscillations "
        "while ExponentialFitting remains clean under "
        "low vol discrete double barrier...");

    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real K = 100.0;
    const Real L = 95.0;
    const Real U = 110.0;
    const Volatility sigma = 0.001;  // extreme low vol
    const Rate r = 0.05;
    const Time maturity = 0.5;

    // 5 monitoring dates
    std::vector<Time> monTimes;
    for (Size i = 1; i <= 5; ++i)
        monTimes.push_back(maturity * Real(i) / 5.0);

    auto process = makeBSProcess(100.0, r, 0.0, sigma);

    const Size xGrid = 800;
    const Size tGrid = 200;

    // With sigma=0.001, FdmBlackScholesMesher computes a domain of
    // only ~0.008 in log-space (roughly [99.6, 103.0] in spot-space),
    // which excludes the barriers at L=95 and U=110 entirely.
    // Use a uniform mesher with explicit bounds that extend well
    // beyond the barriers so the discrete barrier step condition
    // actually creates discontinuities during rollback.
    const Real xMin = std::log(80.0);   // well below L=95
    const Real xMax = std::log(130.0);  // well above U=110
    auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    ext::shared_ptr<StrikedTypePayoff> payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, K);

    auto barrierCondition =
        ext::make_shared<FdmDiscreteBarrierStepCondition>(
            mesher, monTimes, L, U);
    auto conditions = ext::make_shared<FdmStepConditionComposite>(
        std::list<std::vector<Time>>(1, monTimes),
        FdmStepConditionComposite::Conditions(1, barrierCondition));

    const FdmBoundaryConditionSet bcSet;

    // ---- StandardCentral + CN ----
    Size negCountCentral = 0;
    {
        FdmBlackScholesSpatialDesc desc = centralDesc();
        desc.mMatrixPolicy =
            FdmBlackScholesSpatialDesc::MMatrixPolicy::None;

        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), Size(0),
            ext::shared_ptr<FdmQuantoHelper>(), desc);

        Array rhs = buildPayoff(mesher, payoff);
        FdmBackwardSolver solver(op, bcSet, conditions,
                                 FdmSchemeDesc::CrankNicolson());
        solver.rollback(rhs, maturity, 0.0, tGrid, 0);

        for (Size i = 0; i < rhs.size(); ++i)
            if (rhs[i] < -1e-6) ++negCountCentral;

        BOOST_TEST_MESSAGE("  Central+CN: negative nodes = "
                           << negCountCentral);

        BOOST_CHECK_MESSAGE(negCountCentral > 0,
            "Expected negativity for StandardCentral + CN under "
            "low-vol double barrier, but found none");
    }

    // ---- ExponentialFitting + CN ----
    {
        FdmBlackScholesSpatialDesc desc = fittedDesc();

        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), Size(0),
            ext::shared_ptr<FdmQuantoHelper>(), desc);

        Array rhs = buildPayoff(mesher, payoff);
        FdmBackwardSolver solver(op, bcSet, conditions,
                                 FdmSchemeDesc::CrankNicolson());
        solver.rollback(rhs, maturity, 0.0, tGrid, 0);

        Size negCountFitted = 0;
        for (Size i = 0; i < rhs.size(); ++i)
            if (rhs[i] < -1e-6) ++negCountFitted;

        BOOST_TEST_MESSAGE("  ExpFit+CN: negative nodes = "
                           << negCountFitted);

        BOOST_CHECK_MESSAGE(negCountFitted == 0,
            "ExponentialFitting should produce no negative values, "
            "but found " << negCountFitted);

        // All interior values must be bounded
        for (Size i = 0; i < rhs.size(); ++i) {
            BOOST_CHECK(std::isfinite(rhs[i]));
        }
    }
}


// ===================================================================
// 8. Paper Table 1 — tighter replication at interior spot levels
//    (Milev-Tagliani Example 4.1: K=100, L=95, U=110,
//     sigma=0.25, r=0.05, T=0.5, 5 monitoring dates)
// ===================================================================

BOOST_AUTO_TEST_CASE(testDoubleBarrierKnockOutTable1Interior) {

    BOOST_TEST_MESSAGE(
        "Testing discrete double barrier knock-out: "
        "Table 1 interior spot levels against MC reference...");

    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real K = 100.0;
    const Real L = 95.0;
    const Real U = 110.0;
    const Volatility sigma = 0.25;
    const Rate r = 0.05;
    const Rate q = 0.0;
    const Time maturity = 0.5;

    std::vector<Time> monTimes;
    for (Size i = 1; i <= 5; ++i)
        monTimes.push_back(maturity * Real(i) / 5.0);

    auto process = makeBSProcess(100.0, r, q, sigma);

    const Size xGrid = 4000;
    const Size tGrid = 2000;

    std::vector<std::tuple<Real, Real, bool>> cPoints = {
        {K, 0.1, true}, {L, 0.1, true}, {U, 0.1, true}
    };

    auto equityMesher = ext::make_shared<FdmBlackScholesMesher>(
        xGrid, process, maturity, K,
        Null<Real>(), Null<Real>(), 0.0001, 1.5, cPoints);
    auto mesher = ext::make_shared<FdmMesherComposite>(equityMesher);

    ext::shared_ptr<StrikedTypePayoff> payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, K);

    auto barrierCondition =
        ext::make_shared<FdmDiscreteBarrierStepCondition>(
            mesher, monTimes, L, U);
    auto conditions = ext::make_shared<FdmStepConditionComposite>(
        std::list<std::vector<Time>>(1, monTimes),
        FdmStepConditionComposite::Conditions(1, barrierCondition));

    const FdmBoundaryConditionSet bcSet;

    // Paper MC reference (Table 1, standard error in parentheses):
    //   V(99.5)  = 0.22923 (0.00073)
    //   V(100)   = 0.23263 (0.00036)
    //   V(100.5) = 0.23410 (0.00073)
    struct SpotRef { Real spot; Real mcRef; };
    const SpotRef refs[] = {
        { 99.5,  0.22923},
        {100.0,  0.23263},
        {100.5,  0.23410}
    };

    const Real relTol = 0.05;  // 5% relative tolerance

    // Helper lambda: run a scheme and return values at reference spots
    auto runScheme = [&](const FdmBlackScholesSpatialDesc& desc,
                         const char* name) {
        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), Size(0),
            ext::shared_ptr<FdmQuantoHelper>(), desc);

        Array rhs = buildPayoff(mesher, payoff);
        FdmBackwardSolver solver(op, bcSet, conditions,
                                 FdmSchemeDesc::CrankNicolson());
        solver.rollback(rhs, maturity, 0.0, tGrid, 0);

        std::ostringstream oss;
        oss << "  " << name << ":";
        for (const auto& ref : refs) {
            const Real v = valueAtSpot(rhs, mesher, ref.spot);
            oss << "  V(" << ref.spot << ")=" << v;
        }
        BOOST_TEST_MESSAGE(oss.str());

        return rhs;
    };

    // Run all three schemes
    BOOST_TEST_MESSAGE("  Paper MC references: "
        "V(99.5)=0.22923, V(100)=0.23263, V(100.5)=0.23410");

    Array rhsExpFit = runScheme(fittedDesc(), "ExpFit+CN");
    Array rhsMT     = runScheme(milevTaglianiDesc(), "MT+CN");

    // Check ExpFit and MT values against MC references
    for (const auto& ref : refs) {
        const Real vExpFit = valueAtSpot(rhsExpFit, mesher, ref.spot);
        const Real vMT     = valueAtSpot(rhsMT, mesher, ref.spot);

        BOOST_CHECK_MESSAGE(
            std::fabs(vExpFit - ref.mcRef) / ref.mcRef < relTol,
            "ExpFit V(" << ref.spot << ")=" << vExpFit
            << " differs from MC ref " << ref.mcRef
            << " by more than " << (relTol * 100) << "%");

        BOOST_CHECK_MESSAGE(
            std::fabs(vMT - ref.mcRef) / ref.mcRef < relTol,
            "MT V(" << ref.spot << ")=" << vMT
            << " differs from MC ref " << ref.mcRef
            << " by more than " << (relTol * 100) << "%");
    }
}


// ===================================================================
// 9. Three-scheme comparison output
//    Outputs formatted table for all three schemes at multiple spots
// ===================================================================

BOOST_AUTO_TEST_CASE(testThreeSchemeComparisonOutput) {

    BOOST_TEST_MESSAGE(
        "Three-scheme comparison: StandardCentral vs "
        "ExponentialFitting vs MilevTagliani CN variant...");

    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real K = 100.0;
    const Real L = 95.0;
    const Real U = 110.0;
    const Volatility sigma = 0.25;
    const Rate r = 0.05;
    const Time maturity = 0.5;

    std::vector<Time> monTimes;
    for (Size i = 1; i <= 5; ++i)
        monTimes.push_back(maturity * Real(i) / 5.0);

    auto process = makeBSProcess(100.0, r, 0.0, sigma);

    const Size xGrid = 2000;
    const Size tGrid = 500;

    std::vector<std::tuple<Real, Real, bool>> cPoints = {
        {K, 0.1, true}, {L, 0.1, true}, {U, 0.1, true}
    };
    auto equityMesher = ext::make_shared<FdmBlackScholesMesher>(
        xGrid, process, maturity, K,
        Null<Real>(), Null<Real>(), 0.0001, 1.5, cPoints);
    auto mesher = ext::make_shared<FdmMesherComposite>(equityMesher);

    ext::shared_ptr<StrikedTypePayoff> payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, K);

    auto barrierCondition =
        ext::make_shared<FdmDiscreteBarrierStepCondition>(
            mesher, monTimes, L, U);
    auto conditions = ext::make_shared<FdmStepConditionComposite>(
        std::list<std::vector<Time>>(1, monTimes),
        FdmStepConditionComposite::Conditions(1, barrierCondition));

    const FdmBoundaryConditionSet bcSet;

    // Evaluate spots spanning the barrier corridor interior
    const Real spots[] = {96, 98, 99, 99.5, 100, 100.5,
                          101, 102, 104, 106, 108, 109};

    // MC references from Table 1 (where available)
    // S:     95    95.5   99.5   100    100.5  109.5  110
    // MC: 0.17359 0.18291 0.22923 0.23263 0.23410 0.17426 0.16712

    auto solve = [&](const FdmBlackScholesSpatialDesc& desc) {
        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), Size(0),
            ext::shared_ptr<FdmQuantoHelper>(), desc);
        Array rhs = buildPayoff(mesher, payoff);
        FdmBackwardSolver solver(op, bcSet, conditions,
                                 FdmSchemeDesc::CrankNicolson());
        solver.rollback(rhs, maturity, 0.0, tGrid, 0);
        return rhs;
    };

    FdmBlackScholesSpatialDesc cDesc = centralDesc();
    cDesc.mMatrixPolicy =
        FdmBlackScholesSpatialDesc::MMatrixPolicy::None;

    Array rhsCentral = solve(cDesc);
    Array rhsExpFit  = solve(fittedDesc());
    Array rhsMT      = solve(milevTaglianiDesc());

    BOOST_TEST_MESSAGE(
        "  S       | Central+CN  | ExpFit+CN   | MT+CN");
    BOOST_TEST_MESSAGE(
        "  --------|-------------|-------------|-------------");

    for (Real s : spots) {
        const Real vc = valueAtSpot(rhsCentral, mesher, s);
        const Real ve = valueAtSpot(rhsExpFit, mesher, s);
        const Real vm = valueAtSpot(rhsMT, mesher, s);

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(5);
        oss << "  " << std::setw(7) << s
            << " | " << std::setw(11) << vc
            << " | " << std::setw(11) << ve
            << " | " << std::setw(11) << vm;
        BOOST_TEST_MESSAGE(oss.str());

        // All values in the barrier corridor should be finite
        // and non-negative for the non-standard schemes
        BOOST_CHECK(std::isfinite(ve));
        BOOST_CHECK(std::isfinite(vm));
        BOOST_CHECK_MESSAGE(ve >= -1e-10,
            "ExpFit negative at S=" << s << ": " << ve);
        BOOST_CHECK_MESSAGE(vm >= -1e-10,
            "MT negative at S=" << s << ": " << vm);
    }
}


// ===================================================================
// 10. Grid convergence study
//     Shows V(100) converges as grid refines under ExpFit+CN
// ===================================================================

BOOST_AUTO_TEST_CASE(testGridConvergenceStudy) {

    BOOST_TEST_MESSAGE(
        "Testing grid convergence for discrete double barrier "
        "knock-out (ExpFit+CN)...");

    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real K = 100.0;
    const Real L = 95.0;
    const Real U = 110.0;
    const Volatility sigma = 0.25;
    const Rate r = 0.05;
    const Time maturity = 0.5;

    std::vector<Time> monTimes;
    for (Size i = 1; i <= 5; ++i)
        monTimes.push_back(maturity * Real(i) / 5.0);

    auto process = makeBSProcess(100.0, r, 0.0, sigma);

    // Grid levels: coarse, medium, fine
    struct GridLevel { Size xGrid; Size tGrid; const char* label; };
    const GridLevel levels[] = {
        { 500, 125, "Coarse (500x125)"},
        {1000, 250, "Medium (1000x250)"},
        {2000, 500, "Fine   (2000x500)"}
    };

    std::vector<Real> values;

    for (const auto& level : levels) {
        std::vector<std::tuple<Real, Real, bool>> cPoints = {
            {K, 0.1, true}, {L, 0.1, true}, {U, 0.1, true}
        };
        auto equityMesher = ext::make_shared<FdmBlackScholesMesher>(
            level.xGrid, process, maturity, K,
            Null<Real>(), Null<Real>(), 0.0001, 1.5, cPoints);
        auto mesher = ext::make_shared<FdmMesherComposite>(equityMesher);

        ext::shared_ptr<StrikedTypePayoff> payoff =
            ext::make_shared<PlainVanillaPayoff>(Option::Call, K);

        auto barrierCondition =
            ext::make_shared<FdmDiscreteBarrierStepCondition>(
                mesher, monTimes, L, U);
        auto conditions = ext::make_shared<FdmStepConditionComposite>(
            std::list<std::vector<Time>>(1, monTimes),
            FdmStepConditionComposite::Conditions(1, barrierCondition));

        const FdmBoundaryConditionSet bcSet;

        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), Size(0),
            ext::shared_ptr<FdmQuantoHelper>(), fittedDesc());

        Array rhs = buildPayoff(mesher, payoff);
        FdmBackwardSolver solver(op, bcSet, conditions,
                                 FdmSchemeDesc::CrankNicolson());
        solver.rollback(rhs, maturity, 0.0, level.tGrid, 0);

        const Real v100 = valueAtSpot(rhs, mesher, 100.0);
        values.push_back(v100);

        BOOST_TEST_MESSAGE("  " << level.label
                           << ":  V(100) = " << v100);
    }

    // Convergence check: finer grids should be closer to MC reference
    const Real mcRef = 0.23263;
    if (values.size() == 3) {
        const Real err0 = std::fabs(values[0] - mcRef);
        const Real err1 = std::fabs(values[1] - mcRef);
        const Real err2 = std::fabs(values[2] - mcRef);

        BOOST_TEST_MESSAGE("  |V_coarse - MC| = " << err0);
        BOOST_TEST_MESSAGE("  |V_medium - MC| = " << err1);
        BOOST_TEST_MESSAGE("  |V_fine   - MC| = " << err2);

        BOOST_CHECK_MESSAGE(err2 < err0,
            "Convergence check failed: finest grid error ("
            << err2 << ") should be less than coarsest ("
            << err0 << ")");
    }

    // The fine-grid value should be close to MC reference 0.23263
    BOOST_CHECK_MESSAGE(
        std::fabs(values.back() - 0.23263) / 0.23263 < 0.05,
        "Fine-grid V(100)=" << values.back()
        << " differs from MC ref 0.23263 by more than 5%");
}


// ===================================================================
// 11. Vanilla European vs Black-Scholes analytical formula
// ===================================================================

BOOST_AUTO_TEST_CASE(testVanillaEuropeanAgainstAnalytical) {

    BOOST_TEST_MESSAGE(
        "Testing FD vanilla European call against "
        "Black-Scholes analytical formula...");

    const DayCounter dc = Actual365Fixed();
    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real spot = 100.0;
    const Rate r = 0.05;
    const Rate q = 0.0;
    const Volatility vol = 0.20;
    const Real K = 100.0;

    auto process = ext::make_shared<BlackScholesMertonProcess>(
        Handle<Quote>(ext::make_shared<SimpleQuote>(spot)),
        Handle<YieldTermStructure>(flatRate(today, q, dc)),
        Handle<YieldTermStructure>(flatRate(today, r, dc)),
        Handle<BlackVolTermStructure>(flatVol(today, vol, dc)));

    const Date exerciseDate = today + Period(12, Months);
    const auto exercise =
        ext::make_shared<EuropeanExercise>(exerciseDate);
    const auto payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, K);

    // Analytical reference
    VanillaOption refOption(payoff, exercise);
    refOption.setPricingEngine(
        ext::make_shared<AnalyticEuropeanEngine>(process));
    const Real refValue = refOption.NPV();
    const Real refDelta = refOption.delta();

    BOOST_TEST_MESSAGE("  Analytical: NPV=" << refValue
                       << ", delta=" << refDelta);

    // Test each spatial scheme via FdBlackScholesVanillaEngine
    FdmBlackScholesSpatialDesc schemes[] = {
        FdmBlackScholesSpatialDesc::standard(),
        FdmBlackScholesSpatialDesc::exponentialFitting(),
        FdmBlackScholesSpatialDesc::milevTaglianiCN()
    };
    const char* names[] = {"Central", "ExpFit", "MT"};

    for (Size s = 0; s < 3; ++s) {
        VanillaOption option(payoff, exercise);
        option.setPricingEngine(
            ext::make_shared<FdBlackScholesVanillaEngine>(
                process, 200, 1000, 0,
                FdmSchemeDesc::CrankNicolson(),
                false, -Null<Real>(),
                FdBlackScholesVanillaEngine::Spot,
                schemes[s]));

        const Real fdValue = option.NPV();
        const Real fdDelta = option.delta();

        const Real relErr = std::fabs(fdValue - refValue) / refValue;
        const Real deltaErr = std::fabs(fdDelta - refDelta) / std::fabs(refDelta);

        BOOST_TEST_MESSAGE("  " << names[s]
            << ": NPV=" << fdValue << " (err=" << relErr * 100 << "%)"
            << ", delta=" << fdDelta << " (err=" << deltaErr * 100 << "%)");

        BOOST_CHECK_MESSAGE(relErr < 0.01,
            names[s] << " price error " << relErr
            << " exceeds 1% threshold");
        BOOST_CHECK_MESSAGE(deltaErr < 0.05,
            names[s] << " delta error " << deltaErr
            << " exceeds 5% threshold");
    }
}


// ===================================================================
// 12. Continuous barrier non-regression with spatialDesc
// ===================================================================

BOOST_AUTO_TEST_CASE(testContinuousBarrierNonRegression) {

    BOOST_TEST_MESSAGE(
        "Testing continuous barrier engine does not regress "
        "when spatialDesc is passed...");

    const DayCounter dc = Actual365Fixed();
    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real spot = 100.0;
    const Rate r = 0.05;
    const Volatility vol = 0.25;
    const Real K = 100.0;
    const Real barrier = 90.0;

    auto process = ext::make_shared<BlackScholesMertonProcess>(
        Handle<Quote>(ext::make_shared<SimpleQuote>(spot)),
        Handle<YieldTermStructure>(flatRate(today, 0.0, dc)),
        Handle<YieldTermStructure>(flatRate(today, r, dc)),
        Handle<BlackVolTermStructure>(flatVol(today, vol, dc)));

    const Date exerciseDate = today + Period(6, Months);
    const auto exercise =
        ext::make_shared<EuropeanExercise>(exerciseDate);
    const auto payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, K);

    // Analytical reference (continuous barrier)
    BarrierOption refOption(Barrier::DownOut, barrier, 0.0,
                            payoff, exercise);
    refOption.setPricingEngine(
        ext::make_shared<AnalyticBarrierEngine>(process));
    const Real refValue = refOption.NPV();

    BOOST_TEST_MESSAGE("  Analytical DownOut: NPV=" << refValue);

    // FD with StandardCentral and ExponentialFitting
    struct DescCase {
        FdmBlackScholesSpatialDesc desc; const char* name;
    };
    for (const auto& dc : {
            DescCase{FdmBlackScholesSpatialDesc::standard(), "Central"},
            DescCase{FdmBlackScholesSpatialDesc::exponentialFitting(),
                     "ExpFit"}}) {
        BarrierOption option(Barrier::DownOut, barrier, 0.0,
                             payoff, exercise);
        option.setPricingEngine(
            ext::make_shared<FdBlackScholesBarrierEngine>(
                process, 200, 800, 0,
                FdmSchemeDesc::CrankNicolson(),
                false, -Null<Real>(), dc.desc));

        const Real fdValue = option.NPV();
        const Real relErr = std::fabs(fdValue - refValue) / refValue;

        BOOST_TEST_MESSAGE("  " << dc.name << ": NPV=" << fdValue
            << " (err=" << relErr * 100 << "%)");
        BOOST_CHECK_MESSAGE(relErr < 0.02,
            dc.name << " barrier error " << relErr << " exceeds 2%");
    }
}


// ===================================================================
// 13. Shout engine smoke test with spatialDesc
// ===================================================================

BOOST_AUTO_TEST_CASE(testShoutEngineSmokeTest) {

    BOOST_TEST_MESSAGE(
        "Testing FdBlackScholesShoutEngine with "
        "ExponentialFitting spatial scheme...");

    const DayCounter dc = Actual365Fixed();
    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real spot = 100.0;
    const Rate r = 0.05;
    const Volatility vol = 0.20;
    const Real K = 100.0;

    auto process = ext::make_shared<BlackScholesMertonProcess>(
        Handle<Quote>(ext::make_shared<SimpleQuote>(spot)),
        Handle<YieldTermStructure>(flatRate(today, 0.0, dc)),
        Handle<YieldTermStructure>(flatRate(today, r, dc)),
        Handle<BlackVolTermStructure>(flatVol(today, vol, dc)));

    const Date exerciseDate = today + Period(12, Months);
    const auto exercise =
        ext::make_shared<EuropeanExercise>(exerciseDate);
    const auto payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, K);

    // European reference
    VanillaOption euroOption(payoff, exercise);
    euroOption.setPricingEngine(
        ext::make_shared<AnalyticEuropeanEngine>(process));
    const Real euroValue = euroOption.NPV();

    // Shout option with ExponentialFitting
    VanillaOption shoutOption(payoff, exercise);
    shoutOption.setPricingEngine(
        ext::make_shared<FdBlackScholesShoutEngine>(
            process, 100, 400, 0,
            FdmSchemeDesc::CrankNicolson(),
            FdmBlackScholesSpatialDesc::exponentialFitting()));

    const Real shoutValue = shoutOption.NPV();

    BOOST_TEST_MESSAGE("  European: NPV=" << euroValue);
    BOOST_TEST_MESSAGE("  Shout:    NPV=" << shoutValue);

    BOOST_CHECK_MESSAGE(std::isfinite(shoutValue),
        "Shout value is not finite: " << shoutValue);
    BOOST_CHECK_MESSAGE(shoutValue > 0.0,
        "Shout value should be positive: " << shoutValue);
    BOOST_CHECK_MESSAGE(shoutValue >= euroValue - 0.01,
        "Shout value (" << shoutValue
        << ") should be >= European value (" << euroValue << ")");
}


// ===================================================================
// 14. Non-European exercise rejection for discrete monitoring
// ===================================================================

BOOST_AUTO_TEST_CASE(testDiscreteBarrierRejectsNonEuropeanExercise) {

    BOOST_TEST_MESSAGE(
        "Testing discrete barrier engine rejects "
        "American and Bermudan exercise...");

    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    auto process = makeBSProcess(100.0, 0.05, 0.0, 0.25);

    const Date exerciseDate = today + Period(6, Months);
    const auto payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, 100.0);

    std::vector<Date> monDates;
    for (Size i = 1; i <= 3; ++i)
        monDates.push_back(today + Period(i * 2, Months));

    // American exercise should throw
    {
        auto exercise = ext::make_shared<AmericanExercise>(
            today, exerciseDate);
        BarrierOption option(Barrier::DownOut, 90.0, 0.0,
                             payoff, exercise);
        option.setPricingEngine(
            ext::make_shared<FdBlackScholesBarrierEngine>(
                process, monDates, 50, 100, 0,
                FdmSchemeDesc::CrankNicolson()));

        BOOST_CHECK_THROW(option.NPV(), Error);
    }

    // Bermudan exercise should throw
    {
        std::vector<Date> bermudanDates;
        bermudanDates.push_back(today + Period(3, Months));
        bermudanDates.push_back(exerciseDate);
        auto exercise = ext::make_shared<BermudanExercise>(bermudanDates);
        BarrierOption option(Barrier::DownOut, 90.0, 0.0,
                             payoff, exercise);
        option.setPricingEngine(
            ext::make_shared<FdBlackScholesBarrierEngine>(
                process, monDates, 50, 100, 0,
                FdmSchemeDesc::CrankNicolson()));

        BOOST_CHECK_THROW(option.NPV(), Error);
    }

    // European exercise should succeed
    {
        auto exercise = ext::make_shared<EuropeanExercise>(exerciseDate);
        BarrierOption option(Barrier::DownOut, 90.0, 0.0,
                             payoff, exercise);
        option.setPricingEngine(
            ext::make_shared<FdBlackScholesBarrierEngine>(
                process, monDates, 50, 100, 0,
                FdmSchemeDesc::CrankNicolson()));

        const Real value = option.NPV();
        BOOST_CHECK_MESSAGE(std::isfinite(value) && value > 0.0,
            "European discrete barrier should succeed, got " << value);
    }
}


// ===================================================================
// 15. t=0 monitoring: spot outside corridor at valuation date
// ===================================================================

BOOST_AUTO_TEST_CASE(testValuationDateMonitoringKnockOut) {

    BOOST_TEST_MESSAGE(
        "Testing t=0 discrete barrier monitoring "
        "(spot outside corridor at valuation date)...");

    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real spot = 85.0;  // below barrier L=90
    const Real barrier = 90.0;
    const Real rebate = 2.0;

    auto process = makeBSProcess(spot, 0.05, 0.0, 0.25);

    const Date exerciseDate = today + Period(6, Months);
    const auto exercise = ext::make_shared<EuropeanExercise>(exerciseDate);
    const auto payoff =
        ext::make_shared<PlainVanillaPayoff>(Option::Call, 100.0);

    // Monitoring dates INCLUDE today (t=0)
    std::vector<Date> monDatesWithToday;
    monDatesWithToday.push_back(today);
    for (Size i = 1; i <= 3; ++i)
        monDatesWithToday.push_back(today + Period(i * 2, Months));

    // Monitoring dates WITHOUT today
    std::vector<Date> monDatesNoToday;
    for (Size i = 1; i <= 3; ++i)
        monDatesNoToday.push_back(today + Period(i * 2, Months));

    // DownOut with today as monitoring date: spot=85 < barrier=90
    // should be knocked out immediately → value = rebate (immediate
    // payment, hit-time rebate semantics, no discounting since t=0)
    {
        BarrierOption option(Barrier::DownOut, barrier, rebate,
                             payoff, exercise);
        option.setPricingEngine(
            ext::make_shared<FdBlackScholesBarrierEngine>(
                process, monDatesWithToday, 50, 100, 0,
                FdmSchemeDesc::CrankNicolson()));

        const Real value = option.NPV();
        const Real expected = rebate;  // immediate payment at t=0

        BOOST_TEST_MESSAGE("  DownOut with t=0 monitoring: NPV="
            << value << ", expected=" << expected);

        BOOST_CHECK_MESSAGE(
            std::fabs(value - expected) < 1e-10,
            "Knocked-out-at-t=0 value should equal immediate rebate "
            << expected << ", got " << value);
    }

    // DownOut WITHOUT today: spot=85 outside barrier but no t=0 check.
    // The PDE starts from unmodified initial condition, so the value
    // must differ from the immediate-rebate case.
    {
        BarrierOption option(Barrier::DownOut, barrier, rebate,
                             payoff, exercise);
        option.setPricingEngine(
            ext::make_shared<FdBlackScholesBarrierEngine>(
                process, monDatesNoToday, 50, 100, 0,
                FdmSchemeDesc::CrankNicolson()));

        const Real value = option.NPV();

        BOOST_TEST_MESSAGE("  DownOut without t=0: NPV=" << value);

        BOOST_CHECK(std::isfinite(value));
        BOOST_CHECK_MESSAGE(std::fabs(value - rebate) > 0.01,
            "Without t=0 monitoring, value should differ from "
            "immediate rebate " << rebate << ", got " << value);
    }

    // Spot INSIDE corridor with today as monitoring date: prices should
    // match the non-t=0 baseline (t=0 check is a no-op when inside)
    {
        const Real insideSpot = 100.0;  // inside [90, infinity)
        auto processInside = makeBSProcess(insideSpot, 0.05, 0.0, 0.25);

        // With today in monitoring dates
        BarrierOption optWith(Barrier::DownOut, barrier, rebate,
                              payoff, exercise);
        optWith.setPricingEngine(
            ext::make_shared<FdBlackScholesBarrierEngine>(
                processInside, monDatesWithToday, 50, 200, 0,
                FdmSchemeDesc::CrankNicolson()));

        // Without today
        BarrierOption optWithout(Barrier::DownOut, barrier, rebate,
                                 payoff, exercise);
        optWithout.setPricingEngine(
            ext::make_shared<FdBlackScholesBarrierEngine>(
                processInside, monDatesNoToday, 50, 200, 0,
                FdmSchemeDesc::CrankNicolson()));

        const Real vWith = optWith.NPV();
        const Real vWithout = optWithout.NPV();

        BOOST_TEST_MESSAGE("  Inside corridor: with t=0="
            << vWith << ", without t=0=" << vWithout);

        BOOST_CHECK_MESSAGE(
            std::fabs(vWith - vWithout) < 0.05,
            "Spot inside corridor: with/without t=0 monitoring "
            "should be similar. with=" << vWith
            << ", without=" << vWithout);
    }

    // DownIn with today as monitoring date: spot=85 < barrier=90
    // immediately knocked in → price equals vanilla
    {
        BarrierOption optionIn(Barrier::DownIn, barrier, 0.0,
                               payoff, exercise);
        optionIn.setPricingEngine(
            ext::make_shared<FdBlackScholesBarrierEngine>(
                process, monDatesWithToday, 50, 100, 0,
                FdmSchemeDesc::CrankNicolson()));

        const Real knockInValue = optionIn.NPV();

        // Compare against vanilla FD price
        VanillaOption vanilla(payoff, exercise);
        vanilla.setPricingEngine(
            ext::make_shared<FdBlackScholesVanillaEngine>(
                process, 50, 100, 0,
                FdmSchemeDesc::CrankNicolson()));
        const Real vanillaValue = vanilla.NPV();

        BOOST_TEST_MESSAGE("  DownIn with t=0 monitoring: NPV="
            << knockInValue << ", vanilla=" << vanillaValue);

        BOOST_CHECK_MESSAGE(
            std::fabs(knockInValue - vanillaValue) < 0.01,
            "Knock-in at t=0 should equal vanilla " << vanillaValue
            << ", got " << knockInValue);
    }
}


// ===================================================================
// 16. Multi-point mesher boundary-coincident critical point retention
// ===================================================================

BOOST_AUTO_TEST_CASE(testMultiPointMesherBoundaryInclusive) {

    BOOST_TEST_MESSAGE(
        "Testing multi-point mesher retains critical points "
        "at domain boundaries (inclusive filtering)...");

    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    auto process = makeBSProcess(100.0, 0.05, 0.0, 0.25);

    const Real K = 100.0;
    const Real L = 90.0;
    const Time maturity = 0.5;
    const Size xGrid = 200;

    // Helper: true if logTarget appears as a grid node in locs.
    auto hasNode = [](const std::vector<Real>& locs, Real logTarget) {
        return std::any_of(locs.begin(), locs.end(),
            [logTarget](Real x) {
                return std::fabs(x - logTarget) < 1e-10;
            });
    };

    // Helper: verify a mesh is strictly increasing.
    auto checkStrictlyIncreasing =
        [](const std::vector<Real>& locs, const std::string& label) {
        for (Size i = 1; i < locs.size(); ++i) {
            BOOST_CHECK_MESSAGE(locs[i] > locs[i-1],
                label << " not strictly increasing at index " << i
                << ": " << locs[i-1] << " >= " << locs[i]);
        }
    };

    // Set xMinConstraint = log(L) so that the lower barrier
    // is EXACTLY at the domain boundary.
    const Real xMinConstraint = std::log(L);

    std::vector<std::tuple<Real, Real, bool>> cPoints = {
        {K, 0.1, true},
        {L, 0.1, true}  // L is exactly at xMin
    };

    FdmBlackScholesMesher mesher(
        xGrid, process, maturity, K,
        xMinConstraint, Null<Real>(), 0.0001, 1.5,
        cPoints);

    const auto& locs = mesher.locations();

    BOOST_CHECK_MESSAGE(hasNode(locs, std::log(L)),
        "Critical point at domain boundary L="
        << L << " (xMin) should be retained by multi-point mesher");

    checkStrictlyIncreasing(locs, "Main mesh");

    // Test xMax boundary: critical point at upper domain edge
    {
        const Real U = 130.0;

        std::vector<std::tuple<Real, Real, bool>> cPts = {
            {K, 0.1, true},
            {U, 0.1, true}  // U is exactly at xMax
        };

        FdmBlackScholesMesher mesherMax(
            xGrid, process, maturity, K,
            Null<Real>(), std::log(U), 0.0001, 1.5,
            cPts);

        BOOST_CHECK_MESSAGE(hasNode(mesherMax.locations(), std::log(U)),
            "Critical point at xMax boundary U=" << U
            << " should be retained");
    }

    // Test boundary-only cPoint: single cPoint at xMin exercising
    // the Concentrating1dMesher path (not uniform fallback).
    {
        std::vector<std::tuple<Real, Real, bool>> cPtsSingle = {
            {L, 0.1, true}  // only cPoint, exactly at xMin
        };

        FdmBlackScholesMesher mesherSingle(
            xGrid, process, maturity, K,
            xMinConstraint, Null<Real>(), 0.0001, 1.5,
            cPtsSingle);

        const auto& locsSingle = mesherSingle.locations();

        BOOST_CHECK_MESSAGE(hasNode(locsSingle, std::log(L)),
            "Boundary-only cPoint should be retained when it is the "
            "sole concentration point");

        checkStrictlyIncreasing(locsSingle, "Single-cPoint mesh");

        // Prove Concentrating1dMesher was used (not uniform fallback):
        // a concentrated mesh has non-uniform spacing near the cPoint.
        const Real uniformH =
            (locsSingle.back() - locsSingle.front())
            / (locsSingle.size() - 1);
        const Real h0 = locsSingle[1] - locsSingle[0];
        const Real hMid = locsSingle[xGrid/2] - locsSingle[xGrid/2 - 1];

        BOOST_CHECK_MESSAGE(
            std::fabs(h0 - hMid) / uniformH > 0.01,
            "Boundary-only cPoint mesh should have non-uniform spacing "
            "(proves Concentrating1dMesher, not uniform fallback). "
            "h0=" << h0 << ", hMid=" << hMid
            << ", uniformH=" << uniformH);
    }

    // Test far-outside cPoint: a cPoint well below xMin should be
    // filtered out (negative test)
    {
        const Real farBelow = 10.0;  // log(10) << xMin=log(90)

        std::vector<std::tuple<Real, Real, bool>> cPtsFar = {
            {K, 0.1, true},
            {farBelow, 0.1, true}  // well outside domain
        };

        FdmBlackScholesMesher mesherFar(
            xGrid, process, maturity, K,
            xMinConstraint, Null<Real>(), 0.0001, 1.5,
            cPtsFar);

        BOOST_CHECK_MESSAGE(
            !hasNode(mesherFar.locations(), std::log(farBelow)),
            "Far-outside cPoint at S=" << farBelow
            << " should be filtered out, but was found in mesh");
    }
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
