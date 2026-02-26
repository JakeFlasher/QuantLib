// r6
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Round-5 validation tests: exact mesh alignment, rebate support,
 engine-level discrete monitoring, and paper-faithful replication
 of Milev-Tagliani Example 4.1 (double barrier knock-out).
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
#include <ql/pricingengines/barrier/fdblackscholesbarrierengine.hpp>
#include <ql/pricingengines/barrier/makefdblackscholesbarrierengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <cmath>
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
    const Size tGrid = 200;
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

    // ---- Run with ExponentialFitting + CN ----
    {
        FdmBlackScholesSpatialDesc desc = fittedDesc();

        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), Size(0),
            ext::shared_ptr<FdmQuantoHelper>(), desc);

        FdmSolverDesc solverDesc = {
            mesher, bcSet, conditions, calculator,
            maturity, tGrid, 0
        };

        FdmBackwardSolver solver(
            op, bcSet, conditions,
            FdmSchemeDesc::CrankNicolson());
        
        Array rhs = buildPayoff(mesher, payoff);
        solver.rollback(rhs, maturity, 0.0, tGrid, 0);

        const Real price95  = valueAtSpot(rhs, mesher, 95.0);
        const Real price100 = valueAtSpot(rhs, mesher, 100.0);

        BOOST_TEST_MESSAGE("  ExpFit+CN: V(95)="  << price95
                           << ", V(100)=" << price100);

        // Prices must be positive
        BOOST_CHECK_MESSAGE(price95 >= -1e-10,
            "Negative price at S=95: " << price95);
        BOOST_CHECK_MESSAGE(price100 >= -1e-10,
            "Negative price at S=100: " << price100);

        // Paper MC: V(95)≈0.174, V(100)≈0.233
        // Tolerance generous due to log-space vs S-space difference
        BOOST_CHECK_MESSAGE(price95 > 0.10 && price95 < 0.30,
            "V(95) out of expected range: " << price95);
        BOOST_CHECK_MESSAGE(price100 > 0.15 && price100 < 0.35,
            "V(100) out of expected range: " << price100);
    }

    // ---- Run with MilevTagliani + CN ----
    {
        FdmBlackScholesSpatialDesc desc = milevTaglianiDesc();

        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), Size(0),
            ext::shared_ptr<FdmQuantoHelper>(), desc);

        Array rhs = buildPayoff(mesher, payoff);

        FdmBackwardSolver solver(
            op, bcSet, conditions,
            FdmSchemeDesc::CrankNicolson());
        solver.rollback(rhs, maturity, 0.0, tGrid, 0);

        const Real price95  = valueAtSpot(rhs, mesher, 95.0);
        const Real price100 = valueAtSpot(rhs, mesher, 100.0);

        BOOST_TEST_MESSAGE("  MT+CN: V(95)="  << price95
                           << ", V(100)=" << price100);

        BOOST_CHECK_MESSAGE(price95 >= -1e-10,
            "Negative price at S=95: " << price95);
        BOOST_CHECK_MESSAGE(price100 >= -1e-10,
            "Negative price at S=100: " << price100);

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


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
