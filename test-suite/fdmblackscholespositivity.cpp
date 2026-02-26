// r6
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Tests that non-standard spatial discretisations (exponential fitting,
 Milev-Tagliani CN-variant) restore positivity and eliminate spurious
 oscillations for discontinuous payoffs under extreme low-volatility
 conditions.  Also validates the Scheme-2 time-scheme gating in the
 solver layer.
*/

#include "toplevelfixture.hpp"
#include "utilities.hpp"

#include <ql/instruments/payoffs.hpp>
#include <ql/methods/finitedifferences/meshers/uniform1dmesher.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmeshercomposite.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesop.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp>
#include <ql/methods/finitedifferences/solvers/fdmsolverdesc.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp>
#include <ql/methods/finitedifferences/utilities/fdmdirichletboundary.hpp>
#include <ql/methods/finitedifferences/utilities/fdminnervaluecalculator.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <cmath>
#include <list>

using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FdmBlackScholesPositivityTests)

namespace {

    class TruncatedCallPayoff : public Payoff {
      public:
        TruncatedCallPayoff(Real strike, Real upper)
        : strike_(strike), upper_(upper) {}

        std::string name() const override {
            return "TruncatedCallPayoff";
        }
        std::string description() const override {
            return "TruncatedCallPayoff";
        }
        Real operator()(Real s) const override {
            return (s >= strike_ && s <= upper_) ? s - strike_ : 0.0;
        }
      private:
        Real strike_, upper_;
    };

    struct LowVolSetup {
        DayCounter dc;
        Date today;
        Real spot;
        Rate r, q;
        Volatility vol;
        Time maturity;
        ext::shared_ptr<GeneralizedBlackScholesProcess> process;

        LowVolSetup()
        : dc(Actual365Fixed()),
          today(Date(28, March, 2004)),
          spot(60.0), r(0.05), q(0.0), vol(0.001),
          maturity(5.0 / 12.0) {

            Settings::instance().evaluationDate() = today;
            process = ext::make_shared<BlackScholesMertonProcess>(
                Handle<Quote>(ext::make_shared<SimpleQuote>(spot)),
                Handle<YieldTermStructure>(flatRate(today, q, dc)),
                Handle<YieldTermStructure>(flatRate(today, r, dc)),
                Handle<BlackVolTermStructure>(flatVol(today, vol, dc)));
        }
    };

    Array buildPayoffVector(const ext::shared_ptr<FdmMesher>& mesher,
                            const ext::shared_ptr<FdmInnerValueCalculator>& calc,
                            Time maturity) {
        Array rhs(mesher->layout()->size());
        for (const auto& iter : *mesher->layout()) {
            rhs[iter.index()] = calc->avgInnerValue(iter, maturity);
        }
        return rhs;
    }

    void rollbackImplicitEuler(
            const ext::shared_ptr<FdmLinearOpComposite>& op,
            const FdmBoundaryConditionSet& bcSet,
            Array& rhs,
            Time maturity,
            Size timeSteps) {
        ext::shared_ptr<FdmStepConditionComposite> noCondition;
        FdmBackwardSolver(op, bcSet, noCondition,
                          FdmSchemeDesc::ImplicitEuler())
            .rollback(rhs, maturity, 0.0, timeSteps, 0);
    }

    void rollbackCrankNicolson(
            const ext::shared_ptr<FdmLinearOpComposite>& op,
            const FdmBoundaryConditionSet& bcSet,
            Array& rhs,
            Time maturity,
            Size timeSteps) {
        ext::shared_ptr<FdmStepConditionComposite> noCondition;
        FdmBackwardSolver(op, bcSet, noCondition,
                          FdmSchemeDesc::CrankNicolson())
            .rollback(rhs, maturity, 0.0, timeSteps, 0);
    }

    FdmBlackScholesSpatialDesc centralDesc() {
        FdmBlackScholesSpatialDesc d;
        d.scheme = FdmBlackScholesSpatialDesc::Scheme::StandardCentral;
        d.mMatrixPolicy = FdmBlackScholesSpatialDesc::MMatrixPolicy::None;
        return d;
    }

    FdmBlackScholesSpatialDesc fittedDesc() {
        FdmBlackScholesSpatialDesc d =
            FdmBlackScholesSpatialDesc::exponentialFitting();
        d.mMatrixPolicy = FdmBlackScholesSpatialDesc::MMatrixPolicy::None;
        return d;
    }

    FdmBlackScholesSpatialDesc milevTaglianiDesc() {
        FdmBlackScholesSpatialDesc d =
            FdmBlackScholesSpatialDesc::milevTaglianiCN();
        d.mMatrixPolicy = FdmBlackScholesSpatialDesc::MMatrixPolicy::None;
        return d;
    }

    Size countMonotonicityViolations(const Array& v,
                                     Size skip, Size n, Real tol) {
        Size violations = 0;
        for (Size i = skip; i + 1 < n - skip; ++i) {
            if (v[i + 1] > v[i] + tol)
                ++violations;
        }
        return violations;
    }

} // anonymous namespace


BOOST_AUTO_TEST_CASE(testExponentialFittingRestoresPositivity) {

    BOOST_TEST_MESSAGE(
        "Testing exponential fitting restores positivity for "
        "low-vol truncated call (discontinuous payoff)...");

    LowVolSetup setup;

    const Real K = 50.0;
    const Real U = 70.0;
    const Size xGrid = 400;
    const Size tGrid = 25;

    const Real xMin = std::log(1.0);
    const Real xMax = std::log(140.0);

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    const ext::shared_ptr<Payoff> payoff =
        ext::make_shared<TruncatedCallPayoff>(K, U);
    const auto calc =
        ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);

    const Array payoffVec = buildPayoffVector(mesher, calc, setup.maturity);

    const FdmBoundaryConditionSet bcSet = {
        ext::make_shared<FdmDirichletBoundary>(
            mesher, 0.0, 0, FdmDirichletBoundary::Lower),
        ext::make_shared<FdmDirichletBoundary>(
            mesher, 0.0, 0, FdmDirichletBoundary::Upper)
    };

    // ---- Case A: StandardCentral (expect negativity) ----
    {
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), centralDesc());

        Array v = payoffVec;
        rollbackImplicitEuler(op, bcSet, v, setup.maturity, tGrid);

        Size negCount = 0;
        for (Size i = 0; i < v.size(); ++i) {
            BOOST_CHECK(std::isfinite(v[i]));
            if (v[i] < -1e-4)
                ++negCount;
        }
        BOOST_CHECK_MESSAGE(negCount > 0,
            "StandardCentral should produce negative values "
            "for low-vol truncated call, but found none");
        BOOST_TEST_MESSAGE("  StandardCentral: "
            << negCount << " grid points below -1e-4");
    }

    // ---- Case B: ExponentialFitting (expect non-negative & bounded) ----
    {
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), fittedDesc());

        Array v = payoffVec;
        rollbackImplicitEuler(op, bcSet, v, setup.maturity, tGrid);

        const Real maxPayoff = U - K;
        for (Size i = 0; i < v.size(); ++i) {
            BOOST_CHECK_MESSAGE(std::isfinite(v[i]),
                "Non-finite value at node " << i);
            BOOST_CHECK_MESSAGE(v[i] >= -1e-12,
                "Positivity violation at node " << i
                << ": value = " << v[i]);
            BOOST_CHECK_MESSAGE(v[i] <= maxPayoff + 1e-10,
                "Upper bound violation at node " << i
                << ": value = " << v[i]
                << ", max payoff = " << maxPayoff);
        }
    }
}


BOOST_AUTO_TEST_CASE(testExponentialFittingEliminatesOscillations) {

    BOOST_TEST_MESSAGE(
        "Testing exponential fitting eliminates spatial oscillations "
        "for low-vol cash-or-nothing digital put...");

    LowVolSetup setup;

    const Real K = 60.0;
    const Real cash = 10.0;
    const Size xGrid = 400;
    const Size tGrid = 40;

    const Real xMin = std::log(10.0);
    const Real xMax = std::log(200.0);

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    const ext::shared_ptr<StrikedTypePayoff> payoff =
        ext::make_shared<CashOrNothingPayoff>(Option::Put, K, cash);
    const auto calc =
        ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);

    const Array payoffVec = buildPayoffVector(mesher, calc, setup.maturity);

    const FdmBoundaryConditionSet bcSet;

    const Size skip = 5;
    const Real monoTol = 1e-10;

    // ---- StandardCentral (expect monotonicity violations) ----
    Size violCentral;
    {
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), centralDesc());

        Array v = payoffVec;
        rollbackImplicitEuler(op, bcSet, v, setup.maturity, tGrid);

        violCentral = countMonotonicityViolations(v, skip, xGrid, monoTol);
        BOOST_CHECK_MESSAGE(violCentral > 0,
            "StandardCentral should produce monotonicity violations "
            "for low-vol digital put, but found none");
        BOOST_TEST_MESSAGE("  StandardCentral monotonicity violations: "
            << violCentral);
    }

    // ---- ExponentialFitting (expect zero violations, positivity, bounds) ----
    {
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), fittedDesc());

        Array v = payoffVec;
        rollbackImplicitEuler(op, bcSet, v, setup.maturity, tGrid);

        const Size violFitted =
            countMonotonicityViolations(v, skip, xGrid, monoTol);
        BOOST_TEST_MESSAGE("  ExponentialFitting monotonicity violations: "
            << violFitted);

        BOOST_CHECK_MESSAGE(violFitted == 0,
            "ExponentialFitting should eliminate monotonicity "
            "violations, but found " << violFitted);

        for (Size i = skip; i < xGrid - skip; ++i) {
            BOOST_CHECK_MESSAGE(std::isfinite(v[i]),
                "Non-finite value at node " << i);
            BOOST_CHECK_MESSAGE(v[i] >= -1e-12,
                "Positivity violation at node " << i
                << ": value = " << v[i]);
            BOOST_CHECK_MESSAGE(v[i] <= cash + 1e-10,
                "Upper bound violation at node " << i
                << ": value = " << v[i]);
        }
    }
}


BOOST_AUTO_TEST_CASE(testMilevTaglianiEliminatesCNOscillations) {

    BOOST_TEST_MESSAGE(
        "Testing Milev-Tagliani CN effective diffusion eliminates "
        "oscillations under Crank-Nicolson for low-vol digital put...");

    LowVolSetup setup;

    const Real K = 60.0;
    const Real cash = 10.0;
    const Size xGrid = 400;
    const Size tGrid = 40;

    const Real xMin = std::log(10.0);
    const Real xMax = std::log(200.0);

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    const ext::shared_ptr<StrikedTypePayoff> payoff =
        ext::make_shared<CashOrNothingPayoff>(Option::Put, K, cash);
    const auto calc =
        ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);

    const Array payoffVec = buildPayoffVector(mesher, calc, setup.maturity);

    // Natural boundary conditions (no explicit Dirichlet) — skip
    // boundary nodes in value checks to avoid stencil artifacts.
    const FdmBoundaryConditionSet bcSet;

    const Size skip = 5;
    const Real monoTol = 1e-10;

    // ---- Case A: StandardCentral + Crank-Nicolson ----
    // With sigma=0.001, r=0.05 the cell Peclet number is ~375 on this
    // grid.  The explicit-side matrix N = I + 0.5*dt*L has negative
    // off-diagonals, and eigenvalues of P^{-1}N cluster near -1.
    // Expect: spurious negativity, overshoot, or monotonicity violations.
    {
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), centralDesc());

        Array v = payoffVec;
        rollbackCrankNicolson(op, bcSet, v, setup.maturity, tGrid);

        Size negCount = 0, overCount = 0;
        for (Size i = 0; i < v.size(); ++i) {
            BOOST_CHECK(std::isfinite(v[i]));
            if (v[i] < -1e-6)
                ++negCount;
            if (v[i] > cash + 1e-6)
                ++overCount;
        }

        const Size monoViol =
            countMonotonicityViolations(v, skip, xGrid, monoTol);

        BOOST_TEST_MESSAGE("  Central(CN): neg=" << negCount
                           << ", over=" << overCount
                           << ", monoViol=" << monoViol);

        // At least one failure mode must be present to validate the test
        BOOST_CHECK_MESSAGE(
            negCount > 0 || overCount > 0 || monoViol > 0,
            "StandardCentral + CN should produce artifacts "
            "for low-vol digital put, but found none");
    }

    // ---- Case B: MilevTagliani CN-variant + Crank-Nicolson ----
    // The modified reaction-term discretization adds effective diffusion
    // a_eff = v/2 + r^2*h^2/(8*v), which shifts L's off-diagonals to
    // non-negative, moving eigenvalues of P^{-1}N away from -1.
    // Expect: no negativity, no overshoot, no monotonicity violations.
    {
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), milevTaglianiDesc());

        Array v = payoffVec;
        rollbackCrankNicolson(op, bcSet, v, setup.maturity, tGrid);

        Size negCount = 0, overCount = 0;
        for (Size i = skip; i < xGrid - skip; ++i) {
            BOOST_CHECK_MESSAGE(std::isfinite(v[i]),
                "Non-finite value at node " << i);
            if (v[i] < -1e-12)
                ++negCount;
            if (v[i] > cash + 1e-10)
                ++overCount;
        }

        const Size monoViol =
            countMonotonicityViolations(v, skip, xGrid, monoTol);

        BOOST_TEST_MESSAGE("  MT(CN): neg=" << negCount
                           << ", over=" << overCount
                           << ", monoViol=" << monoViol);

        BOOST_CHECK_EQUAL(negCount, Size(0));
        BOOST_CHECK_EQUAL(overCount, Size(0));
        BOOST_CHECK_EQUAL(monoViol, Size(0));
    }
}


BOOST_AUTO_TEST_CASE(testScheme2GatingFallbackMatchesScheme1) {

    BOOST_TEST_MESSAGE(
        "Testing Scheme-2 gating falls back to exponential fitting "
        "when time stepping is Implicit Euler "
        "(values must match explicit Scheme-1 request)...");

    // This test exercises FdmBlackScholesSolver (not the operator
    // directly) because the CN-equivalence gating is enforced in
    // the solver's performCalculations().

    LowVolSetup setup;

    const Real K = 60.0;
    const Real cash = 10.0;
    const Size xGrid = 200;
    const Size tGrid = 50;

    const Real xMin = std::log(10.0);
    const Real xMax = std::log(200.0);

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    const ext::shared_ptr<StrikedTypePayoff> payoff =
        ext::make_shared<CashOrNothingPayoff>(Option::Put, K, cash);
    const auto calculator =
        ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);

    const auto conditions =
        ext::make_shared<FdmStepConditionComposite>(
            std::list<std::vector<Time> >(),
            FdmStepConditionComposite::Conditions());

    const FdmBoundaryConditionSet bcSet;

    const FdmSolverDesc solverDesc = {
        mesher, bcSet, conditions, calculator,
        setup.maturity, tGrid, 0
    };

    // Solver A: ExponentialFitting explicitly requested + Implicit Euler
    FdmBlackScholesSolver solverFitted(
        Handle<GeneralizedBlackScholesProcess>(setup.process),
        K, solverDesc,
        FdmSchemeDesc::ImplicitEuler(),
        false, -Null<Real>(),
        Handle<FdmQuantoHelper>(),
        fittedDesc());

    // Solver B: MilevTagliani requested + Implicit Euler
    // The solver must detect that IE is not CN-equivalent, silently
    // fall back to ExponentialFitting, and produce identical results.
    FdmBlackScholesSolver solverMT(
        Handle<GeneralizedBlackScholesProcess>(setup.process),
        K, solverDesc,
        FdmSchemeDesc::ImplicitEuler(),
        false, -Null<Real>(),
        Handle<FdmQuantoHelper>(),
        milevTaglianiDesc());

    // Evaluate at multiple spots spanning the full domain to make
    // the gating proof robust against accidental partial matches.
    const Real spots[] = { 20.0, 40.0, 55.0, 60.0, 65.0, 80.0, 120.0 };

    for (Real s : spots) {
        const Real v1 = solverFitted.valueAt(s);
        const Real v2 = solverMT.valueAt(s);

        BOOST_CHECK(std::isfinite(v1));
        BOOST_CHECK(std::isfinite(v2));

        // If gating works correctly, the fallback-to-Scheme-1 solver
        // must produce bit-identical results to the explicitly-
        // configured Scheme-1 solver (same operator, same rollback).
        BOOST_CHECK_SMALL(std::fabs(v1 - v2), 1e-10);
    }
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
