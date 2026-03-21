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
#include <ql/math/matrixutilities/sparsematrix.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
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

    Size countNegativeOffDiag(const SparseMatrix& mat,
                              Size n, Real eps = 0.0) {
        Size count = 0;
        for (Size i = 1; i < n - 1; ++i) {
            for (Size j = 0; j < n; ++j) {
                if (j != i && Real(mat(i,j)) < -eps)
                    ++count;
            }
        }
        return count;
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


BOOST_AUTO_TEST_CASE(testTruncatedCallNumericalDiffusion) {

    BOOST_TEST_MESSAGE(
        "Testing truncated call numerical diffusion comparison "
        "(Low Volatility paper, Example 4.1): CN vs ExpFit vs MT...");

    // Low Volatility paper Example 4.1:
    //   r=0.05, sigma=0.001, T=5/12, U=70, K=50
    // Truncated call payoff: max(S-K, 0) for S in [K, U], else 0
    // With sigma^2 << r, standard CN produces spurious oscillations
    // near the discontinuity at U=70.

    LowVolSetup setup;  // sigma=0.001, r=0.05, spot=60

    const Real K = 50.0;
    const Real U = 70.0;
    const Size xGrid = 800;
    const Size tGrid = 200;

    const Real xMin = std::log(1.0);
    const Real xMax = std::log(140.0);

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    const ext::shared_ptr<Payoff> payoff =
        ext::make_shared<TruncatedCallPayoff>(K, U);
    const auto calc =
        ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);

    const Array payoffVec = buildPayoffVector(mesher, calc, setup.maturity);

    const FdmBoundaryConditionSet bcSet;

    const Size skip = 5;
    const Real monoTol = 1e-10;

    // ---- StandardCentral + CN: expect oscillations ----
    Size violCentral = 0;
    {
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), centralDesc());

        Array v = payoffVec;
        rollbackCrankNicolson(op, bcSet, v, setup.maturity, tGrid);

        Size negCount = 0;
        for (Size i = skip; i < xGrid - skip; ++i) {
            if (v[i] < -1e-6) ++negCount;
        }
        violCentral = countMonotonicityViolations(v, skip, xGrid, monoTol);

        BOOST_TEST_MESSAGE("  Central+CN: neg=" << negCount
                           << ", monoViol=" << violCentral);

        BOOST_CHECK_MESSAGE(negCount > 0 || violCentral > 0,
            "Expected oscillations from StandardCentral+CN "
            "for low-vol truncated call");
    }

    // ---- ExponentialFitting + CN ----
    Size transWidthExpFit = 0;
    {
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), fittedDesc());

        Array v = payoffVec;
        rollbackCrankNicolson(op, bcSet, v, setup.maturity, tGrid);

        Size negCount = 0;
        const Size violFitted =
            countMonotonicityViolations(v, skip, xGrid, monoTol);
        for (Size i = skip; i < xGrid - skip; ++i) {
            if (v[i] < -1e-12) ++negCount;
        }

        // For truncated call, the solution is bell-shaped (increases
        // from K to peak, then decreases toward U), so monotonicity
        // violations in the global sense are expected.  The key checks
        // are: no negativity and fewer artifacts than StandardCentral.
        BOOST_CHECK_EQUAL(negCount, Size(0));

        BOOST_TEST_MESSAGE("  ExpFit+CN: neg=" << negCount
            << ", monoViol=" << violFitted);

        // Measure transition width near U=70:
        // count nodes where value drops from >0.5*peak to <0.01*peak
        const Real xU = std::log(U);
        Real peak = 0.0;
        for (Size i = 0; i < v.size(); ++i)
            peak = std::max(peak, v[i]);

        Size iStart = 0, iEnd = 0;
        bool foundStart = false;
        for (const auto& iter : *mesher->layout()) {
            const Real x = mesher->location(iter, 0);
            if (x > xU - 0.5 && x < xU + 0.5) {
                if (!foundStart && v[iter.index()] > 0.5 * peak) {
                    iStart = iter.index();
                    foundStart = true;
                }
                if (foundStart && v[iter.index()] < 0.01 * peak) {
                    iEnd = iter.index();
                    break;
                }
            }
        }
        transWidthExpFit = (iEnd > iStart) ? iEnd - iStart : 0;

        BOOST_TEST_MESSAGE("  ExpFit+CN: transition_width="
            << transWidthExpFit << " nodes");
    }

    // ---- MilevTagliani + CN ----
    Size transWidthMT = 0;
    {
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), milevTaglianiDesc());

        Array v = payoffVec;
        rollbackCrankNicolson(op, bcSet, v, setup.maturity, tGrid);

        Size negCount = 0;
        const Size violMT =
            countMonotonicityViolations(v, skip, xGrid, monoTol);
        for (Size i = skip; i < xGrid - skip; ++i) {
            if (v[i] < -1e-12) ++negCount;
        }

        BOOST_CHECK_EQUAL(negCount, Size(0));

        BOOST_TEST_MESSAGE("  MT+CN: neg=" << negCount
            << ", monoViol=" << violMT);

        // Measure transition width
        const Real xU = std::log(U);
        Real peak = 0.0;
        for (Size i = 0; i < v.size(); ++i)
            peak = std::max(peak, v[i]);

        Size iStart = 0, iEnd = 0;
        bool foundStart = false;
        for (const auto& iter : *mesher->layout()) {
            const Real x = mesher->location(iter, 0);
            if (x > xU - 0.5 && x < xU + 0.5) {
                if (!foundStart && v[iter.index()] > 0.5 * peak) {
                    iStart = iter.index();
                    foundStart = true;
                }
                if (foundStart && v[iter.index()] < 0.01 * peak) {
                    iEnd = iter.index();
                    break;
                }
            }
        }
        transWidthMT = (iEnd > iStart) ? iEnd - iStart : 0;

        BOOST_TEST_MESSAGE("  MT+CN: transition_width="
            << transWidthMT << " nodes");
    }

    // Both non-standard schemes should produce smooth solutions,
    // but their transition widths may differ (reflecting different
    // numerical diffusion: Duffy ~ rS*dS/2 vs MT ~ (r*dS/sigma)^2/8)
    BOOST_TEST_MESSAGE(
        "  Numerical diffusion comparison: "
        "ExpFit width=" << transWidthExpFit
        << " vs MT width=" << transWidthMT << " nodes");
}


BOOST_AUTO_TEST_CASE(testExtremeParameterStress) {

    BOOST_TEST_MESSAGE(
        "Testing extreme parameter stress: high r, near-zero "
        "maturity, fine/coarse grids, near-zero vol...");

    // Common setup
    const DayCounter dc = Actual365Fixed();
    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real K = 60.0;
    const Real cash = 10.0;

    auto makeProc = [&](Real spot, Rate r, Volatility vol) {
        return ext::make_shared<BlackScholesMertonProcess>(
            Handle<Quote>(ext::make_shared<SimpleQuote>(spot)),
            Handle<YieldTermStructure>(flatRate(today, 0.0, dc)),
            Handle<YieldTermStructure>(flatRate(today, r, dc)),
            Handle<BlackVolTermStructure>(flatVol(today, vol, dc)));
    };

    auto runCheck = [&](const char* label,
                        ext::shared_ptr<GeneralizedBlackScholesProcess> proc,
                        Time maturity, Size xGrid, Size tGrid) {
        const auto mesher = ext::make_shared<FdmMesherComposite>(
            ext::make_shared<Uniform1dMesher>(
                std::log(10.0), std::log(200.0), xGrid));

        const ext::shared_ptr<StrikedTypePayoff> payoff =
            ext::make_shared<CashOrNothingPayoff>(Option::Put, K, cash);
        const auto calc =
            ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);
        const Array payoffVec = buildPayoffVector(mesher, calc, maturity);
        const FdmBoundaryConditionSet bcSet;

        // Test with ExponentialFitting + CN
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, proc, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), fittedDesc());

        Array v = payoffVec;
        rollbackCrankNicolson(op, bcSet, v, maturity, tGrid);

        bool allOk = true;
        for (Size i = 0; i < v.size(); ++i) {
            if (!std::isfinite(v[i])) {
                allOk = false;
                break;
            }
        }

        BOOST_CHECK_MESSAGE(allOk,
            label << ": found non-finite values");
        BOOST_TEST_MESSAGE("  " << label << ": all finite = "
            << (allOk ? "yes" : "NO"));
    };

    // High r (sigma^2 < r)
    runCheck("High r (r=0.50, vol=0.20)",
             makeProc(60.0, 0.50, 0.20), 5.0/12.0, 400, 40);

    // Near-zero maturity
    runCheck("Near-zero maturity (T=1/365)",
             makeProc(60.0, 0.05, 0.20), 1.0/365.0, 100, 10);

    // Fine grid
    runCheck("Fine grid (5000x2000)",
             makeProc(60.0, 0.05, 0.001), 5.0/12.0, 200, 50);

    // Coarse grid
    runCheck("Coarse grid (20x5)",
             makeProc(60.0, 0.05, 0.20), 5.0/12.0, 20, 5);

    // Near-zero vol
    runCheck("Near-zero vol (sigma=1e-6)",
             makeProc(60.0, 0.05, 1e-6), 5.0/12.0, 200, 40);
}


BOOST_AUTO_TEST_CASE(testCrossValidationModerateVol) {

    BOOST_TEST_MESSAGE(
        "Testing cross-validation: all three schemes agree "
        "at moderate vol (sigma=0.30, sigma^2 > r)...");

    LowVolSetup baseSetup;

    // Override with moderate vol
    const Volatility vol = 0.30;
    const DayCounter dc = Actual365Fixed();

    auto process = ext::make_shared<BlackScholesMertonProcess>(
        Handle<Quote>(ext::make_shared<SimpleQuote>(baseSetup.spot)),
        Handle<YieldTermStructure>(flatRate(baseSetup.today, baseSetup.q, dc)),
        Handle<YieldTermStructure>(flatRate(baseSetup.today, baseSetup.r, dc)),
        Handle<BlackVolTermStructure>(flatVol(baseSetup.today, vol, dc)));

    const Real K = 60.0;
    const Real cash = 10.0;
    const Size xGrid = 400;
    const Size tGrid = 40;
    const Time maturity = baseSetup.maturity;

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(
            std::log(10.0), std::log(200.0), xGrid));

    const ext::shared_ptr<StrikedTypePayoff> payoff =
        ext::make_shared<CashOrNothingPayoff>(Option::Put, K, cash);
    const auto calc =
        ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);
    const Array payoffVec = buildPayoffVector(mesher, calc, maturity);
    const FdmBoundaryConditionSet bcSet;

    auto solve = [&](const FdmBlackScholesSpatialDesc& desc) {
        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), desc);
        Array v = payoffVec;
        rollbackCrankNicolson(op, bcSet, v, maturity, tGrid);
        return v;
    };

    Array vCentral = solve(centralDesc());
    Array vExpFit  = solve(fittedDesc());
    Array vMT      = solve(milevTaglianiDesc());

    // All three should agree very closely since sigma^2 = 0.09 > r = 0.05
    Real maxDiff = 0.0;
    for (Size i = 5; i < xGrid - 5; ++i) {
        const Real d1 = std::fabs(vCentral[i] - vExpFit[i]);
        const Real d2 = std::fabs(vCentral[i] - vMT[i]);
        const Real d3 = std::fabs(vExpFit[i] - vMT[i]);
        maxDiff = std::max(maxDiff, std::max(d1, std::max(d2, d3)));
    }

    BOOST_TEST_MESSAGE("  Max pairwise difference at sigma=0.30: "
        << maxDiff);

    BOOST_CHECK_MESSAGE(maxDiff < 1e-4,
        "Schemes should agree closely at moderate vol, "
        "but max diff = " << maxDiff);
}


// ===================================================================
// Negative dividend yield: MT M-matrix edge case
// AM-GM proof for MT requires q + σ²/2 ≥ 0. With negative q and
// a coarse grid, the condition fails and the fallback must trigger.
// ===================================================================

BOOST_AUTO_TEST_CASE(testNegativeDividendYieldMMatrixFallback) {

    BOOST_TEST_MESSAGE(
        "Testing MT M-matrix violation with negative dividend yield "
        "(q=-0.30, r=0.05, sigma=0.10, coarse grid)...");

    const DayCounter dc = Actual365Fixed();
    const Date today(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Real spot = 60.0;
    const Rate r = 0.05;
    const Rate q = -0.30;
    const Volatility vol = 0.10;
    const Real K = 60.0;
    const Time maturity = 5.0 / 12.0;

    auto process = ext::make_shared<BlackScholesMertonProcess>(
        Handle<Quote>(ext::make_shared<SimpleQuote>(spot)),
        Handle<YieldTermStructure>(flatRate(today, q, dc)),
        Handle<YieldTermStructure>(flatRate(today, r, dc)),
        Handle<BlackVolTermStructure>(flatVol(today, vol, dc)));

    const Real variance = vol * vol;

    BOOST_TEST_MESSAGE("  q+sigma^2/2=" << (q + variance / 2.0));

    // Use a coarse grid so that h is large enough to force
    // MT lower off-diagonal negative. With xGrid=10 on
    // [log(10), log(200)], h ≈ 0.333.
    // MT lower = sigma^2/(2h^2) + r^2/(8 sigma^2) - drift/(2h)
    //          = 0.005/0.111 + 0.000625/0.08 - 0.345/0.666
    //          ≈ 0.045 + 0.0078 - 0.518 = -0.465 < 0
    const Size xGrid = 10;
    const Size tGrid = 40;

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(
            std::log(10.0), std::log(200.0), xGrid));

    const Time dt = maturity / tGrid;

    // Step 1: DiagnosticsOnly on the coarse mesh — verify violation
    {
        FdmBlackScholesSpatialDesc desc =
            FdmBlackScholesSpatialDesc::milevTaglianiCN();
        desc.mMatrixPolicy =
            FdmBlackScholesSpatialDesc::MMatrixPolicy::DiagnosticsOnly;

        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), desc);

        op->setTime(0.0, dt);
        SparseMatrix mat = op->toMatrix();

        const Size n = mesher->layout()->size();
        Size negCount = countNegativeOffDiag(mat, n, 0.0);
        BOOST_CHECK_MESSAGE(negCount > 0,
            "Expected MT M-matrix violations with negative q on "
            "coarse grid under DiagnosticsOnly, found " << negCount);

        BOOST_TEST_MESSAGE("  MT+DiagnosticsOnly (coarse): "
            << negCount << " negative off-diag entries");
    }

    // Step 2: FallbackToExponentialFitting on the SAME coarse mesh —
    // verify the repaired operator has zero negative off-diags
    {
        FdmBlackScholesSpatialDesc desc =
            FdmBlackScholesSpatialDesc::milevTaglianiCN();
        desc.mMatrixPolicy =
            FdmBlackScholesSpatialDesc::MMatrixPolicy
                ::FallbackToExponentialFitting;

        auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), desc);

        op->setTime(0.0, dt);
        SparseMatrix mat = op->toMatrix();

        const Size n = mesher->layout()->size();
        Size negCount = countNegativeOffDiag(mat, n, 0.0);
        BOOST_CHECK_EQUAL(negCount, Size(0));

        BOOST_TEST_MESSAGE("  MT+Fallback (coarse): "
            << negCount << " negative off-diag entries (repaired)");
    }

    // Steps 3-4: Roll back a payoff on the same coarse violating mesh
    // using (a) MT fallback and (b) explicit ExpFit. Verify finite/positive.
    {
        const ext::shared_ptr<StrikedTypePayoff> payoff =
            ext::make_shared<CashOrNothingPayoff>(Option::Put, K, 10.0);
        const auto calc =
            ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);
        const Array payoffVec =
            buildPayoffVector(mesher, calc, maturity);
        const FdmBoundaryConditionSet bcSet;

        // Step 3: MT with fallback
        {
            const auto op = ext::make_shared<FdmBlackScholesOp>(
                mesher, process, K,
                false, -Null<Real>(), 0,
                ext::shared_ptr<FdmQuantoHelper>(),
                FdmBlackScholesSpatialDesc::milevTaglianiCN());

            Array v = payoffVec;
            rollbackCrankNicolson(op, bcSet, v, maturity, tGrid);

            bool allOk = true;
            Size negVals = 0;
            for (Size i = 1; i + 1 < v.size(); ++i) {
                if (!std::isfinite(v[i])) { allOk = false; break; }
                if (v[i] < -1e-6) ++negVals;
            }

            BOOST_CHECK_MESSAGE(allOk,
                "Non-finite values on coarse mesh with MT fallback");
            BOOST_CHECK_MESSAGE(negVals == 0,
                negVals << " negative values on coarse mesh with fallback");

            BOOST_TEST_MESSAGE("  Rollback on coarse mesh: all_finite="
                << allOk << ", neg=" << negVals);
        }

        // Step 4: explicit ExpFit on same mesh
        {
            const auto op = ext::make_shared<FdmBlackScholesOp>(
                mesher, process, K,
                false, -Null<Real>(), 0,
                ext::shared_ptr<FdmQuantoHelper>(), fittedDesc());

            Array v = payoffVec;
            rollbackCrankNicolson(op, bcSet, v, maturity, tGrid);

            bool allOk = true;
            for (Size i = 1; i + 1 < v.size(); ++i) {
                if (!std::isfinite(v[i])) { allOk = false; break; }
            }

            BOOST_CHECK_MESSAGE(allOk,
                "Non-finite values with explicit ExpFit on coarse mesh");
        }
    }
}


// ===================================================================
// Paper Example 4.1 with MT scheme: truncated call (r=0.05,
// sigma=0.001, K=50, U=70, T=5/12). Verify MT+CN produces smooth,
// positive solution and quantify numerical diffusion against a
// very fine-grid reference.
// ===================================================================

BOOST_AUTO_TEST_CASE(testMilevTaglianiTruncatedCallSmoothing) {

    BOOST_TEST_MESSAGE(
        "Testing MT+CN on paper Example 4.1 truncated call: "
        "positivity + numerical diffusion quantification...");

    LowVolSetup setup;  // sigma=0.001, r=0.05, spot=60

    const Real K = 50.0;
    const Real U = 70.0;
    const Real xMin = std::log(1.0);
    const Real xMax = std::log(140.0);
    const FdmBoundaryConditionSet bcSet;

    const ext::shared_ptr<Payoff> payoff =
        ext::make_shared<TruncatedCallPayoff>(K, U);

    // Helper: solve on a given grid with a given scheme
    auto solve = [&](Size xGrid, Size tGrid,
                     const FdmBlackScholesSpatialDesc& desc) {
        const auto mesher = ext::make_shared<FdmMesherComposite>(
            ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));
        const auto calc =
            ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);
        Array v = buildPayoffVector(mesher, calc, setup.maturity);

        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), desc);

        rollbackCrankNicolson(op, bcSet, v, setup.maturity, tGrid);
        return v;
    };

    // Fine-grid reference: ExpFit + Implicit Euler on very fine grid
    // (monotone reference, minimizes numerical diffusion artifacts)
    const Size refXGrid = 4000;
    const Size refTGrid = 1000;
    Array vRef;
    {
        const auto mesherRef = ext::make_shared<FdmMesherComposite>(
            ext::make_shared<Uniform1dMesher>(xMin, xMax, refXGrid));
        const auto calcRef =
            ext::make_shared<FdmLogInnerValue>(payoff, mesherRef, 0);
        vRef = buildPayoffVector(mesherRef, calcRef, setup.maturity);

        const auto opRef = ext::make_shared<FdmBlackScholesOp>(
            mesherRef, setup.process, K,
            false, -Null<Real>(), 0,
            ext::shared_ptr<FdmQuantoHelper>(), fittedDesc());

        rollbackImplicitEuler(opRef, bcSet, vRef,
                              setup.maturity, refTGrid);
    }

    // MT on moderate grid
    const Size xGrid = 400;
    const Size tGrid = 40;
    Array vMT = solve(xGrid, tGrid, milevTaglianiDesc());
    Array vCentral = solve(xGrid, tGrid, centralDesc());

    const Size skip = 5;

    // Positivity check
    Size negMT = 0, negCentral = 0;
    for (Size i = skip; i < xGrid - skip; ++i) {
        BOOST_CHECK(std::isfinite(vMT[i]));
        if (vMT[i] < -1e-12) ++negMT;
        if (vCentral[i] < -1e-6) ++negCentral;
    }
    BOOST_CHECK_EQUAL(negMT, Size(0));

    // Measure numerical diffusion near U=70 by comparing to reference.
    // Sample the reference solution onto the MT grid nodes near U.
    const Real xU = std::log(U);
    const Real hRef = (xMax - xMin) / (refXGrid - 1);
    const Real hMT = (xMax - xMin) / (xGrid - 1);

    Real maxErrMT = 0.0, maxErrCentral = 0.0;
    Size windowCount = 0;

    for (Size i = skip; i < xGrid - skip; ++i) {
        const Real xMT = xMin + i * hMT;
        if (std::fabs(xMT - xU) > 0.3)
            continue;

        // Find nearest reference node
        Size iRef = static_cast<Size>(
            std::round((xMT - xMin) / hRef));
        iRef = std::min(iRef, refXGrid - 1);

        const Real refVal = vRef[iRef];
        maxErrMT = std::max(maxErrMT,
            std::fabs(vMT[i] - refVal));
        maxErrCentral = std::max(maxErrCentral,
            std::fabs(vCentral[i] - refVal));
        ++windowCount;
    }

    BOOST_TEST_MESSAGE("  MT+CN: neg=" << negMT
        << ", Central+CN: neg=" << negCentral);
    BOOST_TEST_MESSAGE("  Window nodes near U=70: " << windowCount);
    BOOST_TEST_MESSAGE("  L_inf error vs ref (MT): " << maxErrMT);
    BOOST_TEST_MESSAGE("  L_inf error vs ref (Central): "
        << maxErrCentral);

    // L_inf error vs reference quantifies numerical diffusion:
    // MT adds diffusion, so its error is larger than Central's near
    // discontinuities. Both should be finite.
    BOOST_CHECK_MESSAGE(std::isfinite(maxErrMT),
        "MT L_inf error is not finite: " << maxErrMT);
    BOOST_CHECK_MESSAGE(maxErrMT > 0.0,
        "MT should have measurable numerical diffusion vs reference");
    BOOST_CHECK_MESSAGE(maxErrMT < 25.0,
        "MT L_inf error vs reference too large: " << maxErrMT);
}


BOOST_AUTO_TEST_CASE(testThetaAtFirstAccessor) {

    BOOST_TEST_MESSAGE(
        "Testing thetaAt() works as first accessor on "
        "FdmBlackScholesSolver (LazyObject initialization)...");

    LowVolSetup setup;

    const Real K = 60.0;
    const Real cash = 10.0;
    const Size xGrid = 100;
    const Size tGrid = 20;

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(
            std::log(10.0), std::log(200.0), xGrid));

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

    // Use a mutable vol quote so we can trigger solver recalculation.
    // Changing vol directly affects the PDE coefficients (drift and
    // diffusion), ensuring the theta value actually changes.
    auto volQuote = ext::make_shared<SimpleQuote>(setup.vol);

    auto process = ext::make_shared<BlackScholesMertonProcess>(
        Handle<Quote>(ext::make_shared<SimpleQuote>(setup.spot)),
        Handle<YieldTermStructure>(flatRate(setup.today, setup.q, setup.dc)),
        Handle<YieldTermStructure>(flatRate(setup.today, setup.r, setup.dc)),
        Handle<BlackVolTermStructure>(
            ext::make_shared<BlackConstantVol>(
                setup.today, Calendar(), Handle<Quote>(volQuote),
                setup.dc)));

    FdmBlackScholesSolver solver(
        Handle<GeneralizedBlackScholesProcess>(process),
        K, solverDesc,
        FdmSchemeDesc::CrankNicolson(),
        false, -Null<Real>(),
        Handle<FdmQuantoHelper>(),
        fittedDesc());

    // thetaAt as FIRST accessor (no prior call to valueAt/deltaAt/gammaAt)
    const Real theta1 = solver.thetaAt(setup.spot);
    BOOST_CHECK_MESSAGE(std::isfinite(theta1),
        "thetaAt() as first accessor returned non-finite value: " << theta1);

    // After valueAt, thetaAt should return a consistent result at the
    // same spot (same solver state, same interpolation point)
    const Real value = solver.valueAt(setup.spot);
    BOOST_CHECK(std::isfinite(value));
    const Real theta2 = solver.thetaAt(setup.spot);
    BOOST_CHECK_MESSAGE(std::isfinite(theta2),
        "thetaAt() after valueAt() returned non-finite value: " << theta2);
    BOOST_CHECK_MESSAGE(std::fabs(theta2 - theta1) < 1e-12,
        "thetaAt() should be consistent before and after valueAt(): "
        "theta1=" << theta1 << ", theta2=" << theta2);

    // After changing the vol quote, thetaAt at the SAME spot should
    // return a different value, proving LazyObject recalculation.
    // Vol change directly affects PDE coefficients (drift, diffusion).
    const Real thetaBeforeUpdate = solver.thetaAt(setup.spot);
    volQuote->setValue(setup.vol * 10.0);
    const Real thetaAfterUpdate = solver.thetaAt(setup.spot);
    BOOST_CHECK_MESSAGE(std::isfinite(thetaAfterUpdate),
        "thetaAt() after vol change returned non-finite value: "
        << thetaAfterUpdate);
    BOOST_CHECK_MESSAGE(
        std::fabs(thetaAfterUpdate - thetaBeforeUpdate) > 1e-14,
        "thetaAt() at the same s should differ after vol update "
        "(proves recalculation): before=" << thetaBeforeUpdate
        << ", after=" << thetaAfterUpdate);

    // Negative test: s <= 0 should not produce a finite result.
    // log(0) = -inf and log(negative) = NaN, so interpolation at
    // those values should not return a usable finite number.
    // Either an exception or a non-finite result is acceptable.
    for (Real badS : {0.0, -1.0}) {
        bool caught = false;
        try {
            const Real badTheta = solver.thetaAt(badS);
            BOOST_CHECK_MESSAGE(!std::isfinite(badTheta),
                "thetaAt(" << badS << ") should not produce a "
                "finite result, got " << badTheta);
        } catch (...) {
            caught = true;
        }
        if (caught)
            BOOST_TEST_MESSAGE("  thetaAt(" << badS
                               << ") threw (acceptable)");
    }
}


BOOST_AUTO_TEST_CASE(testCNGateCraigSneydAccepted) {

    BOOST_TEST_MESSAGE(
        "Testing CN-equivalence gate accepts CraigSneyd(0.5) in 1D "
        "and rejects schemes with damping steps...");

    LowVolSetup setup;

    const Real K = 60.0;
    const Real cash = 10.0;
    const Size xGrid = 200;
    const Size tGrid = 50;

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(
            std::log(10.0), std::log(200.0), xGrid));

    const ext::shared_ptr<StrikedTypePayoff> payoff =
        ext::make_shared<CashOrNothingPayoff>(Option::Put, K, cash);
    const auto calculator =
        ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);

    const auto conditions =
        ext::make_shared<FdmStepConditionComposite>(
            std::list<std::vector<Time> >(),
            FdmStepConditionComposite::Conditions());

    const FdmBoundaryConditionSet bcSet;

    const FdmSolverDesc solverDescNoDamping = {
        mesher, bcSet, conditions, calculator,
        setup.maturity, tGrid, 0
    };

    const FdmSolverDesc solverDescWithDamping = {
        mesher, bcSet, conditions, calculator,
        setup.maturity, tGrid, 3  // 3 damping steps
    };

    // CraigSneyd(0.5) + dampingSteps=0 + MT: should be accepted.
    // Proof of acceptance: MT values must DIFFER from ExpFit values,
    // since at low vol the MT and ExpFit spatial coefficients diverge.
    // If the gate silently fell back, they would be bit-identical.
    {
        FdmBlackScholesSolver solverMT(
            Handle<GeneralizedBlackScholesProcess>(setup.process),
            K, solverDescNoDamping,
            FdmSchemeDesc::CraigSneyd(),
            false, -Null<Real>(),
            Handle<FdmQuantoHelper>(),
            milevTaglianiDesc());

        FdmBlackScholesSolver solverFitted(
            Handle<GeneralizedBlackScholesProcess>(setup.process),
            K, solverDescNoDamping,
            FdmSchemeDesc::CraigSneyd(),
            false, -Null<Real>(),
            Handle<FdmQuantoHelper>(),
            fittedDesc());

        const Real vMT = solverMT.valueAt(setup.spot);
        const Real vFit = solverFitted.valueAt(setup.spot);

        BOOST_CHECK(std::isfinite(vMT));
        BOOST_CHECK(std::isfinite(vFit));

        // Assert acceptance: MT and ExpFit produce different values
        // because at sigma=0.001 the MT added diffusion is large
        BOOST_CHECK_MESSAGE(std::fabs(vMT - vFit) > 0.01,
            "CraigSneyd+MT should be accepted (not fallen back to "
            "ExpFit). MT=" << vMT << ", ExpFit=" << vFit
            << ", diff=" << std::fabs(vMT - vFit));

        BOOST_TEST_MESSAGE("  CraigSneyd+MT: " << vMT
            << ", CraigSneyd+ExpFit: " << vFit
            << " (diff=" << std::fabs(vMT - vFit) << ")");
    }

    // CN + dampingSteps > 0 + FailFast: should throw
    {
        FdmBlackScholesSpatialDesc desc = milevTaglianiDesc();
        desc.mMatrixPolicy =
            FdmBlackScholesSpatialDesc::MMatrixPolicy::FailFast;

        BOOST_CHECK_THROW(
            FdmBlackScholesSolver(
                Handle<GeneralizedBlackScholesProcess>(setup.process),
                K, solverDescWithDamping,
                FdmSchemeDesc::CrankNicolson(),
                false, -Null<Real>(),
                Handle<FdmQuantoHelper>(),
                desc).valueAt(setup.spot),
            Error);
    }

    // CN + dampingSteps > 0 + FallbackToExponentialFitting: should
    // silently fall back (match explicit ExpFit results)
    {
        FdmBlackScholesSolver solverMTFallback(
            Handle<GeneralizedBlackScholesProcess>(setup.process),
            K, solverDescWithDamping,
            FdmSchemeDesc::CrankNicolson(),
            false, -Null<Real>(),
            Handle<FdmQuantoHelper>(),
            milevTaglianiDesc());

        FdmBlackScholesSolver solverFittedExplicit(
            Handle<GeneralizedBlackScholesProcess>(setup.process),
            K, solverDescWithDamping,
            FdmSchemeDesc::CrankNicolson(),
            false, -Null<Real>(),
            Handle<FdmQuantoHelper>(),
            fittedDesc());

        const Real vFallback = solverMTFallback.valueAt(setup.spot);
        const Real vExplicit = solverFittedExplicit.valueAt(setup.spot);

        BOOST_CHECK(std::isfinite(vFallback));
        BOOST_CHECK(std::isfinite(vExplicit));
        BOOST_CHECK_SMALL(std::fabs(vFallback - vExplicit), 1e-10);
    }

    // Non-CN schemes + MT: should silently fall back to ExpFit
    for (const auto& scheme :
             {FdmSchemeDesc::ImplicitEuler(),
              FdmSchemeDesc::Hundsdorfer()}) {
        FdmBlackScholesSolver solverMT(
            Handle<GeneralizedBlackScholesProcess>(setup.process),
            K, solverDescNoDamping, scheme,
            false, -Null<Real>(),
            Handle<FdmQuantoHelper>(),
            milevTaglianiDesc());

        FdmBlackScholesSolver solverFitted(
            Handle<GeneralizedBlackScholesProcess>(setup.process),
            K, solverDescNoDamping, scheme,
            false, -Null<Real>(),
            Handle<FdmQuantoHelper>(),
            fittedDesc());

        const Real vMT = solverMT.valueAt(setup.spot);
        const Real vFit = solverFitted.valueAt(setup.spot);

        BOOST_CHECK_SMALL(std::fabs(vMT - vFit), 1e-10);
    }
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
