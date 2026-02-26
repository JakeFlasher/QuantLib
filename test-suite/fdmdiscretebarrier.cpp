// r6
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Tests for FdmDiscreteBarrierStepCondition and the paper-faithful
 demonstration that repeated discrete monitoring injections accumulate
 Crank-Nicolson oscillations under standard central differencing, while
 exponential fitting (Scheme 1) and the Milev-Tagliani CN effective
 diffusion (Scheme 2) suppress them.
*/

#include "toplevelfixture.hpp"
#include "utilities.hpp"

#include <ql/methods/finitedifferences/meshers/uniform1dmesher.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmeshercomposite.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesop.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmdiscretebarrierstepcondition.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/utilities/null.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <cmath>
#include <list>
#include <vector>

using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FdmDiscreteBarrierTests)

namespace {

    // ---------------------------------------------------------------
    //  Market-data fixture used by the stress test
    // ---------------------------------------------------------------
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
          spot(100.0), r(0.05), q(0.0), vol(0.001),
          maturity(0.5) {

            Settings::instance().evaluationDate() = today;
            process = ext::make_shared<BlackScholesMertonProcess>(
                Handle<Quote>(ext::make_shared<SimpleQuote>(spot)),
                Handle<YieldTermStructure>(flatRate(today, q, dc)),
                Handle<YieldTermStructure>(flatRate(today, r, dc)),
                Handle<BlackVolTermStructure>(flatVol(today, vol, dc)));
        }
    };

    // ---------------------------------------------------------------
    //  Monitoring-date generator (excludes maturity)
    // ---------------------------------------------------------------
    std::vector<Time> equallySpacedMonitoringTimes(
            Time maturity, Size nIntervals) {
        std::vector<Time> t;
        t.reserve(nIntervals - 1);
        for (Size i = 1; i < nIntervals; ++i)
            t.push_back(maturity * Real(i) / Real(nIntervals));
        return t;
    }

    // ---------------------------------------------------------------
    //  Pre-projected payoff: max(S-K,0) * 1_{[L,U]}(S) in log-space
    // ---------------------------------------------------------------
    Array buildProjectedCallPayoff(
            const ext::shared_ptr<FdmMesher>& mesher,
            Real K, Real L, Real U, Size direction = 0) {
        Array rhs(mesher->layout()->size());
        for (const auto& iter : *mesher->layout()) {
            const Real x = mesher->location(iter, direction);
            const Real s = std::exp(x);
            Real v = std::max(s - K, 0.0);
            if (s < L || s > U) v = 0.0;
            rhs[iter.index()] = v;
        }
        return rhs;
    }

    // ---------------------------------------------------------------
    //  Per-monitoring-date statistics collected by the probe
    // ---------------------------------------------------------------
    struct MonitorStats {
        Time t = 0.0;
        Real minVal = 0.0;
        Real maxVal = 0.0;
        Size negCount = 0;
        Size overCount = 0;
    };

    // ---------------------------------------------------------------
    //  Test-only wrapper: records bound-violation diagnostics at
    //  monitoring times, then delegates to the real barrier condition.
    //
    //  "before" = state after CN rollback, prior to barrier projection
    //  "after"  = state after the projection has zeroed exterior nodes
    // ---------------------------------------------------------------
    class BarrierMonitorProbe : public StepCondition<Array> {
      public:
        BarrierMonitorProbe(
            const ext::shared_ptr<FdmMesher>& mesher,
            const std::vector<Time>& monitoringTimes,
            Real L, Real U, Real maxPayoff,
            Size direction = 0,
            Real negTol = 1e-10,
            Real overTol = 1e-8)
        : barrier_(ext::make_shared<FdmDiscreteBarrierStepCondition>(
              mesher, monitoringTimes, L, U, direction)),
          monitoringTimes_(barrier_->monitoringTimes()),
          maxPayoff_(maxPayoff),
          negTol_(negTol), overTol_(overTol) {}

        void applyTo(Array& a, Time t) const override {
            const bool isMon = std::binary_search(
                monitoringTimes_.begin(), monitoringTimes_.end(), t);

            if (isMon)
                snapshot(a, t, before_);

            barrier_->applyTo(a, t);

            if (isMon)
                snapshot(a, t, after_);
        }

        const std::vector<MonitorStats>& before() const { return before_; }
        const std::vector<MonitorStats>& after()  const { return after_; }

      private:
        void snapshot(const Array& a, Time t,
                      std::vector<MonitorStats>& out) const {
            MonitorStats s;
            s.t = t;
            if (!a.empty()) {
                s.minVal = *std::min_element(a.begin(), a.end());
                s.maxVal = *std::max_element(a.begin(), a.end());
            }
            for (Size i = 0; i < a.size(); ++i) {
                if (a[i] < -negTol_)  ++s.negCount;
                if (a[i] > maxPayoff_ + overTol_) ++s.overCount;
            }
            out.push_back(s);
        }

        ext::shared_ptr<FdmDiscreteBarrierStepCondition> barrier_;
        std::vector<Time> monitoringTimes_;
        Real maxPayoff_;
        Real negTol_, overTol_;

        mutable std::vector<MonitorStats> before_;
        mutable std::vector<MonitorStats> after_;
    };

    // ---------------------------------------------------------------
    //  Spatial-descriptor helpers
    // ---------------------------------------------------------------
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

    // ---------------------------------------------------------------
    //  Run a CN rollback with barrier monitoring and return the probe
    // ---------------------------------------------------------------
    ext::shared_ptr<BarrierMonitorProbe> runBarrierRollback(
            const ext::shared_ptr<FdmMesher>& mesher,
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
            const Array& payoffVec,
            const std::vector<Time>& monTimes,
            Real K, Real L, Real U, Real maxPayoff,
            Time maturity, Size tGrid,
            const FdmBlackScholesSpatialDesc& spatialDesc,
            Array& result) {

        const auto op = ext::make_shared<FdmBlackScholesOp>(
            mesher, process, K,
            false, -Null<Real>(), Size(0),
            ext::shared_ptr<FdmQuantoHelper>(), spatialDesc);

        const auto probe = ext::make_shared<BarrierMonitorProbe>(
            mesher, monTimes, L, U, maxPayoff);

        FdmStepConditionComposite::Conditions condList;
        condList.push_back(probe);

        const auto conditions = ext::make_shared<FdmStepConditionComposite>(
            std::list<std::vector<Time> >(1, monTimes),
            condList);

        result = payoffVec;
        FdmBackwardSolver(op, FdmBoundaryConditionSet(),
                          conditions, FdmSchemeDesc::CrankNicolson())
            .rollback(result, maturity, 0.0, tGrid, 0);

        return probe;
    }

} // anonymous namespace


// ===================================================================
//  Unit test: basic step-condition mechanics
// ===================================================================

BOOST_AUTO_TEST_CASE(testFdmDiscreteBarrierStepConditionBasics) {

    BOOST_TEST_MESSAGE(
        "Testing FdmDiscreteBarrierStepCondition basic operations...");

    const Size xGrid = 100;
    const Real xMin = std::log(50.0);
    const Real xMax = std::log(200.0);
    const Real L = 80.0;
    const Real U = 120.0;
    const Real xL = std::log(L);
    const Real xU = std::log(U);
    const std::vector<Time> monTimes = {0.25, 0.5};

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    FdmDiscreteBarrierStepCondition cond(mesher, monTimes, L, U);

    // Accessors
    BOOST_CHECK_EQUAL(cond.monitoringTimes().size(), Size(2));
    BOOST_CHECK_EQUAL(cond.lowerBarrier(), L);
    BOOST_CHECK_EQUAL(cond.upperBarrier(), U);

    // ---- Non-monitoring time: array must be unchanged ----
    Array a(mesher->layout()->size(), 1.0);
    Array aCopy(a);

    cond.applyTo(a, 0.1);    // t=0.1 is NOT a monitoring time

    for (Size i = 0; i < a.size(); ++i)
        BOOST_CHECK_EQUAL(a[i], aCopy[i]);

    // ---- Monitoring time: exterior nodes zeroed, interior untouched ----
    cond.applyTo(a, 0.25);   // t=0.25 IS a monitoring time

    for (const auto& iter : *mesher->layout()) {
        const Real x = mesher->location(iter, 0);
        if (x < xL || x > xU) {
            BOOST_CHECK_EQUAL(a[iter.index()], 0.0);
        } else {
            BOOST_CHECK_EQUAL(a[iter.index()], 1.0);
        }
    }

    // ---- Second monitoring time works too ----
    // Reset interior to 2.0, leave exterior at 0.0
    for (const auto& iter : *mesher->layout()) {
        const Real x = mesher->location(iter, 0);
        if (!(x < xL || x > xU))
            a[iter.index()] = 2.0;
    }

    cond.applyTo(a, 0.5);

    for (const auto& iter : *mesher->layout()) {
        const Real x = mesher->location(iter, 0);
        if (x < xL || x > xU) {
            BOOST_CHECK_EQUAL(a[iter.index()], 0.0);
        } else {
            BOOST_CHECK_EQUAL(a[iter.index()], 2.0);
        }
    }

    // ---- Duplicate / unsorted times are handled correctly ----
    FdmDiscreteBarrierStepCondition cond2(
        mesher, std::vector<Time>({0.5, 0.25, 0.25}), L, U);
    BOOST_CHECK_EQUAL(cond2.monitoringTimes().size(), Size(2));
}


// ===================================================================
//  Paper-faithful stress test: repeated discrete barrier monitoring
//  with CN + three spatial schemes
// ===================================================================

BOOST_AUTO_TEST_CASE(testDiscreteBarrierOscillationAccumulation) {

    BOOST_TEST_MESSAGE(
        "Testing discrete barrier monitoring: CN + StandardCentral "
        "accumulates oscillations across monitoring dates while "
        "ExponentialFitting and MilevTagliani remain stable...");

    LowVolSetup setup;

    // -- Contract parameters --
    const Real K = 100.0;
    const Real L = 95.0;
    const Real U = 110.0;
    const Real maxPayoff = U - K;   // = 10

    // 12 sub-intervals → 11 interior monitoring dates, all < maturity
    const Size nIntervals = 12;
    const std::vector<Time> monTimes =
        equallySpacedMonitoringTimes(setup.maturity, nIntervals);

    // -- Grid --
    const Size xGrid = 400;
    const Size tGrid = 200;
    const Real xMin = std::log(50.0);
    const Real xMax = std::log(200.0);

    const auto mesher = ext::make_shared<FdmMesherComposite>(
        ext::make_shared<Uniform1dMesher>(xMin, xMax, xGrid));

    // Pre-projected payoff: max(S-K,0) · 1_{[L,U]}(S)
    const Array payoffVec =
        buildProjectedCallPayoff(mesher, K, L, U);

    // ===========================================================
    //  Case A: StandardCentral + CN — expect oscillation accumulation
    // ===========================================================
    {
        Array result;
        const auto probe = runBarrierRollback(
            mesher, setup.process, payoffVec, monTimes,
            K, L, U, maxPayoff,
            setup.maturity, tGrid,
            centralDesc(), result);

        // "after" stats: post-projection; exterior nodes are zeroed,
        // so any remaining violations are interior oscillation artifacts.
        Size violatingDates = 0;
        for (const auto& s : probe->after()) {
            if (s.negCount > 0 || s.overCount > 0)
                ++violatingDates;
        }

        BOOST_TEST_MESSAGE(
            "  StandardCentral+CN: "
            << violatingDates << " / " << probe->after().size()
            << " monitoring dates with post-projection violations");

        // More than one monitoring date with violations proves
        // that oscillations from one discontinuity injection persist
        // and compound across subsequent monitoring intervals.
        BOOST_CHECK_MESSAGE(violatingDates > 1,
            "Expected oscillation accumulation across >1 monitoring "
            "dates for StandardCentral + CN, but only found "
            << violatingDates);
    }

    // ===========================================================
    //  Case B: ExponentialFitting (Scheme 1) + CN — expect stability
    // ===========================================================
    {
        Array result;
        const auto probe = runBarrierRollback(
            mesher, setup.process, payoffVec, monTimes,
            K, L, U, maxPayoff,
            setup.maturity, tGrid,
            fittedDesc(), result);

        BOOST_TEST_MESSAGE(
            "  ExponentialFitting+CN: checking "
            << probe->after().size() << " monitoring dates");

        for (Size d = 0; d < probe->after().size(); ++d) {
            const auto& s = probe->after()[d];
            BOOST_CHECK_MESSAGE(s.negCount == 0,
                "ExponentialFitting: negativity at monitoring date "
                << d << " (t=" << s.t
                << "), count=" << s.negCount
                << ", min=" << s.minVal);
            BOOST_CHECK_MESSAGE(s.overCount == 0,
                "ExponentialFitting: overshoot at monitoring date "
                << d << " (t=" << s.t
                << "), count=" << s.overCount
                << ", max=" << s.maxVal);
        }

        // Final-time sanity: all values must be finite
        for (Size i = 0; i < result.size(); ++i)
            BOOST_CHECK(std::isfinite(result[i]));
    }

    // ===========================================================
    //  Case C: MilevTagliani CN effective diffusion (Scheme 2) + CN
    // ===========================================================
    {
        Array result;
        const auto probe = runBarrierRollback(
            mesher, setup.process, payoffVec, monTimes,
            K, L, U, maxPayoff,
            setup.maturity, tGrid,
            milevTaglianiDesc(), result);

        BOOST_TEST_MESSAGE(
            "  MilevTagliani+CN: checking "
            << probe->after().size() << " monitoring dates");

        for (Size d = 0; d < probe->after().size(); ++d) {
            const auto& s = probe->after()[d];
            BOOST_CHECK_MESSAGE(s.negCount == 0,
                "MilevTagliani: negativity at monitoring date "
                << d << " (t=" << s.t
                << "), count=" << s.negCount
                << ", min=" << s.minVal);
            BOOST_CHECK_MESSAGE(s.overCount == 0,
                "MilevTagliani: overshoot at monitoring date "
                << d << " (t=" << s.t
                << "), count=" << s.overCount
                << ", max=" << s.maxVal);
        }

        for (Size i = 0; i < result.size(); ++i)
            BOOST_CHECK(std::isfinite(result[i]));
    }
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
