// r6
#include <ql/exercise.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmeshercomposite.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmdiscretebarrierstepcondition.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp>
#include <ql/methods/finitedifferences/utilities/fdmdirichletboundary.hpp>
#include <ql/methods/finitedifferences/utilities/fdmdividendhandler.hpp>
#include <ql/methods/finitedifferences/utilities/fdminnervaluecalculator.hpp>
#include <ql/pricingengines/barrier/fdblackscholesbarrierengine.hpp>
#include <ql/pricingengines/barrier/fdblackscholesrebateengine.hpp>
#include <ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp>
#include <utility>

namespace QuantLib {

    // ---- Continuous-monitoring constructors (original) ----

    FdBlackScholesBarrierEngine::FdBlackScholesBarrierEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        Size tGrid, Size xGrid, Size dampingSteps,
        const FdmSchemeDesc& schemeDesc,
        bool localVol, Real illegalLocalVolOverwrite,
        FdmBlackScholesSpatialDesc spatialDesc)
    : process_(std::move(process)),
      tGrid_(tGrid), xGrid_(xGrid), dampingSteps_(dampingSteps),
      schemeDesc_(schemeDesc), localVol_(localVol),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      spatialDesc_(spatialDesc) {
        registerWith(process_);
    }

    FdBlackScholesBarrierEngine::FdBlackScholesBarrierEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        DividendSchedule dividends,
        Size tGrid, Size xGrid, Size dampingSteps,
        const FdmSchemeDesc& schemeDesc,
        bool localVol, Real illegalLocalVolOverwrite,
        FdmBlackScholesSpatialDesc spatialDesc)
    : process_(std::move(process)), dividends_(std::move(dividends)),
      tGrid_(tGrid), xGrid_(xGrid), dampingSteps_(dampingSteps),
      schemeDesc_(schemeDesc), localVol_(localVol),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      spatialDesc_(spatialDesc) {
        registerWith(process_);
    }

    // ---- Discrete-monitoring constructors ----

    FdBlackScholesBarrierEngine::FdBlackScholesBarrierEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        std::vector<Date> monitoringDates,
        Size tGrid, Size xGrid, Size dampingSteps,
        const FdmSchemeDesc& schemeDesc,
        bool localVol, Real illegalLocalVolOverwrite,
        FdmBlackScholesSpatialDesc spatialDesc)
    : process_(std::move(process)),
      monitoringDates_(std::move(monitoringDates)),
      tGrid_(tGrid), xGrid_(xGrid), dampingSteps_(dampingSteps),
      schemeDesc_(schemeDesc), localVol_(localVol),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      spatialDesc_(spatialDesc) {
        registerWith(process_);
    }

    FdBlackScholesBarrierEngine::FdBlackScholesBarrierEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        DividendSchedule dividends,
        std::vector<Date> monitoringDates,
        Size tGrid, Size xGrid, Size dampingSteps,
        const FdmSchemeDesc& schemeDesc,
        bool localVol, Real illegalLocalVolOverwrite,
        FdmBlackScholesSpatialDesc spatialDesc)
    : process_(std::move(process)), dividends_(std::move(dividends)),
      monitoringDates_(std::move(monitoringDates)),
      tGrid_(tGrid), xGrid_(xGrid), dampingSteps_(dampingSteps),
      schemeDesc_(schemeDesc), localVol_(localVol),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      spatialDesc_(spatialDesc) {
        registerWith(process_);
    }

    // ---- Top-level dispatch ----

    void FdBlackScholesBarrierEngine::calculate() const {
        if (monitoringDates_.empty())
            calculateContinuous();
        else
            calculateDiscrete();
    }

    // ---- Continuous path (original logic, verbatim) ----

    void FdBlackScholesBarrierEngine::calculateContinuous() const {

        const ext::shared_ptr<StrikedTypePayoff> payoff =
            ext::dynamic_pointer_cast<StrikedTypePayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-striked type payoff given");
        QL_REQUIRE(payoff->strike() > 0.0, "strike must be positive");
        QL_REQUIRE(arguments_.exercise->type() == Exercise::European,
                   "only european style option are supported");

        const auto spot = process_->x0();
        QL_REQUIRE(spot > 0.0, "negative or null underlying given");
        QL_REQUIRE(!triggered(spot), "barrier touched");

        const Time maturity =
            process_->time(arguments_.exercise->lastDate());

        Real xMin = Null<Real>();
        Real xMax = Null<Real>();
        if (   arguments_.barrierType == Barrier::DownIn
            || arguments_.barrierType == Barrier::DownOut)
            xMin = std::log(arguments_.barrier);
        if (   arguments_.barrierType == Barrier::UpIn
            || arguments_.barrierType == Barrier::UpOut)
            xMax = std::log(arguments_.barrier);

        const ext::shared_ptr<Fdm1dMesher> equityMesher(
            new FdmBlackScholesMesher(
                xGrid_, process_, maturity, payoff->strike(),
                xMin, xMax, 0.0001, 1.5,
                std::make_pair(Null<Real>(), Null<Real>()),
                dividends_));

        const ext::shared_ptr<FdmMesher> mesher(
            ext::make_shared<FdmMesherComposite>(equityMesher));

        ext::shared_ptr<FdmInnerValueCalculator> calculator(
            ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0));

        std::list<ext::shared_ptr<StepCondition<Array>>> stepConditions;
        std::list<std::vector<Time>> stoppingTimes;

        ext::shared_ptr<FdmDividendHandler> dividendCondition(
            ext::make_shared<FdmDividendHandler>(
                dividends_, mesher,
                process_->riskFreeRate()->referenceDate(),
                process_->riskFreeRate()->dayCounter(), 0));

        if (!dividends_.empty()) {
            stepConditions.push_back(dividendCondition);
            std::vector<Time> dividendTimes =
                dividendCondition->dividendTimes();
            for (auto& t : dividendTimes)
                t = std::min(maturity, t);
            stoppingTimes.push_back(dividendTimes);
        }

        ext::shared_ptr<FdmStepConditionComposite> conditions(
            ext::make_shared<FdmStepConditionComposite>(
                stoppingTimes, stepConditions));

        FdmBoundaryConditionSet boundaries;
        if (   arguments_.barrierType == Barrier::DownIn
            || arguments_.barrierType == Barrier::DownOut)
            boundaries.push_back(
                ext::make_shared<FdmDirichletBoundary>(
                    mesher, arguments_.rebate, 0,
                    FdmDirichletBoundary::Lower));
        if (   arguments_.barrierType == Barrier::UpIn
            || arguments_.barrierType == Barrier::UpOut)
            boundaries.push_back(
                ext::make_shared<FdmDirichletBoundary>(
                    mesher, arguments_.rebate, 0,
                    FdmDirichletBoundary::Upper));

        FdmSolverDesc solverDesc = {
            mesher, boundaries, conditions, calculator,
            maturity, tGrid_, dampingSteps_
        };

        ext::shared_ptr<FdmBlackScholesSolver> solver(
            ext::make_shared<FdmBlackScholesSolver>(
                Handle<GeneralizedBlackScholesProcess>(process_),
                payoff->strike(), solverDesc, schemeDesc_,
                localVol_, illegalLocalVolOverwrite_,
                Handle<FdmQuantoHelper>(),
                spatialDesc_));

        results_.value = solver->valueAt(spot);
        results_.delta = solver->deltaAt(spot);
        results_.gamma = solver->gammaAt(spot);
        results_.theta = solver->thetaAt(spot);

        if (   arguments_.barrierType == Barrier::DownIn
            || arguments_.barrierType == Barrier::UpIn) {

            ext::shared_ptr<StrikedTypePayoff> payoff2 =
                ext::dynamic_pointer_cast<StrikedTypePayoff>(
                    arguments_.payoff);
            VanillaOption vanillaOption(payoff2, arguments_.exercise);
            vanillaOption.setPricingEngine(
                ext::make_shared<FdBlackScholesVanillaEngine>(
                    process_, dividends_, tGrid_, xGrid_,
                    0, schemeDesc_, localVol_,
                    illegalLocalVolOverwrite_,
                    FdBlackScholesVanillaEngine::Spot,
                    spatialDesc_));

            BarrierOption rebateOption(
                arguments_.barrierType, arguments_.barrier,
                arguments_.rebate, payoff2, arguments_.exercise);

            const Size min_grid_size = 50;
            const Size rebateDampingSteps =
                (dampingSteps_ > 0)
                    ? std::min(Size(1), dampingSteps_ / 2) : 0;

            rebateOption.setPricingEngine(
                ext::make_shared<FdBlackScholesRebateEngine>(
                    process_, dividends_, tGrid_,
                    std::max(min_grid_size, xGrid_ / 5),
                    rebateDampingSteps, schemeDesc_, localVol_,
                    illegalLocalVolOverwrite_, spatialDesc_));

            results_.value =
                vanillaOption.NPV()   + rebateOption.NPV()
                - results_.value;
            results_.delta =
                vanillaOption.delta() + rebateOption.delta()
                - results_.delta;
            results_.gamma =
                vanillaOption.gamma() + rebateOption.gamma()
                - results_.gamma;
            results_.theta =
                vanillaOption.theta() + rebateOption.theta()
                - results_.theta;
        }
    }

    // ---- Discrete-monitoring path ----

    void FdBlackScholesBarrierEngine::calculateDiscrete() const {

        const ext::shared_ptr<StrikedTypePayoff> payoff =
            ext::dynamic_pointer_cast<StrikedTypePayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-striked type payoff given");
        QL_REQUIRE(payoff->strike() > 0.0, "strike must be positive");
        QL_REQUIRE(arguments_.exercise->type() == Exercise::European,
                   "only european style options are supported "
                   "for discrete-monitoring barriers");

        const bool isKnockIn =
            (arguments_.barrierType == Barrier::DownIn
             || arguments_.barrierType == Barrier::UpIn);

        if (isKnockIn && arguments_.rebate != 0.0) {
            QL_FAIL("Non-zero rebate is not supported for "
                    "discrete-monitoring knock-in barriers");
        }

        const auto spot = process_->x0();
        QL_REQUIRE(spot > 0.0, "negative or null underlying given");

        const Time maturity =
            process_->time(arguments_.exercise->lastDate());

        // ---- Multi-point mesher (grid extends beyond barrier) ----
        std::vector<std::tuple<Real, Real, bool>> cPoints;
        cPoints.emplace_back(payoff->strike(), 0.1, true);

        if (   arguments_.barrierType == Barrier::DownIn
            || arguments_.barrierType == Barrier::DownOut)
            cPoints.emplace_back(arguments_.barrier, 0.1, true);
        if (   arguments_.barrierType == Barrier::UpIn
            || arguments_.barrierType == Barrier::UpOut)
            cPoints.emplace_back(arguments_.barrier, 0.1, true);

        const ext::shared_ptr<Fdm1dMesher> equityMesher =
            ext::make_shared<FdmBlackScholesMesher>(
                xGrid_, process_, maturity, payoff->strike(),
                Null<Real>(), Null<Real>(), 0.0001, 1.5,
                cPoints, dividends_);

        const ext::shared_ptr<FdmMesher> mesher =
            ext::make_shared<FdmMesherComposite>(equityMesher);

        const ext::shared_ptr<FdmInnerValueCalculator> calculator =
            ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);

        // ---- Build step conditions ----
        std::list<std::vector<Time>> stoppingTimeLists;
        std::list<ext::shared_ptr<StepCondition<Array>>> stepConditionsList;

        // Dividends
        if (!dividends_.empty()) {
            auto dividendCondition =
                ext::make_shared<FdmDividendHandler>(
                    dividends_, mesher,
                    process_->riskFreeRate()->referenceDate(),
                    process_->riskFreeRate()->dayCounter(), 0);
            stepConditionsList.push_back(dividendCondition);
            std::vector<Time> divTimes =
                dividendCondition->dividendTimes();
            for (auto& t : divTimes)
                t = std::min(maturity, t);
            stoppingTimeLists.push_back(divTimes);
        }

        // Barrier corridor for a single-barrier instrument
        Real lowerBarrier, upperBarrier;
        if (   arguments_.barrierType == Barrier::DownIn
            || arguments_.barrierType == Barrier::DownOut) {
            lowerBarrier = arguments_.barrier;
            upperBarrier = 1e15;
        } else {
            lowerBarrier = 1e-15;
            upperBarrier = arguments_.barrier;
        }

        // Check for t=0 monitoring (valuation-date barrier check).
        // The paper multiplies the initial condition by 1_[L,U](S),
        // meaning spot outside the corridor at t=0 is an immediate
        // knock-out.
        bool knockedOutAtT0 = false;
        for (const auto& d : monitoringDates_) {
            const Time t = process_->time(d);
            if (std::fabs(t) < 1e-12) {
                if (spot < lowerBarrier || spot > upperBarrier)
                    knockedOutAtT0 = true;
                break;
            }
        }

        if (knockedOutAtT0) {
            if (isKnockIn) {
                // Knock-in: if knocked out at t=0, the knock-in is
                // immediately active → price equals vanilla
                VanillaOption vanillaOption(payoff, arguments_.exercise);
                vanillaOption.setPricingEngine(
                    ext::make_shared<FdBlackScholesVanillaEngine>(
                        process_, dividends_, tGrid_, xGrid_,
                        0, schemeDesc_, localVol_,
                        illegalLocalVolOverwrite_,
                        FdBlackScholesVanillaEngine::Spot,
                        spatialDesc_));
                results_.value = vanillaOption.NPV();
                results_.delta = vanillaOption.delta();
                results_.gamma = vanillaOption.gamma();
                results_.theta = vanillaOption.theta();
            } else {
                // Knock-out: immediately dead → rebate paid today.
                // The step condition uses hit-time rebate semantics
                // (immediate payment), so the value is just the rebate
                // with no further discounting.
                results_.value = arguments_.rebate;
                results_.delta = 0.0;
                results_.gamma = 0.0;
                results_.theta = 0.0;
            }
            return;
        }

        // Monitoring times (t>0 only; t=0 handled above)
        std::vector<Time> monTimes;
        for (const auto& d : monitoringDates_) {
            const Time t = process_->time(d);
            if (t > 0.0 && t <= maturity)
                monTimes.push_back(t);
        }

        auto barrierCondition =
            ext::make_shared<FdmDiscreteBarrierStepCondition>(
                mesher, monTimes,
                lowerBarrier, upperBarrier,
                arguments_.rebate, Size(0));
        stepConditionsList.push_back(barrierCondition);
        stoppingTimeLists.push_back(monTimes);

        auto conditions =
            ext::make_shared<FdmStepConditionComposite>(
                stoppingTimeLists, stepConditionsList);

        // No Dirichlet boundaries for discrete monitoring
        const FdmBoundaryConditionSet boundaries;

        FdmSolverDesc solverDesc = {
            mesher, boundaries, conditions, calculator,
            maturity, tGrid_, dampingSteps_
        };

        auto solver = ext::make_shared<FdmBlackScholesSolver>(
            Handle<GeneralizedBlackScholesProcess>(process_),
            payoff->strike(), solverDesc, schemeDesc_,
            localVol_, illegalLocalVolOverwrite_,
            Handle<FdmQuantoHelper>(),
            spatialDesc_);

        results_.value = solver->valueAt(spot);
        results_.delta = solver->deltaAt(spot);
        results_.gamma = solver->gammaAt(spot);
        results_.theta = solver->thetaAt(spot);

        // Handle knock-in via parity: In = Vanilla - Out
        if (isKnockIn) {
            VanillaOption vanillaOption(payoff, arguments_.exercise);
            vanillaOption.setPricingEngine(
                ext::make_shared<FdBlackScholesVanillaEngine>(
                    process_, dividends_, tGrid_, xGrid_,
                    0, schemeDesc_, localVol_,
                    illegalLocalVolOverwrite_,
                    FdBlackScholesVanillaEngine::Spot,
                    spatialDesc_));

            results_.value = vanillaOption.NPV()   - results_.value;
            results_.delta = vanillaOption.delta()  - results_.delta;
            results_.gamma = vanillaOption.gamma()  - results_.gamma;
            results_.theta = vanillaOption.theta()  - results_.theta;
        }
    }

}
