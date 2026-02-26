// r6
#include <ql/exercise.hpp>
#include <ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp>
#include <ql/methods/finitedifferences/utilities/escroweddividendadjustment.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmeshercomposite.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp>
#include <ql/methods/finitedifferences/utilities/fdminnervaluecalculator.hpp>
#include <ql/methods/finitedifferences/utilities/fdmescrowedloginnervaluecalculator.hpp>
#include <ql/methods/finitedifferences/utilities/fdmquantohelper.hpp>
#include <ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp>
#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {

    FdBlackScholesVanillaEngine::FdBlackScholesVanillaEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        Size tGrid, Size xGrid, Size dampingSteps,
        const FdmSchemeDesc& schemeDesc,
        bool localVol, Real illegalLocalVolOverwrite,
        CashDividendModel cashDividendModel,
        FdmBlackScholesSpatialDesc spatialDesc)
    : process_(std::move(process)), tGrid_(tGrid), xGrid_(xGrid),
      dampingSteps_(dampingSteps), schemeDesc_(schemeDesc), localVol_(localVol),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      cashDividendModel_(cashDividendModel),
      spatialDesc_(spatialDesc) {
        registerWith(process_);
    }

    FdBlackScholesVanillaEngine::FdBlackScholesVanillaEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        DividendSchedule dividends,
        Size tGrid, Size xGrid, Size dampingSteps,
        const FdmSchemeDesc& schemeDesc,
        bool localVol, Real illegalLocalVolOverwrite,
        CashDividendModel cashDividendModel,
        FdmBlackScholesSpatialDesc spatialDesc)
    : process_(std::move(process)), dividends_(std::move(dividends)),
      tGrid_(tGrid), xGrid_(xGrid), dampingSteps_(dampingSteps),
      schemeDesc_(schemeDesc), localVol_(localVol),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      cashDividendModel_(cashDividendModel),
      spatialDesc_(spatialDesc) {
        registerWith(process_);
    }

    FdBlackScholesVanillaEngine::FdBlackScholesVanillaEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        ext::shared_ptr<FdmQuantoHelper> quantoHelper,
        Size tGrid, Size xGrid, Size dampingSteps,
        const FdmSchemeDesc& schemeDesc,
        bool localVol, Real illegalLocalVolOverwrite,
        CashDividendModel cashDividendModel,
        FdmBlackScholesSpatialDesc spatialDesc)
    : process_(std::move(process)),
      tGrid_(tGrid), xGrid_(xGrid), dampingSteps_(dampingSteps),
      schemeDesc_(schemeDesc), localVol_(localVol),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      quantoHelper_(std::move(quantoHelper)),
      cashDividendModel_(cashDividendModel),
      spatialDesc_(spatialDesc) {
        registerWith(process_);
        registerWith(quantoHelper_);
    }

    FdBlackScholesVanillaEngine::FdBlackScholesVanillaEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        DividendSchedule dividends,
        ext::shared_ptr<FdmQuantoHelper> quantoHelper,
        Size tGrid, Size xGrid, Size dampingSteps,
        const FdmSchemeDesc& schemeDesc,
        bool localVol, Real illegalLocalVolOverwrite,
        CashDividendModel cashDividendModel,
        FdmBlackScholesSpatialDesc spatialDesc)
    : process_(std::move(process)), dividends_(std::move(dividends)),
      tGrid_(tGrid), xGrid_(xGrid), dampingSteps_(dampingSteps),
      schemeDesc_(schemeDesc), localVol_(localVol),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      quantoHelper_(std::move(quantoHelper)),
      cashDividendModel_(cashDividendModel),
      spatialDesc_(spatialDesc) {
        registerWith(process_);
        registerWith(quantoHelper_);
    }

    void FdBlackScholesVanillaEngine::calculate() const {

        const Date exerciseDate = arguments_.exercise->lastDate();
        const Time maturity = process_->time(exerciseDate);
        const Date settlementDate = process_->riskFreeRate()->referenceDate();

        Real spotAdjustment = 0.0;
        DividendSchedule dividendSchedule = DividendSchedule();

        ext::shared_ptr<EscrowedDividendAdjustment> escrowedDivAdj;

        switch (cashDividendModel_) {
          case Spot:
            dividendSchedule = dividends_;
            break;
          case Escrowed:
            if  (arguments_.exercise->type() != Exercise::European)
                for (const auto& cf: dividends_)
                    dividendSchedule.push_back(
                        ext::make_shared<FixedDividend>(0.0, cf->date()));

            QL_REQUIRE(quantoHelper_ == nullptr,
                "Escrowed dividend model is not supported for Quanto-Options");

            escrowedDivAdj = ext::make_shared<EscrowedDividendAdjustment>(
                dividends_,
                process_->riskFreeRate(),
                process_->dividendYield(),
                [&](Date d){ return process_->time(d); },
                maturity
            );

            spotAdjustment =
                escrowedDivAdj->dividendAdjustment(process_->time(settlementDate));

            QL_REQUIRE(process_->x0() + spotAdjustment > 0.0,
                    "spot minus dividends becomes negative");

            break;
          default:
              QL_FAIL("unknwon cash dividend model");
        }

        const ext::shared_ptr<StrikedTypePayoff> payoff =
            ext::dynamic_pointer_cast<StrikedTypePayoff>(arguments_.payoff);

        const ext::shared_ptr<Fdm1dMesher> equityMesher =
            ext::make_shared<FdmBlackScholesMesher>(
                    xGrid_, process_, maturity, payoff->strike(),
                    Null<Real>(), Null<Real>(), 0.0001, 1.5,
                    std::pair<Real, Real>(payoff->strike(), 0.1),
                    dividendSchedule, quantoHelper_,
                    spotAdjustment);

        const ext::shared_ptr<FdmMesher> mesher =
            ext::make_shared<FdmMesherComposite>(equityMesher);

        const ext::shared_ptr<FdmInnerValueCalculator> calculator =
            ext::make_shared<FdmLogInnerValue>(payoff, mesher, 0);

        ext::shared_ptr<FdmInnerValueCalculator> earlyExerciseCalculator;
        switch (cashDividendModel_) {
          case Spot:
              earlyExerciseCalculator = calculator;
            break;
          case Escrowed:
              earlyExerciseCalculator = ext::make_shared<FdmEscrowedLogInnerValueCalculator>(
                  escrowedDivAdj, payoff, mesher, 0);
            break;
          default:
              QL_FAIL("unknwon cash dividend model");
        }

        const ext::shared_ptr<FdmStepConditionComposite> conditions =
            FdmStepConditionComposite::vanillaComposite(
                dividendSchedule, arguments_.exercise, mesher,
                earlyExerciseCalculator,
                process_->riskFreeRate()->referenceDate(),
                process_->riskFreeRate()->dayCounter());

        const FdmBoundaryConditionSet boundaries;

        FdmSolverDesc solverDesc = {
            mesher, boundaries, conditions, calculator, maturity, tGrid_, dampingSteps_
        };

        const ext::shared_ptr<FdmBlackScholesSolver> solver(
            ext::make_shared<FdmBlackScholesSolver>(
                Handle<GeneralizedBlackScholesProcess>(process_),
                payoff->strike(), solverDesc, schemeDesc_,
                localVol_, illegalLocalVolOverwrite_,
                Handle<FdmQuantoHelper>(quantoHelper_),
                spatialDesc_));

        const Real spot = process_->x0() + spotAdjustment;

        results_.value = solver->valueAt(spot);
        results_.delta = solver->deltaAt(spot);
        results_.gamma = solver->gammaAt(spot);
        results_.theta = solver->thetaAt(spot);
    }

    MakeFdBlackScholesVanillaEngine::MakeFdBlackScholesVanillaEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process)
    : process_(std::move(process)),
      schemeDesc_(ext::make_shared<FdmSchemeDesc>(FdmSchemeDesc::Douglas())),
      illegalLocalVolOverwrite_(-Null<Real>()) {}

    MakeFdBlackScholesVanillaEngine&
    MakeFdBlackScholesVanillaEngine::withQuantoHelper(
        const ext::shared_ptr<FdmQuantoHelper>& quantoHelper) {
        quantoHelper_ = quantoHelper;
        return *this;
    }

    MakeFdBlackScholesVanillaEngine&
    MakeFdBlackScholesVanillaEngine::withTGrid(Size tGrid) {
        tGrid_ = tGrid;
        return *this;
    }

    MakeFdBlackScholesVanillaEngine&
    MakeFdBlackScholesVanillaEngine::withXGrid(Size xGrid) {
        xGrid_ = xGrid;
        return *this;
    }

    MakeFdBlackScholesVanillaEngine&
    MakeFdBlackScholesVanillaEngine::withDampingSteps(Size dampingSteps) {
        dampingSteps_ = dampingSteps;
        return *this;
    }

    MakeFdBlackScholesVanillaEngine&
    MakeFdBlackScholesVanillaEngine::withFdmSchemeDesc(
        const FdmSchemeDesc& schemeDesc) {
        schemeDesc_ = ext::make_shared<FdmSchemeDesc>(schemeDesc);
        return *this;
    }

    MakeFdBlackScholesVanillaEngine&
    MakeFdBlackScholesVanillaEngine::withLocalVol(bool localVol) {
        localVol_ = localVol;
        return *this;
    }

    MakeFdBlackScholesVanillaEngine&
    MakeFdBlackScholesVanillaEngine::withIllegalLocalVolOverwrite(
        Real illegalLocalVolOverwrite) {
        illegalLocalVolOverwrite_ = illegalLocalVolOverwrite;
        return *this;
    }

    MakeFdBlackScholesVanillaEngine&
    MakeFdBlackScholesVanillaEngine::withCashDividends(
            const std::vector<Date>& dividendDates,
            const std::vector<Real>& dividendAmounts) {
        dividends_ = DividendVector(dividendDates, dividendAmounts);
        return *this;
    }

    MakeFdBlackScholesVanillaEngine&
    MakeFdBlackScholesVanillaEngine::withCashDividendModel(
        FdBlackScholesVanillaEngine::CashDividendModel cashDividendModel) {
        cashDividendModel_ = cashDividendModel;
        return *this;
    }

    MakeFdBlackScholesVanillaEngine&
    MakeFdBlackScholesVanillaEngine::withSpatialDesc(
        const FdmBlackScholesSpatialDesc& spatialDesc) {
        spatialDesc_ = spatialDesc;
        return *this;
    }

    MakeFdBlackScholesVanillaEngine::operator
    ext::shared_ptr<PricingEngine>() const {
        return ext::make_shared<FdBlackScholesVanillaEngine>(
                process_, dividends_, quantoHelper_,
                tGrid_, xGrid_, dampingSteps_,
                *schemeDesc_, localVol_, illegalLocalVolOverwrite_,
                cashDividendModel_, spatialDesc_);
    }

}
