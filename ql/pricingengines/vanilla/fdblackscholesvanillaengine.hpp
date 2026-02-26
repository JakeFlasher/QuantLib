// r6
/*! \file fdblackscholesvanillaengine.hpp
    \brief Finite-differences Black Scholes vanilla option engine
*/

#ifndef quantlib_fd_black_scholes_vanilla_engine_hpp
#define quantlib_fd_black_scholes_vanilla_engine_hpp

#include <ql/pricingengine.hpp>
#include <ql/pricingengines/vanilla/cashdividendeuropeanengine.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>

namespace QuantLib {

    class FdmQuantoHelper;
    class GeneralizedBlackScholesProcess;

    class FdBlackScholesVanillaEngine : public VanillaOption::engine {
      public:
        enum CashDividendModel {
            Spot = CashDividendEuropeanEngine::Spot,
            Escrowed = CashDividendEuropeanEngine::Escrowed
        };

        explicit FdBlackScholesVanillaEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess>,
            Size tGrid = 100,
            Size xGrid = 100,
            Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            CashDividendModel cashDividendModel = Spot,
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        FdBlackScholesVanillaEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess>,
            DividendSchedule dividends,
            Size tGrid = 100,
            Size xGrid = 100,
            Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            CashDividendModel cashDividendModel = Spot,
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        FdBlackScholesVanillaEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess>,
            ext::shared_ptr<FdmQuantoHelper> quantoHelper,
            Size tGrid = 100,
            Size xGrid = 100,
            Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            CashDividendModel cashDividendModel = Spot,
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        FdBlackScholesVanillaEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess>,
            DividendSchedule dividends,
            ext::shared_ptr<FdmQuantoHelper> quantoHelper,
            Size tGrid = 100,
            Size xGrid = 100,
            Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            CashDividendModel cashDividendModel = Spot,
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        void calculate() const override;

      private:
        ext::shared_ptr<GeneralizedBlackScholesProcess> process_;
        DividendSchedule dividends_;
        Size tGrid_, xGrid_, dampingSteps_;
        FdmSchemeDesc schemeDesc_;
        bool localVol_;
        Real illegalLocalVolOverwrite_;
        ext::shared_ptr<FdmQuantoHelper> quantoHelper_;
        CashDividendModel cashDividendModel_;
        FdmBlackScholesSpatialDesc spatialDesc_;
    };


    class MakeFdBlackScholesVanillaEngine {
      public:
        explicit MakeFdBlackScholesVanillaEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process);

        MakeFdBlackScholesVanillaEngine& withQuantoHelper(
            const ext::shared_ptr<FdmQuantoHelper>& quantoHelper);

        MakeFdBlackScholesVanillaEngine& withTGrid(Size tGrid);
        MakeFdBlackScholesVanillaEngine& withXGrid(Size xGrid);
        MakeFdBlackScholesVanillaEngine& withDampingSteps(
            Size dampingSteps);

        MakeFdBlackScholesVanillaEngine& withFdmSchemeDesc(
            const FdmSchemeDesc& schemeDesc);

        MakeFdBlackScholesVanillaEngine& withLocalVol(bool localVol);
        MakeFdBlackScholesVanillaEngine& withIllegalLocalVolOverwrite(
            Real illegalLocalVolOverwrite);

        MakeFdBlackScholesVanillaEngine& withCashDividends(
            const std::vector<Date>& dividendDates,
            const std::vector<Real>& dividendAmounts);

        MakeFdBlackScholesVanillaEngine& withCashDividendModel(
            FdBlackScholesVanillaEngine::CashDividendModel cashDividendModel);

        MakeFdBlackScholesVanillaEngine& withSpatialDesc(
            const FdmBlackScholesSpatialDesc& spatialDesc);

        operator ext::shared_ptr<PricingEngine>() const;
      private:
        ext::shared_ptr<GeneralizedBlackScholesProcess> process_;
        DividendSchedule dividends_;
        Size tGrid_ = 100, xGrid_ = 100, dampingSteps_ = 0;
        ext::shared_ptr<FdmSchemeDesc> schemeDesc_;
        bool localVol_ = false;
        Real illegalLocalVolOverwrite_;
        ext::shared_ptr<FdmQuantoHelper> quantoHelper_;
        FdBlackScholesVanillaEngine::CashDividendModel cashDividendModel_
            = FdBlackScholesVanillaEngine::Spot;
        FdmBlackScholesSpatialDesc spatialDesc_;
    };

}

#endif
