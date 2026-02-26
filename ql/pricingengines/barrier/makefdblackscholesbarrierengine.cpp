// r6
#include <ql/pricingengines/barrier/makefdblackscholesbarrierengine.hpp>
#include <ql/pricingengines/barrier/fdblackscholesbarrierengine.hpp>
#include <ql/instruments/dividendschedule.hpp>
#include <ql/utilities/null.hpp>
#include <utility>

namespace QuantLib {

    MakeFdBlackScholesBarrierEngine::MakeFdBlackScholesBarrierEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process)
    : process_(std::move(process)),
      schemeDesc_(ext::make_shared<FdmSchemeDesc>(
          FdmSchemeDesc::Douglas())),
      illegalLocalVolOverwrite_(-Null<Real>()) {}

    MakeFdBlackScholesBarrierEngine&
    MakeFdBlackScholesBarrierEngine::withTGrid(Size tGrid) {
        tGrid_ = tGrid;
        return *this;
    }

    MakeFdBlackScholesBarrierEngine&
    MakeFdBlackScholesBarrierEngine::withXGrid(Size xGrid) {
        xGrid_ = xGrid;
        return *this;
    }

    MakeFdBlackScholesBarrierEngine&
    MakeFdBlackScholesBarrierEngine::withDampingSteps(Size dampingSteps) {
        dampingSteps_ = dampingSteps;
        return *this;
    }

    MakeFdBlackScholesBarrierEngine&
    MakeFdBlackScholesBarrierEngine::withFdmSchemeDesc(
            const FdmSchemeDesc& schemeDesc) {
        schemeDesc_ = ext::make_shared<FdmSchemeDesc>(schemeDesc);
        return *this;
    }

    MakeFdBlackScholesBarrierEngine&
    MakeFdBlackScholesBarrierEngine::withLocalVol(bool localVol) {
        localVol_ = localVol;
        return *this;
    }

    MakeFdBlackScholesBarrierEngine&
    MakeFdBlackScholesBarrierEngine::withIllegalLocalVolOverwrite(
            Real illegalLocalVolOverwrite) {
        illegalLocalVolOverwrite_ = illegalLocalVolOverwrite;
        return *this;
    }

    MakeFdBlackScholesBarrierEngine&
    MakeFdBlackScholesBarrierEngine::withCashDividends(
            const std::vector<Date>& dividendDates,
            const std::vector<Real>& dividendAmounts) {
        dividends_ = DividendVector(dividendDates, dividendAmounts);
        return *this;
    }

    MakeFdBlackScholesBarrierEngine&
    MakeFdBlackScholesBarrierEngine::withSpatialDesc(
            const FdmBlackScholesSpatialDesc& spatialDesc) {
        spatialDesc_ = spatialDesc;
        return *this;
    }

    MakeFdBlackScholesBarrierEngine&
    MakeFdBlackScholesBarrierEngine::withDiscreteMonitoring(
            const std::vector<Date>& monitoringDates) {
        monitoringDates_ = monitoringDates;
        return *this;
    }

    MakeFdBlackScholesBarrierEngine::operator
    ext::shared_ptr<PricingEngine>() const {
        return ext::make_shared<FdBlackScholesBarrierEngine>(
            process_, dividends_, monitoringDates_,
            tGrid_, xGrid_, dampingSteps_,
            *schemeDesc_, localVol_, illegalLocalVolOverwrite_,
            spatialDesc_);
    }

}
