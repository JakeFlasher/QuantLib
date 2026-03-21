// r6
/*! \file makefdblackscholesbarrierengine.hpp
    \brief Fluent builder for FdBlackScholesBarrierEngine
*/

#ifndef quantlib_make_fd_black_scholes_barrier_engine_hpp
#define quantlib_make_fd_black_scholes_barrier_engine_hpp

#include <ql/pricingengine.hpp>
#include <ql/instruments/barrieroption.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>
#include <ql/time/date.hpp>
#include <vector>

namespace QuantLib {

    class GeneralizedBlackScholesProcess;

    //! Fluent builder for FdBlackScholesBarrierEngine.
    //! Note: the underlying engine handles single-barrier instruments only.
    //! For double-barrier knock-out options, build the FD solver manually.
    class MakeFdBlackScholesBarrierEngine {
      public:
        explicit MakeFdBlackScholesBarrierEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process);

        MakeFdBlackScholesBarrierEngine& withTGrid(Size tGrid);
        MakeFdBlackScholesBarrierEngine& withXGrid(Size xGrid);
        MakeFdBlackScholesBarrierEngine& withDampingSteps(Size dampingSteps);

        MakeFdBlackScholesBarrierEngine& withFdmSchemeDesc(
            const FdmSchemeDesc& schemeDesc);

        MakeFdBlackScholesBarrierEngine& withLocalVol(bool localVol);
        MakeFdBlackScholesBarrierEngine& withIllegalLocalVolOverwrite(
            Real illegalLocalVolOverwrite);

        MakeFdBlackScholesBarrierEngine& withCashDividends(
            const std::vector<Date>& dividendDates,
            const std::vector<Real>& dividendAmounts);

        MakeFdBlackScholesBarrierEngine& withSpatialDesc(
            const FdmBlackScholesSpatialDesc& spatialDesc);

        MakeFdBlackScholesBarrierEngine& withDiscreteMonitoring(
            const std::vector<Date>& monitoringDates);

        operator ext::shared_ptr<PricingEngine>() const;

      private:
        ext::shared_ptr<GeneralizedBlackScholesProcess> process_;
        DividendSchedule dividends_;
        std::vector<Date> monitoringDates_;
        Size tGrid_ = 100, xGrid_ = 100, dampingSteps_ = 0;
        ext::shared_ptr<FdmSchemeDesc> schemeDesc_;
        bool localVol_ = false;
        Real illegalLocalVolOverwrite_;
        FdmBlackScholesSpatialDesc spatialDesc_;
    };

}

#endif
