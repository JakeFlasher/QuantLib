// r6
/*! \file fdblackscholesbarrierengine.hpp
    \brief Finite-differences Black/Scholes barrier-option engine
*/

#ifndef quantlib_fd_black_scholes_barrier_engine_hpp
#define quantlib_fd_black_scholes_barrier_engine_hpp

#include <ql/processes/blackscholesprocess.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>
#include <ql/instruments/barrieroption.hpp>
#include <vector>

namespace QuantLib {

    class FdBlackScholesBarrierEngine : public BarrierOption::engine {
      public:
        // --- Continuous monitoring constructors (unchanged API) ---

        explicit FdBlackScholesBarrierEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        explicit FdBlackScholesBarrierEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            DividendSchedule dividends,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        // --- Discrete monitoring constructors ---

        FdBlackScholesBarrierEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            std::vector<Date> monitoringDates,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        FdBlackScholesBarrierEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            DividendSchedule dividends,
            std::vector<Date> monitoringDates,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        void calculate() const override;

      private:
        void calculateContinuous() const;
        void calculateDiscrete() const;

        ext::shared_ptr<GeneralizedBlackScholesProcess> process_;
        DividendSchedule dividends_;
        std::vector<Date> monitoringDates_;
        Size tGrid_, xGrid_, dampingSteps_;
        FdmSchemeDesc schemeDesc_;
        bool localVol_;
        Real illegalLocalVolOverwrite_;
        FdmBlackScholesSpatialDesc spatialDesc_;
    };

}

#endif
