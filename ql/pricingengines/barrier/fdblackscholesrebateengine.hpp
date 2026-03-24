// spatialDesc is threaded through to FdmBlackScholesSolver so that
// the rebate pricing uses the same spatial discretization scheme
// (StandardCentral / ExponentialFitting / MilevTaglianiCN) as the
// main barrier engine that delegates to this helper.
// cf. [Ballabio20, Ch. 11] for the rebate engine framework.

// r6
/*! \file fdblackscholesrebateengine.hpp
    \brief Finite-differences Black/Scholes barrier option rebate helper engine
*/

#ifndef quantlib_fd_black_scholes_rebate_engine_hpp
#define quantlib_fd_black_scholes_rebate_engine_hpp

#include <ql/processes/blackscholesprocess.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>
#include <ql/instruments/barrieroption.hpp>

namespace QuantLib {

    class FdBlackScholesRebateEngine : public BarrierOption::engine {
      public:
        explicit FdBlackScholesRebateEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        explicit FdBlackScholesRebateEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            DividendSchedule dividends,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
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
        FdmBlackScholesSpatialDesc spatialDesc_;
    };

}

#endif
