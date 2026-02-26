// r6
/*! \file fdblackscholesshoutengine.hpp
    \brief Finite-Differences Black Scholes shout option engine
*/

#ifndef quantlib_fd_black_scholes_shout_engine_hpp
#define quantlib_fd_black_scholes_shout_engine_hpp

#include <ql/pricingengine.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>

namespace QuantLib {

    class GeneralizedBlackScholesProcess;

    class FdBlackScholesShoutEngine : public VanillaOption::engine {
      public:
        explicit FdBlackScholesShoutEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess>,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        FdBlackScholesShoutEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess>,
            DividendSchedule dividends,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        void calculate() const override;

      private:
        const ext::shared_ptr<GeneralizedBlackScholesProcess> process_;
        DividendSchedule dividends_;
        const Size tGrid_, xGrid_, dampingSteps_;
        const FdmSchemeDesc schemeDesc_;
        const FdmBlackScholesSpatialDesc spatialDesc_;
    };

}

#endif
