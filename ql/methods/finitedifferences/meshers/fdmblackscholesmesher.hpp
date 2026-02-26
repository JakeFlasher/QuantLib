// r6
/*! \file fdmblackscholesmesher.hpp
    \brief 1-d mesher for the Black-Scholes process (in ln(S))
*/

#ifndef quantlib_fdm_black_scholes_mesher_hpp
#define quantlib_fdm_black_scholes_mesher_hpp

#include <ql/instruments/dividendschedule.hpp>
#include <ql/methods/finitedifferences/meshers/fdm1dmesher.hpp>

#include <ql/handle.hpp>
#include <ql/quote.hpp>
#include <tuple>

namespace QuantLib {

    class FdmQuantoHelper;
    class YieldTermStructure;
    class GeneralizedBlackScholesProcess;

    class FdmBlackScholesMesher : public Fdm1dMesher {
      public:
        FdmBlackScholesMesher(
            Size size,
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
            Time maturity, Real strike,
            Real xMinConstraint = Null<Real>(),
            Real xMaxConstraint = Null<Real>(),
            Real eps = 0.0001,
            Real scaleFactor = 1.5,
            const std::pair<Real, Real>& cPoint = { Null<Real>(), Null<Real>() },
            const DividendSchedule& dividendSchedule = {},
            const ext::shared_ptr<FdmQuantoHelper>& fdmQuantoHelper = {},
            Real spotAdjustment = 0.0);

        /*! Multi-point constructor.  Each tuple is
            (spot-space level, density, require-exact-node).
            Spot-space levels are converted to log-space internally. */
        FdmBlackScholesMesher(
            Size size,
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
            Time maturity, Real strike,
            Real xMinConstraint,
            Real xMaxConstraint,
            Real eps,
            Real scaleFactor,
            const std::vector<std::tuple<Real, Real, bool>>& cPoints,
            const DividendSchedule& dividendSchedule = {},
            const ext::shared_ptr<FdmQuantoHelper>& fdmQuantoHelper = {},
            Real spotAdjustment = 0.0);

        static ext::shared_ptr<GeneralizedBlackScholesProcess> processHelper(
             const Handle<Quote>& s0,
             const Handle<YieldTermStructure>& rTS,
             const Handle<YieldTermStructure>& qTS,
             Volatility vol);
    };
}

#endif
