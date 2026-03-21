// r6
/*! \file fdmblackscholesmesher.cpp
    \brief 1-d mesher for the Black-Scholes process (in ln(S))
*/

#include <ql/processes/blackscholesprocess.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/termstructures/yield/quantotermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/methods/finitedifferences/utilities/fdmquantohelper.hpp>
#include <ql/methods/finitedifferences/meshers/uniform1dmesher.hpp>
#include <ql/methods/finitedifferences/meshers/concentrating1dmesher.hpp>
#include <ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp>

namespace QuantLib {

    namespace {
        // Shared domain-bound computation used by both constructors.
        // Returns (xMin, xMax).
        std::pair<Real, Real> computeDomainBounds(
                Size /*size*/,
                const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                Time maturity, Real strike,
                Real xMinConstraint, Real xMaxConstraint,
                Real eps, Real scaleFactor,
                const DividendSchedule& dividendSchedule,
                const ext::shared_ptr<FdmQuantoHelper>& fdmQuantoHelper,
                Real spotAdjustment) {

            const Real S = process->x0();
            QL_REQUIRE(S > 0.0, "negative or null underlying given");

            std::vector<std::pair<Time, Real> > intermediateSteps;
            for (const auto& i : dividendSchedule) {
                const Time t = process->time(i->date());
                if (t <= maturity && t >= 0.0)
                    intermediateSteps.emplace_back(
                        process->time(i->date()), i->amount());
            }

            const Size intermediateTimeSteps =
                std::max<Size>(2, Size(24.0*maturity));
            for (Size i=0; i < intermediateTimeSteps; ++i)
                intermediateSteps.emplace_back(
                    (i + 1) * (maturity / intermediateTimeSteps), 0.0);

            std::sort(intermediateSteps.begin(), intermediateSteps.end());

            const Handle<YieldTermStructure> rTS =
                process->riskFreeRate();

            const Handle<YieldTermStructure> qTS =
                (fdmQuantoHelper) != nullptr
                    ? Handle<YieldTermStructure>(
                          ext::make_shared<QuantoTermStructure>(
                              process->dividendYield(),
                              process->riskFreeRate(),
                              Handle<YieldTermStructure>(fdmQuantoHelper->fTS_),
                              process->blackVolatility(), strike,
                              Handle<BlackVolTermStructure>(
                                  fdmQuantoHelper->fxVolTS_),
                              fdmQuantoHelper->exchRateATMlevel_,
                              fdmQuantoHelper->equityFxCorrelation_))
                    : process->dividendYield();

            Time lastDivTime = 0.0;
            Real fwd = S + spotAdjustment;
            Real mi = fwd, ma = fwd;

            for (auto& intermediateStep : intermediateSteps) {
                const Time divTime  = intermediateStep.first;
                const Real divAmount = intermediateStep.second;

                fwd = fwd / rTS->discount(divTime)
                          * rTS->discount(lastDivTime)
                          * qTS->discount(divTime)
                          / qTS->discount(lastDivTime);

                mi = std::min(mi, fwd);
                ma = std::max(ma, fwd);

                fwd -= divAmount;

                mi = std::min(mi, fwd);
                ma = std::max(ma, fwd);

                lastDivTime = divTime;
            }

            const Real normInvEps =
                InverseCumulativeNormal()(1 - eps);
            const Real sigmaSqrtT =
                process->blackVolatility()->blackVol(maturity, strike)
                * std::sqrt(maturity);

            Real xMin =
                std::log(mi) - sigmaSqrtT * normInvEps * scaleFactor;
            Real xMax =
                std::log(ma) + sigmaSqrtT * normInvEps * scaleFactor;

            if (xMinConstraint != Null<Real>())
                xMin = xMinConstraint;
            if (xMaxConstraint != Null<Real>())
                xMax = xMaxConstraint;

            return {xMin, xMax};
        }
    }

    // ---------------------------------------------------------------
    //  Original single-point constructor (unchanged behaviour)
    // ---------------------------------------------------------------
    FdmBlackScholesMesher::FdmBlackScholesMesher(
        Size size,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Time maturity, Real strike,
        Real xMinConstraint, Real xMaxConstraint,
        Real eps, Real scaleFactor,
        const std::pair<Real, Real>& cPoint,
        const DividendSchedule& dividendSchedule,
        const ext::shared_ptr<FdmQuantoHelper>& fdmQuantoHelper,
        Real spotAdjustment)
    : Fdm1dMesher(size) {

        const auto bounds = computeDomainBounds(
            size, process, maturity, strike,
            xMinConstraint, xMaxConstraint,
            eps, scaleFactor,
            dividendSchedule, fdmQuantoHelper, spotAdjustment);

        const Real xMin = bounds.first;
        const Real xMax = bounds.second;

        ext::shared_ptr<Fdm1dMesher> helper;
        if (   cPoint.first != Null<Real>()
            && std::log(cPoint.first) >= xMin
            && std::log(cPoint.first) <= xMax) {

            helper = ext::shared_ptr<Fdm1dMesher>(
                new Concentrating1dMesher(
                    xMin, xMax, size,
                    std::pair<Real, Real>(
                        std::log(cPoint.first), cPoint.second)));
        }
        else {
            helper = ext::shared_ptr<Fdm1dMesher>(
                new Uniform1dMesher(xMin, xMax, size));
        }

        locations_ = helper->locations();
        for (Size i = 0; i < locations_.size(); ++i) {
            dplus_[i]  = helper->dplus(i);
            dminus_[i] = helper->dminus(i);
        }
    }

    // ---------------------------------------------------------------
    //  Multi-point constructor (spot-space → log-space conversion)
    // ---------------------------------------------------------------
    FdmBlackScholesMesher::FdmBlackScholesMesher(
        Size size,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Time maturity, Real strike,
        Real xMinConstraint, Real xMaxConstraint,
        Real eps, Real scaleFactor,
        const std::vector<std::tuple<Real, Real, bool>>& cPoints,
        const DividendSchedule& dividendSchedule,
        const ext::shared_ptr<FdmQuantoHelper>& fdmQuantoHelper,
        Real spotAdjustment)
    : Fdm1dMesher(size) {

        const auto bounds = computeDomainBounds(
            size, process, maturity, strike,
            xMinConstraint, xMaxConstraint,
            eps, scaleFactor,
            dividendSchedule, fdmQuantoHelper, spotAdjustment);

        const Real xMin = bounds.first;
        const Real xMax = bounds.second;

        // Convert spot-space critical points to log-space,
        // filtering to those inside the domain.
        std::vector<std::tuple<Real, Real, bool>> logCPoints;
        for (const auto& cp : cPoints) {
            const Real spotLevel = std::get<0>(cp);
            const Real density   = std::get<1>(cp);
            const bool required  = std::get<2>(cp);

            if (spotLevel > 0.0) {
                const Real logLevel = std::log(spotLevel);
                if (logLevel >= xMin && logLevel <= xMax) {
                    logCPoints.emplace_back(logLevel, density, required);
                }
            }
        }

        ext::shared_ptr<Fdm1dMesher> helper;
        if (!logCPoints.empty()) {
            helper = ext::make_shared<Concentrating1dMesher>(
                xMin, xMax, size, logCPoints, 1e-8);
        }
        else {
            helper = ext::make_shared<Uniform1dMesher>(xMin, xMax, size);
        }

        locations_ = helper->locations();
        for (Size i = 0; i < locations_.size(); ++i) {
            dplus_[i]  = helper->dplus(i);
            dminus_[i] = helper->dminus(i);
        }
    }

    ext::shared_ptr<GeneralizedBlackScholesProcess>
    FdmBlackScholesMesher::processHelper(
            const Handle<Quote>& s0,
            const Handle<YieldTermStructure>& rTS,
            const Handle<YieldTermStructure>& qTS,
            Volatility vol) {

        return ext::make_shared<GeneralizedBlackScholesProcess>(
            s0, qTS, rTS,
            Handle<BlackVolTermStructure>(
                ext::shared_ptr<BlackVolTermStructure>(
                    new BlackConstantVol(rTS->referenceDate(),
                                         Calendar(),
                                         vol,
                                         rTS->dayCounter()))));
    }
}
