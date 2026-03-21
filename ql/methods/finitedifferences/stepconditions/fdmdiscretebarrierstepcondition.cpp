// r6
#include <ql/methods/finitedifferences/stepconditions/fdmdiscretebarrierstepcondition.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/errors.hpp>
#include <algorithm>
#include <cmath>
#include <utility>

namespace QuantLib {

    FdmDiscreteBarrierStepCondition::FdmDiscreteBarrierStepCondition(
        ext::shared_ptr<FdmMesher> mesher,
        std::vector<Time> monitoringTimes,
        Real lowerBarrier,
        Real upperBarrier,
        Size direction)
    : FdmDiscreteBarrierStepCondition(
          std::move(mesher), std::move(monitoringTimes),
          lowerBarrier, upperBarrier, 0.0, direction) {}

    FdmDiscreteBarrierStepCondition::FdmDiscreteBarrierStepCondition(
        ext::shared_ptr<FdmMesher> mesher,
        std::vector<Time> monitoringTimes,
        Real lowerBarrier,
        Real upperBarrier,
        Real rebate,
        Size direction)
    : mesher_(std::move(mesher)),
      monitoringTimes_(std::move(monitoringTimes)),
      lowerBarrier_(lowerBarrier),
      upperBarrier_(upperBarrier),
      rebate_(rebate),
      xLower_(std::log(lowerBarrier)),
      xUpper_(std::log(upperBarrier)),
      direction_(direction) {

        QL_REQUIRE(mesher_ != nullptr, "null mesher");
        QL_REQUIRE(lowerBarrier_ > 0.0,
                   "lower barrier (" << lowerBarrier_
                   << ") must be positive");
        QL_REQUIRE(upperBarrier_ > lowerBarrier_,
                   "upper barrier (" << upperBarrier_
                   << ") must exceed lower barrier ("
                   << lowerBarrier_ << ")");
        QL_REQUIRE(direction_ < mesher_->layout()->dim().size(),
                   "direction (" << direction_ << ") out of range");

        std::sort(monitoringTimes_.begin(), monitoringTimes_.end());
        monitoringTimes_.erase(
            std::unique(monitoringTimes_.begin(),
                        monitoringTimes_.end()),
            monitoringTimes_.end());

        for (const auto& iter : *mesher_->layout()) {
            const Real x = mesher_->location(iter, direction_);
            if (x < xLower_ || x > xUpper_)
                knockOutIndices_.push_back(iter.index());
        }
    }

    void FdmDiscreteBarrierStepCondition::applyTo(Array& a, Time t) const {
        if (monitoringTimes_.empty())
            return;

        // Tolerant time matching: use lower_bound and check both the
        // found element and its predecessor for closeness within a
        // relative tolerance.  This is more robust than exact
        // binary_search when monitoring times suffer floating-point
        // representation or rounding differences.
        {
            const Real tol = 1e-10;
            auto closeEnough = [&](Real a, Real b) {
                return std::fabs(a - b)
                    <= tol * std::max(std::fabs(a), std::fabs(b));
            };

            auto it = std::lower_bound(
                monitoringTimes_.begin(), monitoringTimes_.end(), t);

            bool matched = (it != monitoringTimes_.end()
                            && closeEnough(*it, t));
            if (!matched && it != monitoringTimes_.begin())
                matched = closeEnough(*std::prev(it), t);
            if (!matched)
                return;
        }

        QL_REQUIRE(a.size() == mesher_->layout()->size(),
                   "array size (" << a.size()
                   << ") does not match layout size ("
                   << mesher_->layout()->size() << ")");

        for (Size idx : knockOutIndices_)
            a[idx] = rebate_;
    }

    const std::vector<Time>&
    FdmDiscreteBarrierStepCondition::monitoringTimes() const {
        return monitoringTimes_;
    }

    Real FdmDiscreteBarrierStepCondition::lowerBarrier() const {
        return lowerBarrier_;
    }

    Real FdmDiscreteBarrierStepCondition::upperBarrier() const {
        return upperBarrier_;
    }

    Real FdmDiscreteBarrierStepCondition::rebate() const {
        return rebate_;
    }
}
