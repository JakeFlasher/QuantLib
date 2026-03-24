// cf. [Ballabio20, Ch. 11] for the QuantLib step-condition framework.
// r6
#include <ql/methods/finitedifferences/stepconditions/fdmdiscretebarrierstepcondition.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/errors.hpp>
#include <algorithm>
#include <cmath>
#include <utility>

namespace QuantLib {

    // ── Backward-compatible constructor (rebate = 0) ─────────────
    // Delegates to the full constructor with zero rebate, matching
    // the standard knock-out convention where the option becomes
    // worthless upon barrier breach.
    FdmDiscreteBarrierStepCondition::FdmDiscreteBarrierStepCondition(
        ext::shared_ptr<FdmMesher> mesher,
        std::vector<Time> monitoringTimes,
        Real lowerBarrier,
        Real upperBarrier,
        Size direction)
    : FdmDiscreteBarrierStepCondition(
          std::move(mesher), std::move(monitoringTimes),
          lowerBarrier, upperBarrier, 0.0, direction) {}

    // ── Full constructor ─────────────────────────────────────────
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
      // Store log-barriers to compare against the log-space mesh.
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

        // Sort and deduplicate monitoring times so that applyTo()
        // can use std::lower_bound for O(log n) matching.
        std::sort(monitoringTimes_.begin(), monitoringTimes_.end());
        monitoringTimes_.erase(
            std::unique(monitoringTimes_.begin(),
                        monitoringTimes_.end()),
            monitoringTimes_.end());

        // Precompute which mesh nodes lie outside the [xLower, xUpper]
        // corridor.  These are the nodes that get overwritten with the
        // rebate value at each monitoring date (the knock-out indicator).
        for (const auto& iter : *mesher_->layout()) {
            const Real x = mesher_->location(iter, direction_);
            if (x < xLower_ || x > xUpper_)
                knockOutIndices_.push_back(iter.index());
        }
    }

    // ── applyTo — the step condition callback ────────────────────
    // Called by the FDM solver at each time step.  If the current
    // time t matches a monitoring date (within tolerance), overwrite
    // all knocked-out nodes with the rebate value.
    void FdmDiscreteBarrierStepCondition::applyTo(Array& a, Time t) const {
        if (monitoringTimes_.empty())
            return;

        // ── Tolerant time matching ───────────────────────────────
        // Uses std::lower_bound + relative tolerance instead of exact
        // binary_search.  The relative tolerance 1e-10 with
        // max(|candidate|, |query|) scaling avoids two failure modes:
        //   1. Exact matching fails when monitoring times have
        //      floating-point rounding from day-count conversions.
        //   2. Absolute tolerance can spuriously match ultra-short
        //      maturities (t near zero).
        // cf. BL-20260321-relative-tolerance-scaling for the scale
        // factor choice.
        {
            const Real tol = 1e-10;
            auto closeEnough = [&](Real a, Real b) {
                return std::fabs(a - b)
                    <= tol * std::max(std::fabs(a), std::fabs(b));
            };

            auto it = std::lower_bound(
                monitoringTimes_.begin(), monitoringTimes_.end(), t);

            // Check the lower_bound result and its predecessor.
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

        // Overwrite knocked-out nodes with rebate (0 for standard
        // knock-out, or a specified rebate amount).
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
