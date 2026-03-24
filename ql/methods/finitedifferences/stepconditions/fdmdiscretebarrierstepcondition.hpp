// ══════════════════════════════════════════════════════════════════
// FdmDiscreteBarrierStepCondition
//
// Step condition for discrete barrier monitoring on a 1-D
// log-space grid.  At each monitoring date the condition
// overwrites grid values outside the [L, U] corridor with the
// knock-out rebate, implementing the indicator-function approach:
//
//   V(x, t_k) = rebate   if  x < ln(L)  or  x > ln(U)
//   V(x, t_k) = V(x, t_k)  otherwise
//
// where t_k are the discrete monitoring times.
//
// Key design choices:
//   - Monitoring times are sorted and deduplicated at construction
//     to guarantee O(log n) lookup via std::lower_bound.
//   - Time matching uses a relative tolerance of 1e-10 with
//     max(|candidate|, |query|) scaling, avoiding spurious matches
//     for short maturities — cf. BL-20260321-relative-tolerance-scaling.
//   - Knock-out indices are precomputed from the mesher at
//     construction time (O(n) once) so that applyTo() is O(k)
//     where k = number of knocked-out nodes.
//   - Two constructors: backward-compatible (rebate = 0) and full
//     (explicit rebate + direction).
// ══════════════════════════════════════════════════════════════════

// r6
/*! \file fdmdiscretebarrierstepcondition.hpp
    \brief Discrete barrier monitoring step condition for 1D log-space grids
*/

#ifndef quantlib_fdm_discrete_barrier_step_condition_hpp
#define quantlib_fdm_discrete_barrier_step_condition_hpp

#include <ql/methods/finitedifferences/stepcondition.hpp>
#include <ql/shared_ptr.hpp>
#include <vector>

namespace QuantLib {

    class FdmMesher;

    class FdmDiscreteBarrierStepCondition : public StepCondition<Array> {
      public:
        /*! Backward-compatible constructor (rebate = 0). */
        FdmDiscreteBarrierStepCondition(
            ext::shared_ptr<FdmMesher> mesher,
            std::vector<Time> monitoringTimes,
            Real lowerBarrier,
            Real upperBarrier,
            Size direction = 0);

        /*! Full constructor with rebate.
            \note direction has no default to avoid ambiguity with the
                  backward-compatible overload. */
        FdmDiscreteBarrierStepCondition(
            ext::shared_ptr<FdmMesher> mesher,
            std::vector<Time> monitoringTimes,
            Real lowerBarrier,
            Real upperBarrier,
            Real rebate,
            Size direction);

        void applyTo(Array& a, Time t) const override;

        const std::vector<Time>& monitoringTimes() const;
        Real lowerBarrier() const;
        Real upperBarrier() const;
        Real rebate() const;

      private:
        const ext::shared_ptr<FdmMesher> mesher_;
        std::vector<Time> monitoringTimes_;
        const Real lowerBarrier_;
        const Real upperBarrier_;
        const Real rebate_;
        const Real xLower_;       //!< ln(lowerBarrier) — precomputed
        const Real xUpper_;       //!< ln(upperBarrier) — precomputed
        const Size direction_;
        std::vector<Size> knockOutIndices_;  //!< mesh nodes outside corridor
    };

}

#endif
