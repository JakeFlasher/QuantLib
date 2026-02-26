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
        const Real xLower_;
        const Real xUpper_;
        const Size direction_;
        std::vector<Size> knockOutIndices_;
    };

}

#endif
