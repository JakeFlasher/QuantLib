// ══════════════════════════════════════════════════════════════════
// FdmBlackScholesSolver — 1-D BS FDM solver with CN-equivalence
// gating for the Milev-Tagliani scheme [MT10, Thm 3.2].
//
// spatialDesc is threaded through from the engine; the solver
// applies gating logic (checking time-stepping scheme and damping
// steps) before constructing the operator.  If gating triggers,
// the scheme is downgraded to ExponentialFitting — now observable
// via solverGatingTriggered().
// ══════════════════════════════════════════════════════════════════

// r6
/*! \file fdmblackscholessolver.hpp
*/

#ifndef quantlib_fdm_black_scholes_solver_hpp
#define quantlib_fdm_black_scholes_solver_hpp

#include <ql/handle.hpp>
#include <ql/patterns/lazyobject.hpp>
#include <ql/methods/finitedifferences/solvers/fdmsolverdesc.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/utilities/fdmquantohelper.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>

namespace QuantLib {

    class Fdm1DimSolver;
    class FdmSnapshotCondition;
    class GeneralizedBlackScholesProcess;

    class FdmBlackScholesSolver : public LazyObject {
      public:
        FdmBlackScholesSolver(
            Handle<GeneralizedBlackScholesProcess> process,
            Real strike,
            FdmSolverDesc solverDesc,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            Handle<FdmQuantoHelper> quantoHelper = Handle<FdmQuantoHelper>(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        Real valueAt(Real s) const;
        Real deltaAt(Real s) const;
        Real gammaAt(Real s) const;
        Real thetaAt(Real s) const;

        //! True if solver gating downgraded the requested scheme
        //  (e.g. MT -> ExponentialFitting due to non-CN time stepping
        //  or nonzero damping steps).
        bool solverGatingTriggered() const {
            calculate();
            return solverGatingTriggered_;
        }

      protected:
        void performCalculations() const override;

      private:
        static bool isCrankNicolsonEquivalent1D(
            const FdmSchemeDesc& schemeDesc, Real tol = 1e-12);

        Handle<GeneralizedBlackScholesProcess> process_;
        const Real strike_;
        const FdmSolverDesc solverDesc_;
        const FdmSchemeDesc schemeDesc_;
        const bool localVol_;
        const Real illegalLocalVolOverwrite_;
        const Handle<FdmQuantoHelper> quantoHelper_;
        const FdmBlackScholesSpatialDesc spatialDesc_;

        mutable ext::shared_ptr<Fdm1DimSolver> solver_;
        mutable bool solverGatingTriggered_ = false;
    };
}

#endif
