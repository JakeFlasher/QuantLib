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
    };
}

#endif
