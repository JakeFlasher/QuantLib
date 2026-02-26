// r6
#include <ql/methods/finitedifferences/operators/fdmblackscholesop.hpp>
#include <ql/methods/finitedifferences/solvers/fdm1dimsolver.hpp>
#include <ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <utility>
#include <cmath>

namespace QuantLib {

    FdmBlackScholesSolver::FdmBlackScholesSolver(
        Handle<GeneralizedBlackScholesProcess> process,
        Real strike,
        FdmSolverDesc solverDesc,
        const FdmSchemeDesc& schemeDesc,
        bool localVol,
        Real illegalLocalVolOverwrite,
        Handle<FdmQuantoHelper> quantoHelper,
        FdmBlackScholesSpatialDesc spatialDesc)
    : process_(std::move(process)),
      strike_(strike),
      solverDesc_(std::move(solverDesc)),
      schemeDesc_(schemeDesc),
      localVol_(localVol),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      quantoHelper_(std::move(quantoHelper)),
      spatialDesc_(spatialDesc) {

        registerWith(process_);
        registerWith(quantoHelper_);
    }

    bool FdmBlackScholesSolver::isCrankNicolsonEquivalent1D(
            const FdmSchemeDesc& schemeDesc, Real tol) {
        switch (schemeDesc.type) {
          case FdmSchemeDesc::CrankNicolsonType:
            return std::fabs(schemeDesc.theta - 0.5) <= tol;
          case FdmSchemeDesc::DouglasType:
            // Douglas reduces to CN in one dimension
            return std::fabs(schemeDesc.theta - 0.5) <= tol;
          default:
            return false;
        }
    }

    void FdmBlackScholesSolver::performCalculations() const {

        FdmBlackScholesSpatialDesc effectiveDesc = spatialDesc_;

        // Gating: Scheme-2 requires CN-equivalent time stepping in 1-D
        if (effectiveDesc.scheme ==
                FdmBlackScholesSpatialDesc::Scheme
                    ::MilevTaglianiCNEffectiveDiffusion) {
            if (!isCrankNicolsonEquivalent1D(schemeDesc_)) {
                if (effectiveDesc.mMatrixPolicy ==
                        FdmBlackScholesSpatialDesc::MMatrixPolicy::FailFast) {
                    QL_REQUIRE(false,
                        "Scheme-2 (Milev-Tagliani CN effective diffusion) "
                        "requires a Crank-Nicolson-equivalent time scheme "
                        "in 1D (CrankNicolsonType or DouglasType with "
                        "theta=0.5)");
                }
                effectiveDesc.scheme =
                    FdmBlackScholesSpatialDesc::Scheme::ExponentialFitting;
            }
        }

        const ext::shared_ptr<FdmBlackScholesOp> op(
            ext::make_shared<FdmBlackScholesOp>(
                solverDesc_.mesher, process_.currentLink(), strike_,
                localVol_, illegalLocalVolOverwrite_, 0,
                (quantoHelper_.empty())
                    ? ext::shared_ptr<FdmQuantoHelper>()
                    : quantoHelper_.currentLink(),
                effectiveDesc));

        solver_ = ext::make_shared<Fdm1DimSolver>(
            solverDesc_, schemeDesc_, op);
    }

    Real FdmBlackScholesSolver::valueAt(Real s) const {
        calculate();
        return solver_->interpolateAt(std::log(s));
    }

    Real FdmBlackScholesSolver::deltaAt(Real s) const {
        calculate();
        return solver_->derivativeX(std::log(s))/s;
    }

    Real FdmBlackScholesSolver::gammaAt(Real s) const {
        calculate();
        return (solver_->derivativeXX(std::log(s))
                -solver_->derivativeX(std::log(s)))/(s*s);
    }

    Real FdmBlackScholesSolver::thetaAt(Real s) const {
        return solver_->thetaAt(std::log(s));
    }
}
