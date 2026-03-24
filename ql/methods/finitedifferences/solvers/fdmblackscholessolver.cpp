// ══════════════════════════════════════════════════════════════════
// FdmBlackScholesSolver — 1-D Black-Scholes FDM solver
//
// Wraps FdmBlackScholesOp + Fdm1DimSolver with two additional
// responsibilities:
//
//   1. CN-equivalence gating [MT10, Theorem 3.2]:
//      The Milev-Tagliani scheme requires pure Crank-Nicolson
//      time stepping (theta = 0.5) with zero damping steps.
//      If either condition fails, the solver silently downgrades
//      to ExponentialFitting.  This is now observable via
//      solverGatingTriggered().
//
//   2. Lazy recalculation:
//      Inherits from LazyObject so that the operator and 1-D solver
//      are rebuilt only when market data changes.
//
// Design rationale for gating at the solver level:
//   The solver is the only component that knows *both* the spatial
//   descriptor AND the time-stepping scheme.  The operator only
//   sees the spatial descriptor; the time-stepping scheme is set
//   by Fdm1DimSolver.  Therefore gating must happen here, before
//   the operator is constructed, so that the operator receives the
//   effective (possibly downgraded) spatial descriptor.
// ══════════════════════════════════════════════════════════════════

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

    // ── CN-equivalence test [MT10, Theorem 3.2] ─────────────────
    // In 1D, Douglas and CraigSneyd both reduce to CN because
    // apply_mixed() returns zero (no cross-derivative term).
    // The only requirement is theta = 0.5 within numerical tolerance.
    bool FdmBlackScholesSolver::isCrankNicolsonEquivalent1D(
            const FdmSchemeDesc& schemeDesc, Real tol) {
        switch (schemeDesc.type) {
          case FdmSchemeDesc::CrankNicolsonType:
          case FdmSchemeDesc::DouglasType:
          case FdmSchemeDesc::CraigSneydType:
            // Douglas and CraigSneyd reduce to CN in 1D because
            // the mixed-derivative operator is zero.
            return std::fabs(schemeDesc.theta - 0.5) <= tol;
          default:
            return false;
        }
    }

    void FdmBlackScholesSolver::performCalculations() const {

        FdmBlackScholesSpatialDesc effectiveDesc = spatialDesc_;
        solverGatingTriggered_ = false;

        // ── Gating: MT scheme requires CN-equivalent time stepping ──
        // [MT10, Theorem 3.2] proves non-negative solutions only under
        // pure Crank-Nicolson (theta = 0.5).  Damping steps prepend
        // Implicit Euler sub-steps before the main scheme, breaking
        // CN-equivalence even when theta = 0.5.
        //
        // The paper also requires dt < 1/(r*M) for positive distinct
        // eigenvalues.  The exact bound in log-space may differ from
        // the S-space result; no runtime check is enforced here.
        if (effectiveDesc.scheme ==
                FdmBlackScholesSpatialDesc::Scheme
                    ::MilevTaglianiCNEffectiveDiffusion) {
            if (!isCrankNicolsonEquivalent1D(schemeDesc_)
                || solverDesc_.dampingSteps > 0) {
                if (effectiveDesc.mMatrixPolicy ==
                        FdmBlackScholesSpatialDesc::MMatrixPolicy::FailFast) {
                    QL_REQUIRE(false,
                        "Scheme-2 (Milev-Tagliani CN effective diffusion) "
                        "requires a Crank-Nicolson-equivalent time scheme "
                        "in 1D (CrankNicolsonType, DouglasType, or "
                        "CraigSneydType with theta=0.5) and "
                        "dampingSteps=0");
                }
                // Downgrade to ExponentialFitting — record for
                // observability via solverGatingTriggered().
                effectiveDesc.scheme =
                    FdmBlackScholesSpatialDesc::Scheme::ExponentialFitting;
                solverGatingTriggered_ = true;
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
        calculate();
        return solver_->thetaAt(std::log(s));
    }
}
