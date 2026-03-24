// ══════════════════════════════════════════════════════════════════
// FdmBlackScholesSpatialDesc — spatial discretization descriptor
// for 1-D Black-Scholes finite-difference operators.
//
// This POD-like struct bundles all knobs that select *how* the
// spatial (x = ln S) discretization modifies the standard
// central-difference operator.  Three schemes are supported:
//
//   StandardCentral              — vanilla 2nd-order central FD
//   ExponentialFitting           — Duffy (2004) coth-fitting
//                                  [Duffy04, §4, Eq. 12-13]
//   MilevTaglianiCNEffectiveDiffusion
//                                — Milev & Tagliani (2010)
//                                  added-diffusion Scheme-2
//                                  [MT10, §3.2]
//
// The struct also controls the M-matrix diagnostic policy and the
// mesh-spacing rule used when computing Peclet numbers.
//
// Usage: construct via the named factories standard(),
//        exponentialFitting(), or milevTaglianiCN(), then
//        optionally override individual fields.
// ══════════════════════════════════════════════════════════════════

/*! \file fdmblackscholesspatialdesc.hpp
    \brief spatial discretization descriptor for 1D Black-Scholes FDM operators
*/

#ifndef quantlib_fdm_black_scholes_spatial_desc_hpp
#define quantlib_fdm_black_scholes_spatial_desc_hpp

#include <ql/types.hpp>

namespace QuantLib {

    struct FdmBlackScholesSpatialDesc {

        // ── Spatial discretization scheme ────────────────────────
        // Selects the diffusion-coefficient modification applied at
        // each interior node before assembling the tridiagonal
        // operator.  StandardCentral leaves a(x) = sigma^2/2;
        // the two non-standard schemes augment it to improve
        // positivity or suppress spurious oscillations.
        enum class Scheme {
            StandardCentral,
            ExponentialFitting,
            MilevTaglianiCNEffectiveDiffusion
        };

        // ── Mesh-spacing rule for Peclet / fitting computations ──
        // On a non-uniform log-space mesh the left half-spacing
        // h_minus and right half-spacing h_plus may differ.  HPolicy
        // determines the single representative spacing h used in the
        // Peclet number Pe = mu*h / sigma^2 and in the MT added
        // diffusion r^2 h^2 / (8 sigma^2).
        //   Average   — h = (h_minus + h_plus) / 2   (default)
        //   Min       — h = min(h_minus, h_plus)
        //   Harmonic  — h = 2 h_minus h_plus / (h_minus + h_plus)
        enum class HPolicy {
            Average,
            Min,
            Harmonic
        };

        // ── M-matrix diagnostic / enforcement policy ─────────────
        // After assembling the tridiagonal operator the code can
        // optionally check whether the off-diagonal entries are
        // non-negative (the M-matrix property, guaranteeing
        // positivity preservation — cf. [MT10, Thm 3.1]).
        //   None                         — skip the check entirely
        //   FailFast                     — QL_REQUIRE failure on
        //                                  first violation
        //   FallbackToExponentialFitting — silently recompute with
        //                                  ExponentialFitting if the
        //                                  current scheme violates
        //                                  the M-matrix condition
        enum class MMatrixPolicy {
            None,
            FailFast,
            FallbackToExponentialFitting
        };

        // ── Scheme selection (default: StandardCentral) ──────────
        Scheme scheme = Scheme::StandardCentral;

        // ── Mesh-spacing rule (default: arithmetic mean) ─────────
        HPolicy hPolicy = HPolicy::Average;

        // ── Peclet-number thresholds for xCothx evaluation ───────
        // peSmall: below this |Pe|, use Taylor expansion
        //          1 + x^2/3 to avoid 0/0 in x/tanh(x)
        // peLarge: above this |Pe|, use asymptotic |x| since
        //          coth(x) ~ sign(x) for large arguments
        // cf. [Duffy04, §4] for the exponential fitting factor.
        Real peSmall  = 1e-6;
        Real peLarge  = 50.0;

        // ── Numerical safety floors ──────────────────────────────
        // minVariance: lower clamp on per-node variance to prevent
        //              division by zero in a_add = r^2 h^2/(8 v)
        // maxAddedDiffusionRatio: caps a_add at this multiple of
        //              the base diffusion sigma^2/2, preventing
        //              extreme added diffusion under low-vol /
        //              high-rate scenarios
        Real minVariance           = 1e-12;
        Real maxAddedDiffusionRatio = 1e6;

        // ── M-matrix diagnostic settings ─────────────────────────
        // mMatrixPolicy: what to do when off-diagonals are negative
        //                (default: fall back to exponential fitting)
        // mMatrixEps: tolerance for near-zero negatives; entries
        //             in [-eps, 0) are treated as non-negative
        // checkBoundaries: when true, include boundary rows in the
        //                  M-matrix check (normally only interior)
        MMatrixPolicy mMatrixPolicy = MMatrixPolicy::FallbackToExponentialFitting;
        Real mMatrixEps     = 0.0;
        bool checkBoundaries = false;

        // ── Named factory methods ────────────────────────────────
        // Return pre-configured descriptors for common use cases.
        // Callers can further customize individual fields after
        // construction.

        //! Standard central differences — baseline, no fitting.
        static FdmBlackScholesSpatialDesc standard() {
            return {};
        }

        //! Exponential fitting per [Duffy04, §4].
        static FdmBlackScholesSpatialDesc exponentialFitting() {
            FdmBlackScholesSpatialDesc d;
            d.scheme = Scheme::ExponentialFitting;
            return d;
        }

        //! Milev-Tagliani CN effective diffusion per [MT10, §3.2].
        static FdmBlackScholesSpatialDesc milevTaglianiCN() {
            FdmBlackScholesSpatialDesc d;
            d.scheme = Scheme::MilevTaglianiCNEffectiveDiffusion;
            return d;
        }
    };

}

#endif
