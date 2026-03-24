// ══════════════════════════════════════════════════════════════════
// xCothx — numerically stable evaluation of f(x) = x * coth(x)
//
// This is the exponential fitting factor used in the
// Duffy (2004) diffusion-coefficient modification for
// central-difference operators on the Black-Scholes PDE.
//
// In the fitted scheme the diffusion coefficient becomes
//   a_fitted = (sigma^2 / 2) * rho
// where rho = xCothx(Pe) and Pe is the cell Peclet number
//   Pe = mu * h / sigma^2
// with mu = r - q - sigma^2/2 (log-space drift) and h = Delta x.
//
// cf. [Duffy04, §4, Eq. 12-13].
//
// The function is even (f(-x) = f(x)) and always >= 1, which
// means fitting never reduces the base diffusion — it can only
// increase it, improving positivity.
//
// Three evaluation regimes avoid numerical hazards:
//   1. |x| < xSmall  — Taylor series (avoids 0/0 in x/tanh(x))
//   2. xSmall <= |x| <= xLarge — direct std::tanh
//   3. |x| > xLarge  — asymptotic |x| (avoids overflow in
//                       exp(2x) inside tanh for large |x|)
// ══════════════════════════════════════════════════════════════════

/*! \file fdmhyperboliccot.hpp
    \brief numerically stable evaluation of x * coth(x)
*/

#ifndef quantlib_fdm_hyperbolic_cot_hpp
#define quantlib_fdm_hyperbolic_cot_hpp

#include <ql/types.hpp>
#include <cmath>

namespace QuantLib {
namespace detail {

    /*! Evaluates \f$ f(x) = x \coth(x) = x / \tanh(x) \f$ in a
        numerically stable way across the full real line.

        - For \f$|x| < x_{\text{small}}\f$, uses the Taylor series
          \f$ 1 + x^2/3 + O(x^4) \f$.
        - For \f$|x| > x_{\text{large}}\f$, uses the asymptotic
          approximation \f$ |x| \f$.
        - Otherwise computes directly via \c std::tanh.

        The function is even and always \f$\ge 1\f$.
    */
    inline Real xCothx(Real x, Real xSmall, Real xLarge) {
        const Real ax = std::fabs(x);

        if (ax < xSmall) {
            // Taylor: x*coth(x) = 1 + x^2/3 - x^4/45 + ...
            // Two-term truncation; relative error < (x^4/45)/(1) ~ 2e-25
            // when xSmall = 1e-6.  [Duffy04, §4]
            const Real x2 = x * x;
            return 1.0 + x2 / 3.0;
        }
        if (ax > xLarge) {
            // For |x| >> 1, tanh(x) -> sign(x), so
            // x * coth(x) = x / tanh(x) -> |x|.
            // At xLarge = 50 the relative error is ~exp(-100) ~ 0.
            return ax;
        }
        // General case: tanh is well-behaved for moderate |x|.
        return x / std::tanh(x);
    }

} // namespace detail
} // namespace QuantLib

#endif
