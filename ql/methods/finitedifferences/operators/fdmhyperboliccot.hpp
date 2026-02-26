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
            const Real x2 = x * x;
            return 1.0 + x2 / 3.0;
        }
        if (ax > xLarge) {
            // coth(x) -> sign(x), so x*coth(x) -> |x|
            return ax;
        }
        // General case: stable for moderate arguments
        return x / std::tanh(x);
    }

} // namespace detail
} // namespace QuantLib

#endif
