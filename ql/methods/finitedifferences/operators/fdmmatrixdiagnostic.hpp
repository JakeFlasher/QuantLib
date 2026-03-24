// ══════════════════════════════════════════════════════════════════
// FdmMMatrixReport / checkOffDiagonalNonNegative
//
// M-matrix diagnostic for tridiagonal finite-difference operators.
//
// A tridiagonal matrix A is an M-matrix (Z-matrix with non-negative
// inverse) if and only if all off-diagonal entries are non-positive
// in A, i.e. non-negative in -A (the negated convention used by
// QuantLib's operator representation).  Positivity of the solution
// vector under explicit or Crank-Nicolson time stepping requires
// this property — cf. [MT10, Thm 3.1], [Duffy04, §3].
//
// checkOffDiagonalNonNegative() scans the lower and upper bands of
// a ModTripleBandLinearOp and reports how many interior rows (and
// optionally boundary rows) violate the non-negativity condition.
//
// The caller (FdmBlackScholesOp::setTime) uses this report to
// decide on corrective action according to MMatrixPolicy:
//   - None:                         skip the check entirely
//   - FailFast:                     abort on first violation
//   - FallbackToExponentialFitting: silently recompute diffusion
//                                   using exponential fitting
// ══════════════════════════════════════════════════════════════════

/*! \file fdmmatrixdiagnostic.hpp
    \brief M-matrix off-diagonal sign diagnostics for tridiagonal operators
*/

#ifndef quantlib_fdm_matrix_diagnostic_hpp
#define quantlib_fdm_matrix_diagnostic_hpp

#include <ql/types.hpp>
#include <ql/utilities/null.hpp>
#include <ql/methods/finitedifferences/operators/modtriplebandlinearop.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <cmath>
#include <string>
#include <sstream>

namespace QuantLib {

    //! Diagnostic report produced by checkOffDiagonalNonNegative().
    struct FdmMMatrixReport {
        bool ok = true;                 //!< true if all checked off-diagonals >= -eps

        Size checkedRows   = 0;         //!< number of rows inspected
        Size negativeLower = 0;         //!< rows where lower band < -eps
        Size negativeUpper = 0;         //!< rows where upper band < -eps

        Real minLower = 0.0;            //!< most-negative lower-band value
        Real minUpper = 0.0;            //!< most-negative upper-band value

        Size idxMinLower = Null<Size>();  //!< layout index of worst lower
        Size idxMinUpper = Null<Size>();  //!< layout index of worst upper

        bool hasNonFinite  = false;     //!< any NaN / Inf detected
        Size nonFiniteCount = 0;        //!< count of non-finite entries

        //! Human-readable one-line summary for error messages.
        std::string summary() const {
            std::ostringstream oss;
            oss << "M-matrix check: " << (ok ? "PASS" : "FAIL")
                << " (" << checkedRows << " rows checked";
            if (!ok) {
                oss << ", neg-lower=" << negativeLower
                    << " (min=" << minLower << " @" << idxMinLower << ")"
                    << ", neg-upper=" << negativeUpper
                    << " (min=" << minUpper << " @" << idxMinUpper << ")";
                if (hasNonFinite)
                    oss << ", non-finite=" << nonFiniteCount;
            }
            oss << ")";
            return oss.str();
        }
    };

    /*! Scan the off-diagonal bands of a tridiagonal operator and
        report rows where the lower or upper entry is more negative
        than \p eps.

        \param op               tridiagonal operator to inspect
        \param mesher           FDM mesher (provides layout and dim)
        \param direction        spatial dimension to check
        \param includeBoundaries  when true, check boundary rows too
        \param eps              tolerance — entries in [-eps, 0) are
                                treated as non-negative (default 0)
    */
    inline FdmMMatrixReport checkOffDiagonalNonNegative(
            const ModTripleBandLinearOp& op,
            const ext::shared_ptr<FdmMesher>& mesher,
            Size direction,
            bool includeBoundaries,
            Real eps = 0.0) {

        FdmMMatrixReport report;
        const Size dimDir = mesher->layout()->dim()[direction];

        for (const auto& iter : *mesher->layout()) {
            const Size i  = iter.index();
            const Size co = iter.coordinates()[direction];

            // Skip boundary nodes unless explicitly requested.
            if (!includeBoundaries && (co == 0 || co == dimDir - 1))
                continue;

            report.checkedRows++;

            const Real lo = op.lower(i);
            const Real di = op.diag(i);
            const Real up = op.upper(i);

            // Guard against degenerate operators (NaN / Inf).
            if (!std::isfinite(lo) || !std::isfinite(di) || !std::isfinite(up)) {
                report.hasNonFinite = true;
                report.nonFiniteCount++;
            }

            // Check lower band: must be >= -eps for M-matrix property.
            if (lo < -eps) {
                report.negativeLower++;
                if (lo < report.minLower) {
                    report.minLower    = lo;
                    report.idxMinLower = i;
                }
            }

            // Check upper band: same condition.
            if (up < -eps) {
                report.negativeUpper++;
                if (up < report.minUpper) {
                    report.minUpper    = up;
                    report.idxMinUpper = i;
                }
            }
        }

        report.ok = (report.negativeLower == 0
                     && report.negativeUpper == 0
                     && !report.hasNonFinite);
        return report;
    }

}

#endif
