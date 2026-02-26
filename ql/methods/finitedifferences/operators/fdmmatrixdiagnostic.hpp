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

    struct FdmMMatrixReport {
        bool ok = true;

        Size checkedRows   = 0;
        Size negativeLower = 0;
        Size negativeUpper = 0;

        Real minLower = 0.0;
        Real minUpper = 0.0;

        Size idxMinLower = Null<Size>();
        Size idxMinUpper = Null<Size>();

        bool hasNonFinite  = false;
        Size nonFiniteCount = 0;

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

            if (!includeBoundaries && (co == 0 || co == dimDir - 1))
                continue;

            report.checkedRows++;

            const Real lo = op.lower(i);
            const Real di = op.diag(i);
            const Real up = op.upper(i);

            if (!std::isfinite(lo) || !std::isfinite(di) || !std::isfinite(up)) {
                report.hasNonFinite = true;
                report.nonFiniteCount++;
            }

            if (lo < -eps) {
                report.negativeLower++;
                if (lo < report.minLower) {
                    report.minLower    = lo;
                    report.idxMinLower = i;
                }
            }

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
