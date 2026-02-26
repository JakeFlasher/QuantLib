// r6
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#include <ql/instruments/payoffs.hpp>
#include <ql/math/functional.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesop.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/operators/secondderivativeop.hpp>
#include <ql/methods/finitedifferences/operators/fdmhyperboliccot.hpp>
#include <ql/methods/finitedifferences/operators/fdmmatrixdiagnostic.hpp>
#include <ql/methods/finitedifferences/operators/modtriplebandlinearop.hpp>
#include <utility>
#include <vector>
#include <algorithm>
#include <cmath>

namespace QuantLib {

    namespace {

        Array computeEffectiveDiffusion(
                const Array& v,
                const Array& drift,
                const Array& h,
                const std::vector<bool>& isInterior,
                const FdmBlackScholesSpatialDesc& desc,
                FdmBlackScholesSpatialDesc::Scheme scheme,
                Rate r) {

            const Size n = v.size();
            Array aUsed(n);

            for (Size i = 0; i < n; ++i) {
                const Real vEff = std::max(v[i], desc.minVariance);
                const Real aBase = 0.5 * vEff;

                if (!isInterior[i]) {
                    aUsed[i] = aBase;
                    continue;
                }

                // Defensive: skip fitting when mesh spacing is
                // non-positive (should not happen on valid meshers).
                if (!(h[i] > 0.0)) {
                    aUsed[i] = aBase;
                    continue;
                }

                switch (scheme) {
                  case FdmBlackScholesSpatialDesc::Scheme::ExponentialFitting: {
                    const Real Pe = drift[i] * h[i] / vEff;
                    const Real rho = detail::xCothx(Pe, desc.peSmall,
                                                         desc.peLarge);
                    aUsed[i] = aBase * rho;
                    break;
                  }
                  case FdmBlackScholesSpatialDesc::Scheme::MilevTaglianiCNEffectiveDiffusion: {
                    Real aAdd = r * r * h[i] * h[i] / (8.0 * vEff);
                    aAdd = std::min(aAdd, desc.maxAddedDiffusionRatio * aBase);
                    aUsed[i] = aBase + aAdd;
                    break;
                  }
                  default:
                    aUsed[i] = aBase;
                    break;
                }
            }
            return aUsed;
        }

    } // anonymous namespace

    FdmBlackScholesOp::FdmBlackScholesOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& bsProcess,
        Real strike,
        bool localVol,
        Real illegalLocalVolOverwrite,
        Size direction,
        ext::shared_ptr<FdmQuantoHelper> quantoHelper,
        FdmBlackScholesSpatialDesc spatialDesc)
    : mesher_(mesher), rTS_(bsProcess->riskFreeRate().currentLink()),
      qTS_(bsProcess->dividendYield().currentLink()),
      volTS_(bsProcess->blackVolatility().currentLink()),
      localVol_((localVol) ? bsProcess->localVolatility().currentLink()
                           : ext::shared_ptr<LocalVolTermStructure>()),
      x_((localVol) ? Array(Exp(mesher->locations(direction))) : Array()),
      dxMap_(FirstDerivativeOp(direction, mesher)),
      dxxMap_(SecondDerivativeOp(direction, mesher)),
      mapT_(direction, mesher),
      strike_(strike),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      direction_(direction),
      quantoHelper_(std::move(quantoHelper)),
      spatialDesc_(spatialDesc) {}

    void FdmBlackScholesOp::setTime(Time t1, Time t2) {
        const Rate r = rTS_->forwardRate(t1, t2, Continuous).rate();
        const Rate q = qTS_->forwardRate(t1, t2, Continuous).rate();

        // ---- StandardCentral: original code path, verbatim ----
        if (spatialDesc_.scheme
                == FdmBlackScholesSpatialDesc::Scheme::StandardCentral) {

            if (localVol_ != nullptr) {
                Array v(mesher_->layout()->size());
                for (const auto& iter : *mesher_->layout()) {
                    const Size i = iter.index();
                    if (illegalLocalVolOverwrite_ < 0.0) {
                        v[i] = squared(
                            localVol_->localVol(0.5*(t1+t2), x_[i], true));
                    } else {
                        try {
                            v[i] = squared(
                                localVol_->localVol(0.5*(t1+t2), x_[i], true));
                        } catch (Error&) {
                            v[i] = squared(illegalLocalVolOverwrite_);
                        }
                    }
                }
                if (quantoHelper_ != nullptr) {
                    mapT_.axpyb(r - q - 0.5*v
                        - quantoHelper_->quantoAdjustment(Sqrt(v), t1, t2),
                        dxMap_, dxxMap_.mult(0.5*v), Array(1, -r));
                } else {
                    mapT_.axpyb(r - q - 0.5*v, dxMap_,
                                dxxMap_.mult(0.5*v), Array(1, -r));
                }
            } else {
                const Real v
                    = volTS_->blackForwardVariance(t1, t2, strike_)/(t2-t1);
                if (quantoHelper_ != nullptr) {
                    mapT_.axpyb(
                        Array(1, r - q - 0.5*v)
                            - quantoHelper_->quantoAdjustment(
                                Array(1, std::sqrt(v)), t1, t2),
                        dxMap_,
                        dxxMap_.mult(0.5*Array(mesher_->layout()->size(), v)),
                        Array(1, -r));
                } else {
                    mapT_.axpyb(Array(1, r - q - 0.5*v), dxMap_,
                        dxxMap_.mult(0.5*Array(mesher_->layout()->size(), v)),
                        Array(1, -r));
                }
            }
            return;  // done — no diagnostics needed for baseline
        }

        // ---- Non-standard scheme paths (ExponentialFitting / Scheme-2) ----
        const Size n = mesher_->layout()->size();

        // 1. Per-node variance
        Array v(n);
        if (localVol_ != nullptr) {
            for (const auto& iter : *mesher_->layout()) {
                const Size i = iter.index();
                if (illegalLocalVolOverwrite_ < 0.0) {
                    v[i] = squared(
                        localVol_->localVol(0.5*(t1+t2), x_[i], true));
                } else {
                    try {
                        v[i] = squared(
                            localVol_->localVol(0.5*(t1+t2), x_[i], true));
                    } catch (Error&) {
                        v[i] = squared(illegalLocalVolOverwrite_);
                    }
                }
            }
        } else {
            const Real vScalar =
                volTS_->blackForwardVariance(t1, t2, strike_) / (t2 - t1);
            std::fill(v.begin(), v.end(), vScalar);
        }

        // 2. Per-node drift (including quanto adjustment)
        Array drift(n);
        if (quantoHelper_ != nullptr) {
            drift = (r - q) - 0.5 * v
                    - quantoHelper_->quantoAdjustment(Sqrt(v), t1, t2);
        } else {
            drift = (r - q) - 0.5 * v;
        }

        // 3. Per-node mesh spacing (interior only, boundary-safe)
        Array h(n, 0.0);
        std::vector<bool> isInterior(n, false);
        const Size dimDir = mesher_->layout()->dim()[direction_];

        for (const auto& iter : *mesher_->layout()) {
            const Size i  = iter.index();
            const Size co = iter.coordinates()[direction_];

            if (co > 0 && co < dimDir - 1) {
                isInterior[i] = true;
                const Real hm = mesher_->dminus(iter, direction_);
                const Real hp = mesher_->dplus(iter, direction_);

                switch (spatialDesc_.hPolicy) {
                  case FdmBlackScholesSpatialDesc::HPolicy::Average:
                    h[i] = 0.5 * (hm + hp);
                    break;
                  case FdmBlackScholesSpatialDesc::HPolicy::Min:
                    h[i] = std::min(hm, hp);
                    break;
                  case FdmBlackScholesSpatialDesc::HPolicy::Harmonic:
                    h[i] = (hm > 0.0 && hp > 0.0)
                                ? 2.0 * hm * hp / (hm + hp)
                                : 0.5 * (hm + hp);
                    break;
                }
            }
        }

        // 4. Compute modified diffusion coefficient array
        Array aUsed = computeEffectiveDiffusion(
            v, drift, h, isInterior, spatialDesc_,
            spatialDesc_.scheme, r);

        // 5. Assemble tridiagonal operator
        mapT_.axpyb(drift, dxMap_, dxxMap_.mult(aUsed), Array(1, -r));

        // 6. Diagnostics & fallback
        if (spatialDesc_.mMatrixPolicy
                != FdmBlackScholesSpatialDesc::MMatrixPolicy::None) {

            ModTripleBandLinearOp probe(mapT_);
            FdmMMatrixReport report = checkOffDiagonalNonNegative(
                probe, mesher_, direction_,
                spatialDesc_.checkBoundaries, spatialDesc_.mMatrixEps);

            if (!report.ok) {
                using Policy = FdmBlackScholesSpatialDesc::MMatrixPolicy;
                using Scheme = FdmBlackScholesSpatialDesc::Scheme;

                switch (spatialDesc_.mMatrixPolicy) {
                  case Policy::FailFast:
                    QL_REQUIRE(false,
                        "M-matrix violation: " << report.summary());
                    break;

                  case Policy::FallbackToExponentialFitting:
                    if (spatialDesc_.scheme != Scheme::ExponentialFitting) {
                        aUsed = computeEffectiveDiffusion(
                            v, drift, h, isInterior, spatialDesc_,
                            Scheme::ExponentialFitting, r);
                        mapT_.axpyb(drift, dxMap_,
                                    dxxMap_.mult(aUsed), Array(1, -r));
                    }
                    break;

                  case Policy::DiagnosticsOnly:
                    break;

                  default:
                    break;
                }
            }
        }
    }

    Size FdmBlackScholesOp::size() const { return 1U; }

    Array FdmBlackScholesOp::apply(const Array& u) const {
        return mapT_.apply(u);
    }

    Array FdmBlackScholesOp::apply_direction(Size direction,
                                             const Array& r) const {
        if (direction == direction_)
            return mapT_.apply(r);
        else {
            return Array(r.size(), 0.0);
        }
    }

    Array FdmBlackScholesOp::apply_mixed(const Array& r) const {
        return Array(r.size(), 0.0);
    }

    Array FdmBlackScholesOp::solve_splitting(Size direction,
                                             const Array& r, Real dt) const {
        if (direction == direction_)
            return mapT_.solve_splitting(r, dt, 1.0);
        else {
            return r;
        }
    }

    Array FdmBlackScholesOp::preconditioner(const Array& r,
                                            Real dt) const {
        return solve_splitting(direction_, r, dt);
    }

    std::vector<SparseMatrix> FdmBlackScholesOp::toMatrixDecomp() const {
        return std::vector<SparseMatrix>(1, mapT_.toMatrix());
    }

}
