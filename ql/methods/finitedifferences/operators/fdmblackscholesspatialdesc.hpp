/*! \file fdmblackscholesspatialdesc.hpp
    \brief spatial discretization descriptor for 1D Black-Scholes FDM operators
*/

#ifndef quantlib_fdm_black_scholes_spatial_desc_hpp
#define quantlib_fdm_black_scholes_spatial_desc_hpp

#include <ql/types.hpp>

namespace QuantLib {

    struct FdmBlackScholesSpatialDesc {

        enum class Scheme {
            StandardCentral,
            ExponentialFitting,
            MilevTaglianiCNEffectiveDiffusion
        };

        enum class HPolicy {
            Average,
            Min,
            Harmonic
        };

        enum class BoundaryPolicy {
            InteriorOnly
        };

        enum class MMatrixPolicy {
            None,
            DiagnosticsOnly,
            FailFast,
            FallbackToExponentialFitting
        };

        Scheme scheme = Scheme::StandardCentral;

        HPolicy hPolicy = HPolicy::Average;
        BoundaryPolicy boundaryPolicy = BoundaryPolicy::InteriorOnly;

        Real peSmall  = 1e-6;
        Real peLarge  = 50.0;

        Real minVariance           = 1e-12;
        Real maxAddedDiffusionRatio = 1e6;

        MMatrixPolicy mMatrixPolicy = MMatrixPolicy::FallbackToExponentialFitting;
        Real mMatrixEps     = 0.0;
        bool checkBoundaries = false;

        static FdmBlackScholesSpatialDesc standard() {
            return {};
        }

        static FdmBlackScholesSpatialDesc exponentialFitting() {
            FdmBlackScholesSpatialDesc d;
            d.scheme = Scheme::ExponentialFitting;
            return d;
        }

        static FdmBlackScholesSpatialDesc milevTaglianiCN() {
            FdmBlackScholesSpatialDesc d;
            d.scheme = Scheme::MilevTaglianiCNEffectiveDiffusion;
            return d;
        }
    };

}

#endif
