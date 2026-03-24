// ══════════════════════════════════════════════════════════════════
// FdBlackScholesBarrierEngine — 1-D FDM barrier-option engine
//
// Supports two monitoring modes:
//   Continuous  — Dirichlet BCs at barrier level(s)
//   Discrete    — FdmDiscreteBarrierStepCondition at monitoring dates
//
// The spatial descriptor (FdmBlackScholesSpatialDesc) is threaded
// through to the solver/operator, allowing the caller to select
// StandardCentral, ExponentialFitting, or MilevTaglianiCN schemes.
//
// Knock-in pricing uses parity: V_in = V_vanilla - V_out.
// cf. [Ballabio20, Ch. 11] for the FDM barrier framework.
// ══════════════════════════════════════════════════════════════════

// r6
/*! \file fdblackscholesbarrierengine.hpp
    \brief Finite-differences Black/Scholes barrier-option engine
*/

#ifndef quantlib_fd_black_scholes_barrier_engine_hpp
#define quantlib_fd_black_scholes_barrier_engine_hpp

#include <ql/processes/blackscholesprocess.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp>
#include <ql/instruments/barrieroption.hpp>
#include <vector>

namespace QuantLib {

    class FdmBlackScholesSolver;

    //! Finite-differences barrier engine (single-barrier only).
    //! For double-barrier knock-out options (Milev-Tagliani Example 4.1),
    //! use FdmDiscreteBarrierStepCondition directly with both barriers,
    //! or a dedicated double-barrier engine such as FdHestonDoubleBarrierEngine.
    class FdBlackScholesBarrierEngine : public BarrierOption::engine {
      public:
        // --- Continuous monitoring constructors (unchanged API) ---

        explicit FdBlackScholesBarrierEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        explicit FdBlackScholesBarrierEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            DividendSchedule dividends,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        // --- Discrete monitoring constructors ---

        FdBlackScholesBarrierEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            std::vector<Date> monitoringDates,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        FdBlackScholesBarrierEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            DividendSchedule dividends,
            std::vector<Date> monitoringDates,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            FdmBlackScholesSpatialDesc spatialDesc
                = FdmBlackScholesSpatialDesc());

        void calculate() const override;

      private:
        void calculateContinuous() const;
        void calculateDiscrete() const;

        //! Populate the 4 fallback-observability additionalResults keys
        //  by aggregating state from the main solver and any sub-option
        //  instruments (vanilla, rebate) that contributed to the final
        //  result.  If any contributing solve triggered a fallback, the
        //  top-level engine reports it.
        void reportSpatialScheme(
            const ext::shared_ptr<FdmBlackScholesSolver>& mainSolver,
            const std::vector<const Instrument*>& subInstruments
                = std::vector<const Instrument*>()) const;

        ext::shared_ptr<GeneralizedBlackScholesProcess> process_;
        DividendSchedule dividends_;
        std::vector<Date> monitoringDates_;
        Size tGrid_, xGrid_, dampingSteps_;
        FdmSchemeDesc schemeDesc_;
        bool localVol_;
        Real illegalLocalVolOverwrite_;
        FdmBlackScholesSpatialDesc spatialDesc_;
    };

}

#endif
