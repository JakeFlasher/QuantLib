# Future Work Roadmap — Nonstandard FD Schemes for QuantLib

Based on review of the Milev-Tagliani and Duffy papers, Codex (GPT-5.4 xhigh) review,
and the current r6 implementation of exponential fitting + CN-variant effective diffusion.

Suggested order: 1 → 2 → 3 → 4 → 6, with 5 pulled forward if runtime is a pain point.

---

## 1. Generalized Event Engine (American/Bermudan + Time-Dependent Barriers)

**What:** Refactor the discrete-barrier path into a reusable event layer that can combine
monitoring dates, barrier schedules B(t_i), and early-exercise conditions, with an explicit
same-date ordering rule (e.g., exercise checked before barrier).

**Why:** Opens callable/barrier hybrids, Bermudan knock-outs, step-up barriers, and window
barriers — common in structured products desks.

**Complexity:** Medium

**Leverage:** `FdmDiscreteBarrierStepCondition`, `FdmStepConditionComposite`,
`FdmAmericanStepCondition`, `FdmBermudanStepCondition`, `FdmTimeDepDirichletBoundary`,
and the current barrier-engine builder.

---

## 2. Event-Aware Accuracy Layer (Rannacher + Richardson)

**What:** Extend existing `dampingSteps` so implicit-Euler/Rannacher damping can restart
after each discontinuity injection (monitoring date), not just at rollback start. Add an
optional two-grid Richardson extrapolation wrapper for O(dt^2) recovery.

**Why:** Provides a clean benchmark against the Milev-Tagliani/Duffy fixes and improves
post-monitoring Greeks stability with minimal model risk. Rannacher time-stepping is an
established alternative to CN-variant schemes for handling non-smooth data.

**Complexity:** Medium

**Leverage:** `FdmBackwardSolver`, `FdmSolverDesc`, `FdmStepConditionComposite::stoppingTimes`,
`FdmSchemeDesc`, and `MakeFdBlackScholesBarrierEngine`.

---

## 3. Local-Vol Surface Hardening

**What:** Make the new spatial schemes first-class under calibrated local-vol surfaces, with
explicit tests on interpolation/extrapolation edge cases and barrier-event stress scenarios.

**Why:** Desks usually have a smile surface before a full stochastic-vol calibration, so
this is the lowest-risk path to smile-consistent barrier pricing. The node-wise variance
`v[i]` already supports local vol — but it needs testing under the modified schemes.

**Complexity:** Medium

**Leverage:** The `localVol` path in `FdBlackScholesBarrierEngine`, node-wise local-vol
assembly in `FdmBlackScholesOp`, `LocalVolSurface`, and `MMatrixPolicy`/diagnostics.

---

## 4. 2D/ND Multi-Asset Barrier Engines via ADI

**What:** Build a 2D barrier engine first (worst-of/best-of, min/max barriers, dual-asset
no-touch), then generalize to N dimensions. The exponential fitting and effective-diffusion
ideas extend naturally to each spatial dimension in ADI splitting.

**Why:** Biggest product-coverage jump. Multi-asset barriers are commonly traded but lack
robust FD support in QuantLib.

**Complexity:** High

**Leverage:** `Fdm2dBlackScholesOp`, `Fdm2dBlackScholesSolver`, `Fdm2DimSolver`,
`FdmNdimSolver`, and ADI schemes already exposed through `FdmSchemeDesc` (Douglas, Craig-Sneyd,
Hundsdorfer-Verwer, Modified Craig-Sneyd).

---

## 5. Adaptive Mesh Refinement Near Discontinuities

**What:** Move from exact strike/barrier alignment (current multi-point mesher) to adaptive
density or remeshing around renewed kinks after monitoring/exercise dates. Use M-matrix
diagnostics as refinement/error indicators.

**Why:** Runtime savings become material before going 2D/Heston. This is the cleanest way
to reduce numerical diffusion without always raising global `xGrid`.

**Complexity:** High

**Leverage:** `Concentrating1dMesher`, multi-point `FdmBlackScholesMesher`,
`FdmMesherComposite`, and M-matrix diagnostic infrastructure as refinement triggers.

---

## 6. Stochastic-Vol Barrier Roadmap (Heston → SLV/SABR)

**What:** Port the same barrier/event layer to stochastic-vol PDEs, starting with Heston
(2D PDE). Once stable, use leverage-function hooks for SLV and the existing SABR PDE engine
as a second 2D testbed.

**Why:** Skew dynamics matter significantly for knock-out probabilities and long-dated exotics.
A flat-vol barrier engine misses the smile, which can easily be 10-30% of the price.

**Complexity:** High

**Leverage:** `FdHestonBarrierEngine`, `FdmHestonOp` (including leverage-function support),
`FdSabrVanillaEngine`, `FdmSabrOp`, and the same 2D ADI solver stack from item 4.

---

## Notes

- **GPU/SIMD:** Deferred until items 4/6 create enough compute pressure. Obvious hot spots
  are per-node coefficient assembly in `FdmBlackScholesOp::setTime` and
  `Fdm2dBlackScholesOp::setTime`. Payoff is much larger once batching 2D/stoch-vol solves.

- **The semi-implicit scheme (Paper 1, Section 3.2)** with b = -M/2 and O(dS^2, dt)
  accuracy is NOT implemented. It would require a custom time-stepping scheme rather than
  a spatial-operator modification. The CN variant (Paper 3) was chosen instead because it
  preserves O(dt^2) accuracy and integrates cleanly with QuantLib's existing time-stepping.

- **Non-uniform mesh h-policy:** Codex noted that for monotonicity on non-uniform grids,
  an upwind-side / max-spacing choice for h might be more defensible than plain averaging.
  This deserves investigation if adaptive meshes (item 5) are pursued.
