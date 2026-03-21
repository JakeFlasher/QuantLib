# Comprehensive Code Audit: Nonstandard FD Schemes in QuantLib

## Methodology

This audit was conducted by:
1. Exhaustive exploration of all 22 r6-tagged source files
2. Sharing the three reference papers (Milev-Tagliani 2010 Nonstandard FD, Duffy 2004, Milev-Tagliani 2010 Low Volatility) with gpt-5.4:xhigh via Codex
3. Iterative Codex queries on each component (operator, mesher, step condition, solver, engine) with full code snippets
4. Manual verification of each identified pitfall against the actual codebase

---

## CONFIRMED FINDINGS

### Finding 1: MT Scheme Missing Drift Correction in Log-Space [SEVERITY: HIGH]

**File:** `ql/methods/finitedifferences/operators/fdmblackscholesop.cpp:58-63`

**Description:** The Milev-Tagliani CN-variant paper derives its modified equation in S-space. The artificial numerical diffusion term is `(1/8)(r*dS/sigma)^2 * V_SS`. When this is transformed to log-space `x = ln(S)`, the chain rule gives:

```
V_SS = S^{-2} (V_xx - V_x)
```

So the added diffusion becomes `aAdd * (V_xx - V_x)`, meaning there should be BOTH:
- An added diffusion term: `+aAdd` on the `V_xx` coefficient
- An added drift correction: `-aAdd` on the `V_x` coefficient

The current code ONLY modifies the diffusion coefficient:
```cpp
case Scheme::MilevTaglianiCNEffectiveDiffusion: {
    Real aAdd = r * r * h[i] * h[i] / (8.0 * vEff);
    aAdd = std::min(aAdd, desc.maxAddedDiffusionRatio * aBase);
    aUsed[i] = aBase + aAdd;  // Only diffusion modified
    break;
}
```

And uses the unmodified drift in `axpyb`:
```cpp
mapT_.axpyb(drift, dxMap_, dxxMap_.mult(aUsed), Array(1, -r));
// drift is NOT corrected by -aAdd
```

**Impact:** The missing drift correction means the MT scheme in log-space is NOT a faithful translation of the paper's S-space scheme. The operator solves a slightly different modified equation than intended, which can lead to systematic pricing bias, especially when `aAdd` is large (low volatility, large grid spacing, high interest rate).

**Test Scenario:**
```
Parameters: S=100, K=100, r=0.5, sigma=0.001, T=5/12, q=0
Grid: xGrid=200, tGrid=100
Scheme: MilevTaglianiCNEffectiveDiffusion

With these extreme parameters (sigma^2 << r), aAdd = r^2*h^2/(8*sigma^2)
becomes very large relative to aBase = sigma^2/2 = 5e-7.

Compare:
  (a) Current implementation (diffusion-only correction)
  (b) Full correction (diffusion + drift)
  (c) Analytical Black-Scholes price

Expected: (a) shows systematic bias vs analytical; (b) reduces bias.
The drift correction term -aAdd would be comparable to the original drift
r-q-sigma^2/2 = 0.4999995, creating a visible price difference.
```

---

### Finding 2: thetaAt() Missing calculate() Call [SEVERITY: HIGH]

**File:** `ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp:97-99`

**Description:** All other accessor methods (`valueAt`, `deltaAt`, `gammaAt`) call `calculate()` before accessing `solver_`. But `thetaAt()` does not:

```cpp
Real FdmBlackScholesSolver::thetaAt(Real s) const {
    return solver_->thetaAt(std::log(s));  // No calculate()!
}
```

This is a classic LazyObject bug. If `thetaAt()` is called before any other accessor, `solver_` is null and this crashes. Even if called after another accessor, if the process/curves have changed, `solver_` may hold stale results.

**Test Scenario:**
```cpp
// Create a solver but only access theta:
auto solver = ext::make_shared<FdmBlackScholesSolver>(...);
Real theta = solver->thetaAt(100.0);  // CRASH: null pointer dereference
```

---

### Finding 3: Missing European-Only Guard in calculateDiscrete() [SEVERITY: MEDIUM]

**File:** `ql/pricingengines/barrier/fdblackscholesbarrierengine.cpp:238`

**Description:** `calculateContinuous()` (line 104) explicitly requires European exercise. `calculateDiscrete()` has no such guard and does NOT add American/Bermudan step conditions. A non-European discrete barrier option would silently compute European-style pricing.

**Test Scenario:**
```cpp
// Create an American barrier option with discrete monitoring:
auto exercise = ext::make_shared<AmericanExercise>(today, maturityDate);
BarrierOption option(Barrier::DownOut, 90.0, 0.0, payoff, exercise);
option.setPricingEngine(
    ext::make_shared<FdBlackScholesBarrierEngine>(
        process, monitoringDates, tGrid, xGrid, 0,
        FdmSchemeDesc::CrankNicolson()));
Real price = option.NPV();  // Silently computes European, not American
```

---

### Finding 4: Inconsistent Boundary Filtering in Multi-Point Mesher [SEVERITY: MEDIUM]

**File:** `ql/methods/finitedifferences/meshers/fdmblackscholesmesher.cpp:192`

**Description:** The single-point constructor uses inclusive comparisons (`>= xMin && <= xMax`), but the multi-point constructor uses strict comparisons (`> xMin && < xMax`). This means barrier critical points exactly at the domain boundary are silently dropped in multi-point mode.

```cpp
// Single-point (line 137-138): INCLUSIVE
if (std::log(cPoint.first) >= xMin && std::log(cPoint.first) <= xMax)

// Multi-point (line 192): STRICT
if (logLevel > xMin && logLevel < xMax)
```

**Test Scenario:**
```
Set xMinConstraint = log(90) for a barrier at S=90.
Single-point constructor: cPoint at S=90 is included (log(90) >= log(90)).
Multi-point constructor: cPoint at S=90 is DROPPED (log(90) > log(90) is false).
Result: barrier pricing with multi-point mesher loses mesh concentration
at the barrier, degrading accuracy near the barrier.
```

---

### Finding 5: CN-Equivalence Gate Too Narrow and Too Permissive [SEVERITY: MEDIUM]

**File:** `ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp:33-44`

**Description:** The gate has two problems:

**(a) Too narrow:** In 1D, `CraigSneydType` with `theta=0.5` also reduces to CN because `apply_mixed()` returns zero. The gate only accepts `CrankNicolsonType` and `DouglasType`.

**(b) Too permissive:** The gate ignores `dampingSteps`. When `dampingSteps > 0`, `FdmBackwardSolver` prepends Implicit Euler damping steps (fdmbackwardsolver.cpp:100-106), so the actual time stepping is NOT purely CN even when `schemeDesc` says CN.

**Test Scenario:**
```cpp
// (a) CraigSneyd in 1D wrongly rejected:
FdmSchemeDesc desc = FdmSchemeDesc::CraigSneyd(); // theta=0.5
FdmBlackScholesSolver solver(..., desc, ..., spatialDesc_mt);
// MT silently falls back to ExponentialFitting despite being CN-equivalent in 1D

// (b) CN with damping wrongly accepted:
FdmSchemeDesc desc = FdmSchemeDesc::CrankNicolson();
// dampingSteps=3 in solverDesc
FdmBlackScholesSolver solver(..., desc, ..., spatialDesc_mt);
// MT is accepted but first 3 steps are Implicit Euler, not CN
```

---

### Finding 6: Single-Barrier Engine Used for Double-Barrier Paper [SEVERITY: LOW-MEDIUM]

**File:** `ql/pricingengines/barrier/fdblackscholesbarrierengine.cpp:311-319`

**Description:** The discrete-monitoring path treats single-barrier instruments by faking the missing barrier with `1e-15`/`1e15` sentinels. The Milev-Tagliani papers specifically address double-barrier knock-out options. For a true double barrier, the engine would need both barriers from the instrument arguments, but `BarrierOption::arguments` only has one barrier.

**Impact:** This is a design limitation, not a bug per se. The engine works correctly for single-barrier discrete monitoring. Double-barrier requires a separate `DoubleBarrierOption` instrument.

**Test Scenario:**
```
For a double-barrier knock-out (L=90, U=110):
- Current engine cannot represent this as a BarrierOption (only one barrier)
- User must use test-suite code that manually constructs FdmDiscreteBarrierStepCondition
  with both barriers, bypassing the engine
- The paper's primary example is double-barrier, so the engine is not paper-faithful
```

---

### Finding 7: Floating-Point binary_search on Monitoring Times [SEVERITY: LOW]

**File:** `ql/methods/finitedifferences/stepconditions/fdmdiscretebarrierstepcondition.cpp:66-67`

**Description:** The `applyTo` method uses exact `binary_search` on floating-point `Time` values. Within the current pipeline this works because the same time values are passed to the FD stopping-time scheduler and the step condition. However, this is brittle as a reusable component.

```cpp
if (!std::binary_search(monitoringTimes_.begin(),
                        monitoringTimes_.end(), t))
    return;
```

**Test Scenario:**
```
If monitoring times are recomputed in a slightly different way (e.g.,
different day-count convention, or different DayCounter rounding),
the binary_search can fail to match:

Time t_scheduled = 0.50000000000001;  // from solver
Time t_stored    = 0.49999999999999;  // in monitoringTimes_
binary_search fails -> barrier is never applied at this monitoring date
```

---

### Finding 8: No MT Time-Step Positivity Restriction Enforced [SEVERITY: LOW-MEDIUM]

**File:** `ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp`

**Description:** The MT paper (Theorem 3.2) requires `dt < 1/(r*M)` for the scheme to have real positive distinct eigenvalues. The code enforces no such restriction. The M-matrix diagnostic in the operator only checks off-diagonal signs of the spatial operator, NOT the full time-discrete positivity condition from the paper.

**Test Scenario:**
```
Parameters: r=0.5, M=200 (grid nodes), dt = maturity/tGrid
Paper requires: dt < 1/(0.5 * 200) = 0.01
If tGrid=10, maturity=1.0, then dt=0.1 >> 0.01
Result: Eigenvalues of iteration matrix may not be in (0,1),
potentially causing oscillatory behavior despite using MT scheme.
```

---

### Finding 9: Discrete Barrier t=0 Monitoring Omitted [SEVERITY: LOW]

**File:** `ql/pricingengines/barrier/fdblackscholesbarrierengine.cpp:306`

**Description:** The discrete path filters monitoring times with `t > 0.0 && t <= maturity`, deliberately excluding t=0. The paper's formulation uses `0 = t_0 < t_1 < ... < t_F = T` with the initial condition already multiplied by the indicator function.

**Impact:** If the spot starts outside [L,U], the engine does not immediately knock out, which is inconsistent with paper equation (5). In practice, for standard barrier option setups, the spot is usually inside the corridor.

**Test Scenario:**
```
S0=85 (below barrier L=90), monitoring dates include today
Paper: option should be worth 0 (knocked out at t=0)
Code: option returns non-zero value because t=0 monitoring is skipped
```

---

## VERIFIED NON-ISSUES / CORRECTLY IMPLEMENTED

1. **Peclet number computation (ExponentialFitting):** `Pe = drift*h/vEff` where drift = r-q-sigma^2/2 and vEff = sigma^2. This correctly maps Duffy's formula to log-space.

2. **xCothx numerical stability:** Taylor truncation at xSmall=1e-6 with 2 terms is more than adequate (omitted term ~1e-26). Large-x cutoff at 50 is conservative but correct.

3. **Gamma formula:** `(V_xx - V_x)/(s*s)` correctly converts log-space second derivative to S-space gamma.

4. **ExponentialFitting scheme:** Correctly implements Duffy's fitting factor as `aUsed = aBase * rho` where `rho = xCothx(Pe, ...)`.

5. **FallbackToExponentialFitting safety net:** Conceptually sound, though could be more complete (re-check after fallback).

6. **Knock-in parity (In = Vanilla - Out):** Mathematically correct for European exercise with zero rebate.
