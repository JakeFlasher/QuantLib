# Results Analysis: Nonstandard Finite Difference Schemes for Black-Scholes Option Pricing in QuantLib

## 1. Introduction

The Crank-Nicolson (CN) finite difference scheme is the most widely used method for solving the Black-Scholes partial differential equation in computational finance. While CN offers second-order accuracy in both space and time, it suffers from a fundamental weakness: **spurious oscillations** near discontinuities in the initial or boundary conditions. These oscillations produce non-physical negative option prices and inaccurate Greeks, particularly in the low-volatility regime where σ² ≪ r.

This analysis evaluates three spatial discretization schemes implemented in QuantLib's 1D Black-Scholes FDM operator:

1. **StandardCentral** — Baseline centered differences with Crank-Nicolson time stepping
2. **ExponentialFitting** (Duffy, 2004) — Implicit scheme with a fitting factor that guarantees M-matrix property
3. **MilevTaglianiCNEffectiveDiffusion** (Milev-Tagliani, 2010) — CN variant with artificial diffusion from a modified reaction-term discretization

## 2. Mathematical Framework

### 2.1 The Black-Scholes PDE

The option price V(S,t) satisfies:

$$-\frac{\partial V}{\partial t} + (r-q)S\frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} - rV = 0$$

In log-space (x = ln S), this becomes:

$$-\frac{\partial V}{\partial t} + \mu\frac{\partial V}{\partial x} + \frac{\sigma^2}{2}\frac{\partial^2 V}{\partial x^2} - rV = 0$$

where μ = r − q − σ²/2 is the risk-neutral drift.

### 2.2 StandardCentral Scheme

Centered differences for the convection and diffusion terms yield the tridiagonal operator:
- Lower off-diagonal: a_{i,i-1} = σ²/(2h²) − μ/(2h)
- Diagonal: a_{i,i} = −σ²/h² − r
- Upper off-diagonal: a_{i,i+1} = σ²/(2h²) + μ/(2h)

When σ² ≪ r, the convection-dominated Péclet number Pe = μh/σ² becomes large, and the off-diagonal entries can become negative, violating the M-matrix property. This causes spurious oscillations.

### 2.3 Exponential Fitting (Duffy)

Duffy's scheme replaces the standard diffusion coefficient σ²/2 with:

$$a_{\text{eff}} = \frac{\sigma^2}{2} \cdot \rho(Pe), \quad \rho(Pe) = \frac{Pe}{\tanh(Pe)} = x\coth(x)$$

where Pe = μh/σ² is the mesh Péclet number. The fitting factor ρ satisfies ρ ≥ 1 and transitions:
- ρ → 1 as Pe → 0 (recovers standard central differences)
- ρ → |Pe| as |Pe| → ∞ (becomes first-order upwind)

For typical parameter regimes, this produces non-negative off-diagonal entries (M-matrix property), yielding oscillation-free, positive solutions. The worst-case convergence bound is uniform: |V(S_j, t_n) − U_j^n| ≤ c(h + k), independent of h, k, and σ, though observed convergence on smooth data is typically O(h²).

### 2.4 Milev-Tagliani CN Variant

The MT scheme discretizes the reaction term −rV using a 6-point stencil with parameters ω₁ = ω₂ = −r/(16σ²). In typical parameter regimes, this introduces sufficient artificial diffusion to restore the M-matrix property (though edge cases such as negative dividend yield or very coarse grids can still violate it). In the code's log-space formulation, the effective diffusion is:

$$a_{\text{eff}} = \frac{\sigma^2}{2} + \frac{r^2 h^2}{8\sigma^2}$$

**Important implementation note:** The shipped QuantLib operator implements a log-space *diffusion-only adaptation* of the MT scheme. The paper's full S-space → log-space translation would also include a drift correction term of −r²h²/(8σ²). This correction is O(h²) and has been verified (via the `testMilevTaglianiDriftCorrectionAudit` test) to be bounded by the grid convergence error. See Section 3 for details.

### 2.5 CN-Equivalence Requirement

The MT scheme is designed for Crank-Nicolson time stepping. The solver enforces this: any time scheme other than `CrankNicolson` (or `Douglas` with θ=0.5 in 1D) triggers a solver-level fallback to ExponentialFitting, independent of the operator's M-matrix diagnostic. With `FailFast` policy, a non-CN time scheme causes the solver to throw.

## 3. Implementation Notes

### 3.1 Log-Space Formulation

All three schemes operate in log-space (x = ln S), where the Black-Scholes PDE has constant coefficients for flat-vol models. The mesh spacing h is in log-space; the physical ΔS varies across the grid as ΔS ≈ S·h.

### 3.2 MT Drift Correction Omission

The paper's full S-space → log-space translation of the MT artificial diffusion yields both:
- Diffusion addition: +r²h²/(8σ²) to the coefficient of ∂²V/∂x²
- Drift correction: −r²h²/(8σ²) to the coefficient of ∂V/∂x

The shipped code applies only the diffusion addition (using `vEff = max(σ², minVariance)` and capping the added diffusion via `maxAddedDiffusionRatio`). The drift correction is O(h²) — the same order as the artificial diffusion itself — and the test `testMilevTaglianiDriftCorrectionAudit` empirically confirms that on the audited mesh (xGrid=401, r=0.50, σ=0.001) the scheme difference / grid error ratio is < 1.0, and the scheme difference shrinks on finer grids. This is an empirical audit, not a general proof.

### 3.3 Fallback Mechanism

When `mMatrixPolicy=FallbackToExponentialFitting` is set, the operator checks off-diagonal signs after assembly. If any interior off-diagonal is negative (M-matrix violation), it automatically falls back to the ExponentialFitting scheme. This safety net is useful for non-CN time schemes or edge-case parameters.

### 3.4 Numerically Stable xCothx

The fitting factor ρ = x·coth(x) is evaluated using a three-regime approach:
- |x| < 10⁻⁶: Taylor series 1 + x²/3
- |x| > 50: asymptotic approximation |x|
- Otherwise: direct computation x/tanh(x)

This avoids the 0/0 indeterminate form at x = 0 and overflow for large arguments.

## 4. Numerical Experiments

### Experiment 1: Truncated Call — Spurious Oscillations (Figures 1–2)

**Parameters:** r = 0.05, σ = 0.001, K = 50, U = 70, T = 5/12 (paper Example 4.1 from Milev-Tagliani, *Low Volatility Options*, 2010).

**Results:** Figure 1 shows severe oscillations in the StandardCentral solution near the upper barrier U = 70, with negative values. A fine-grid (8× refinement) ExponentialFitting solution serves as the reference, confirming the smooth, correct solution shape. Figure 2 provides a side-by-side comparison of all three schemes with the fine-grid reference, plus an inset zoom near U = 70.

### Experiment 2: Discrete Double Barrier — Moderate Volatility (Figure 3)

**Parameters:** K = 100, σ = 0.25, r = 0.05, L = 95, U = 110, T = 0.5, 5 monitoring dates (paper Example 4.1 from Milev-Tagliani, *Nonstandard FD Schemes*, 2010).

**Results:** At moderate volatility (σ² > r), all three schemes produce nearly identical, positive solutions — confirming that the nonstandard schemes do not distort pricing when CN already works correctly (Theorem 3.1 in the paper). The Monte Carlo reference (5×10⁶ paths, SE < 0.27%) provides an independent check. On the 2000-node log-space grid, FD prices are systematically ~1.5–2% higher than MC, reflecting grid convergence error rather than scheme error (all three FD schemes agree to 14 digits). The paper's S-space results (ΔS=0.05, ~4000 nodes) show similar but not identical values, as expected given the different coordinate system. The barrier implementation uses `FdmDiscreteBarrierStepCondition` with `FdmStepConditionComposite` and a single `FdmBackwardSolver::rollback()` call, matching the validated test-suite pattern.

### Experiment 3: Discrete Double Barrier — Low Volatility (Figure 4)

**Parameters:** K = 100, σ = 0.001, r = 0.05, L = 95, U = 110, T = 1.0, 5 monitoring dates. At σ = 0.001, the `FdmBlackScholesMesher` auto-domain collapses to roughly [99.6, 103.0] in spot-space, excluding both barriers. This experiment uses a `Uniform1dMesher` with explicit bounds [ln(80), ln(130)] to ensure the barriers are represented, matching the validated test pattern.

**Results:** On the uniform log-mesh that includes both barriers, StandardCentral produces negative prices near the barriers — a direct M-matrix violation (37 negative grid nodes). ExponentialFitting and MilevTaglianiCN maintain positivity throughout. The Monte Carlo reference (5×10⁶ paths, SE < 0.005% of price) is plotted alongside the FD curves in Figure 4.

### Experiment 4: Grid Convergence (Figure 5)

**Parameters:** European call, S = 100, K = 100, r = 0.05, q = 0.02, σ = 0.20, T = 1.0. Reference: analytical Black-Scholes formula. Joint space-time refinement with tGrid = 4·xGrid.

**Results:** All three schemes converge at O(h²) rate on the log-log plot. In the moderate-volatility regime (σ² > r), the convergence behavior is virtually identical across schemes, confirming that the nonstandard schemes do not degrade accuracy when the standard scheme already works well.

### Experiment 5: Volatility-Sweep Diagnostics (Figures 6–7)

**Parameters:** Volatility sweep from σ = 0.001 to σ = 0.5 (50 log-spaced values), r = 0.05, q = 0, uniform log-mesh with 200 nodes, mMatrixPolicy = None (raw scheme behavior).

**Results:** Figure 6 plots effective diffusion a_eff vs σ on log-log axes for all three schemes. At low σ (≤ 0.01), the schemes differ by orders of magnitude: MilevTaglianiCN >> ExponentialFitting >> StandardCentral. At high σ (≥ 0.2), all three converge as the base diffusion σ²/2 dominates. This illustrates the papers' core insight: nonstandard schemes add artificial diffusion precisely in the low-volatility regime where the standard scheme produces spurious oscillations.

Figure 7 shows the lower and upper off-diagonal entries of the operator matrix vs σ. The left subplot (lower off-diagonal) uses a symlog scale to accommodate both positive and negative values. StandardCentral's lower off-diagonal becomes negative below σ ≈ 0.02 — an M-matrix violation. ExponentialFitting and MilevTaglianiCN maintain non-negative lower off-diagonals across the entire σ range, confirming their M-matrix compliance.

### Experiment 6: Performance Benchmark (Figure 8)

**Parameters:** European call, same as convergence study. Uses `mMatrixPolicy=None` (no fallback) so scheme identity is unambiguous. The cost metric is the deterministic proxy N_x × N_t (proportional to the number of tridiagonal solves), which avoids the non-reproducibility of wall-clock timing.

**Results:** The cost-accuracy tradeoff is essentially identical across all three schemes. Since σ=0.20 is well into the σ² > r regime, the nonstandard schemes produce the same coefficients as StandardCentral and impose no additional computational cost. The tridiagonal solve dominates, and the scheme coefficient computation (xCothx or r²h²/(8σ²)) is negligible.

### Experiment 7: Péclet Number Dependence (Figure 9)

**Results:** The fitting factor ρ = x·coth(x) smoothly transitions from ρ ≈ 1 (standard central) at low Pe to ρ ≈ |Pe| (upwind) at high Pe. The three-regime numerical implementation (Taylor, direct, asymptotic) provides seamless coverage.

## 5. Discussion

### When to Use Each Scheme

| Scheme | Best For | Limitations |
|--------|----------|-------------|
| **StandardCentral** | Moderate/high volatility, smooth payoffs | Oscillates and produces negative prices when σ² ≪ r or near sharp discontinuities |
| **ExponentialFitting** | All regimes, particularly low volatility | Worst-case first-order accuracy; observed O(h²) convergence on smooth data. Trades guaranteed positivity for potential diffusion |
| **MilevTaglianiCN** | Low volatility with CN time stepping | Requires CN-equivalent time scheme; more diffusive than ExponentialFitting in some regimes |

### Key Tradeoffs

1. **Accuracy vs. Positivity:** StandardCentral is O(h², Δt²) but can produce non-physical negative prices. The nonstandard schemes sacrifice some accuracy (through artificial diffusion) to achieve positivity in typical parameter regimes (though edge cases like negative dividend yield can still cause violations).

2. **Artificial Diffusion:** The ExponentialFitting scheme introduces diffusion driven by the drift μ = r − q − σ²/2 (reducing to approximately (1/2)rSΔS only in the low-vol, q=0 limit). The MT variant introduces diffusion proportional to (1/8)(rΔS/σ)². The latter grows quadratically with r/σ, making MT more diffusive in extreme low-vol regimes.

3. **Fallback Safety:** The `FallbackToExponentialFitting` policy provides a safety net when the MT scheme's M-matrix property is not guaranteed (e.g., with non-CN time schemes).

## 6. Conclusion

The implementation of Duffy's exponential fitting and Milev-Tagliani's CN variant in QuantLib provides robust alternatives to standard Crank-Nicolson for pricing options with discontinuous payoffs or in low-volatility regimes. In typical parameter regimes, both nonstandard schemes eliminate spurious oscillations and produce positive prices (though edge cases such as negative dividend yield can still cause violations). The choice between them depends on the specific problem:

- Use **ExponentialFitting** for general robustness across all parameter regimes
- Use **MilevTaglianiCN** when CN time stepping is already in use and the problem is specifically in the σ² ≪ r regime
- Use **StandardCentral** only when the problem is well-conditioned (σ² > r, smooth payoffs)

## References

1. Milev, M., Tagliani, A. (2010). "Low Volatility Options and Numerical Diffusion of Finite Difference Schemes." *Serdica Mathematical Journal* 36: 223–236.
2. Duffy, D.J. (2004). "A Critique of the Crank Nicolson Scheme Strengths and Weaknesses for Financial Instrument Pricing." *Wilmott Magazine*.
3. Milev, M., Tagliani, A. (2010). "Nonstandard Finite Difference Schemes with Application to Finance: Option Pricing." *Serdica Mathematical Journal* 36: 75–88.
