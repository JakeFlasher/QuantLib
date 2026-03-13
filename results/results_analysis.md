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

$$-\frac{\partial V}{\partial t} + rS\frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} - rV = 0$$

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

This guarantees non-negative off-diagonal entries (M-matrix property), producing oscillation-free, positive solutions. Convergence is uniform: |V(S_j, t_n) − U_j^n| ≤ c(h + k), independent of h, k, and σ.

### 2.4 Milev-Tagliani CN Variant

The MT scheme discretizes the reaction term −rV using a 6-point stencil with parameters ω₁ = ω₂ = −r/(16σ²). This introduces artificial diffusion that restores the M-matrix property. In the code's log-space formulation, the effective diffusion is:

$$a_{\text{eff}} = \frac{\sigma^2}{2} + \frac{r^2 h^2}{8\sigma^2}$$

**Important implementation note:** The shipped QuantLib operator implements a log-space *diffusion-only adaptation* of the MT scheme. The paper's full S-space → log-space translation would also include a drift correction term of −r²h²/(8σ²). This correction is O(h²) and has been verified (via the `testMilevTaglianiDriftCorrectionAudit` test) to be bounded by the grid convergence error. See Section 3 for details.

### 2.5 CN-Equivalence Requirement

The MT scheme is designed for Crank-Nicolson time stepping. Using it with other time schemes (e.g., implicit Euler) may trigger the M-matrix fallback to ExponentialFitting when `mMatrixPolicy=FallbackToExponentialFitting` is set.

## 3. Implementation Notes

### 3.1 Log-Space Formulation

All three schemes operate in log-space (x = ln S), where the Black-Scholes PDE has constant coefficients for flat-vol models. The mesh spacing h is in log-space; the physical ΔS varies across the grid as ΔS ≈ S·h.

### 3.2 MT Drift Correction Omission

The paper's full S-space → log-space translation of the MT artificial diffusion yields both:
- Diffusion addition: +r²h²/(8σ²) to the coefficient of ∂²V/∂x²
- Drift correction: −r²h²/(8σ²) to the coefficient of ∂V/∂x

The shipped code applies only the diffusion addition. The drift correction is O(h²) — the same order as the artificial diffusion itself — but testing confirms that the scheme difference (with vs. without drift correction) is dominated by the grid convergence error. On the audited mesh (xGrid=401, r=0.50, σ=0.001), the scheme difference / grid error ratio is < 1.0, and the scheme difference shrinks on finer grids.

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

**Results:** Figure 1 shows severe oscillations in the StandardCentral solution near the upper barrier U = 70, with negative values reaching below −5. The ExponentialFitting and MilevTaglianiCN solutions are smooth and positive throughout.

Figure 2 provides a side-by-side comparison with an inset zoom near U = 70, clearly showing the three schemes' behavior at the discontinuity.

### Experiment 2: Discrete Double Barrier — Moderate Volatility (Figure 3)

**Parameters:** K = 100, σ = 0.25, r = 0.05, L = 95, U = 110, T = 0.5, 5 monitoring dates (paper Example 4.1 from Milev-Tagliani, *Nonstandard FD Schemes*, 2010).

**Results:** At moderate volatility (σ² > r), all three schemes produce similar, positive solutions. The Monte Carlo reference (10⁶ paths) confirms the FD results within statistical error bars. This is the regime where CN works correctly, as shown by Theorem 3.1 in the paper.

### Experiment 3: Discrete Double Barrier — Low Volatility (Figure 4)

**Parameters:** K = 100, σ = 0.001, r = 0.05, L = 95, U = 110, T = 1.0, 5 monitoring dates.

**Results:** In the σ² ≪ r regime, StandardCentral produces negative prices near the barriers — a direct M-matrix violation. ExponentialFitting and MilevTaglianiCN maintain positivity throughout. The MC reference confirms the correct solution shape.

### Experiment 4: Grid Convergence (Figure 5)

**Parameters:** European call, S = 100, K = 100, r = 0.05, q = 0.02, σ = 0.20, T = 1.0. Reference: analytical Black-Scholes formula. Joint space-time refinement with tGrid = 4·xGrid.

**Results:** All three schemes converge at O(h²) rate on the log-log plot. In the moderate-volatility regime (σ² > r), the convergence behavior is virtually identical across schemes, confirming that the nonstandard schemes do not degrade accuracy when the standard scheme already works well.

### Experiment 5: Operator Diagnostics (Figures 6–7)

**Parameters:** σ = 0.001, r = 0.05, uniform log-mesh with 200 nodes.

**Results:** Figure 7 shows the off-diagonal entries of the operator matrix. StandardCentral has negative lower off-diagonals across most of the grid — a clear M-matrix violation. Both ExponentialFitting and MilevTaglianiCN maintain non-negative off-diagonals everywhere.

Figure 6 compares the effective diffusion coefficients. StandardCentral uses σ²/2 = 5×10⁻⁷ (barely visible on log scale). ExponentialFitting inflates this by the fitting factor ρ, while MilevTaglianiCN adds the r²h²/(8σ²) term. Both achieve values orders of magnitude larger than the base diffusion, which is the mechanism that ensures M-matrix compliance.

### Experiment 6: Performance Benchmark (Figure 8)

**Parameters:** European call, same as convergence study. Median of 5 runs per configuration.

**Results:** Runtime scales similarly for all three schemes. The nonstandard schemes impose minimal overhead — the coefficient computation (xCothx or r²h²/(8σ²)) is negligible compared to the tridiagonal solve. The runtime-accuracy tradeoff is nearly identical across schemes.

### Experiment 7: Péclet Number Dependence (Figure 9)

**Results:** The fitting factor ρ = x·coth(x) smoothly transitions from ρ ≈ 1 (standard central) at low Pe to ρ ≈ |Pe| (upwind) at high Pe. The three-regime numerical implementation (Taylor, direct, asymptotic) provides seamless coverage.

## 5. Discussion

### When to Use Each Scheme

| Scheme | Best For | Limitations |
|--------|----------|-------------|
| **StandardCentral** | Moderate/high volatility, smooth payoffs | Oscillates and produces negative prices when σ² ≪ r or near sharp discontinuities |
| **ExponentialFitting** | All regimes, particularly low volatility | First-order accuracy in space (trades accuracy for unconditional stability) |
| **MilevTaglianiCN** | Low volatility with CN time stepping | Requires CN-equivalent time scheme; more diffusive than ExponentialFitting in some regimes |

### Key Tradeoffs

1. **Accuracy vs. Positivity:** StandardCentral is O(h², Δt²) but can produce non-physical negative prices. The nonstandard schemes sacrifice some accuracy (through artificial diffusion) to guarantee positivity.

2. **Artificial Diffusion:** The ExponentialFitting scheme introduces diffusion proportional to (1/2)rSΔS, while the MT variant introduces diffusion proportional to (1/8)(rΔS/σ)². The latter grows quadratically with r/σ, making MT more diffusive in extreme low-vol regimes.

3. **Fallback Safety:** The `FallbackToExponentialFitting` policy provides a safety net when the MT scheme's M-matrix property is not guaranteed (e.g., with non-CN time schemes).

## 6. Conclusion

The implementation of Duffy's exponential fitting and Milev-Tagliani's CN variant in QuantLib provides robust alternatives to standard Crank-Nicolson for pricing options with discontinuous payoffs or in low-volatility regimes. Both nonstandard schemes eliminate spurious oscillations and guarantee positive prices, with minimal runtime overhead. The choice between them depends on the specific problem:

- Use **ExponentialFitting** for general robustness across all parameter regimes
- Use **MilevTaglianiCN** when CN time stepping is already in use and the problem is specifically in the σ² ≪ r regime
- Use **StandardCentral** only when the problem is well-conditioned (σ² > r, smooth payoffs)

## References

1. Milev, M., Tagliani, A. (2010). "Low Volatility Options and Numerical Diffusion of Finite Difference Schemes." *Serdica Mathematical Journal* 36: 223–236.
2. Duffy, D.J. (2004). "A Critique of the Crank Nicolson Scheme Strengths and Weaknesses for Financial Instrument Pricing." *Wilmott Magazine*.
3. Milev, M., Tagliani, A. (2010). "Nonstandard Finite Difference Schemes with Application to Finance: Option Pricing." *Serdica Mathematical Journal* 36: 75–88.
