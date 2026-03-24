# Draft: Academic Paper on Nonstandard Finite Difference Schemes for Black-Scholes Option Pricing

## Paper Objective

Write a detailed, comprehensive, professional academic paper in LaTeX (*.tex) presenting the implementation, analysis, and numerical validation of three spatial discretization schemes for the Black-Scholes PDE: StandardCentral (baseline Crank-Nicolson), ExponentialFitting (Duffy 2004), and MilevTaglianiCNEffectiveDiffusion (Milev-Tagliani 2010). The paper should be suitable for submission to a computational finance or numerical methods journal.

## Target Venue Style

- Professional academic formatting (article class or similar)
- Numbered equations, theorems, definitions
- Proper bibliography with BibTeX
- Publication-quality figures (the 9 existing PDFs from results/figures/)
- Tables with numerical results
- Algorithm pseudocode where appropriate

## Paper Structure

### 1. Introduction
- Motivation: CN scheme produces spurious oscillations when σ² ≪ r (convection-dominated regime)
- Literature context: Duffy (2004) critique of CN, Milev-Tagliani (2010) nonstandard FD schemes
- Contribution: QuantLib implementation of three schemes with comprehensive numerical validation
- Paper organization overview

### 2. Mathematical Framework

#### 2.1 Black-Scholes PDE
- Original S-space formulation
- Log-space transformation x = ln S with constant coefficients for flat-vol models
- Drift μ = r − q − σ²/2

#### 2.2 Spatial Discretization
- General tridiagonal system PV^{n+1} = NV^n
- StandardCentral: centered differences, O(h²) accuracy
  - Coefficients: a_{i,i-1} = σ²/(2h²) − μ/(2h), etc.
- ExponentialFitting: Péclet-dependent fitting factor ρ(Pe) = Pe·coth(Pe)
  - Effective diffusion: a_eff = (σ²/2)·ρ(Pe)
  - Guaranteed M-matrix property
- MilevTaglianiCN: artificial diffusion from 6-point reaction-term stencil
  - Added diffusion: r²h²/(8σ²) with cap
  - Requires CN time scheme (CN-equivalence gate)
  - Deliberate omission of drift correction term (audited, O(h²) bounded)

#### 2.3 Péclet Number and Regime Analysis
- Definition: Pe = μh/σ²
- Three regimes: |Pe| < 10⁻⁶ (Taylor), 10⁻⁶ ≤ |Pe| ≤ 50 (direct), |Pe| > 50 (asymptotic)
- Stable evaluation of x·coth(x) via three-branch implementation
- Regime boundary σ_* where Pe(σ) = 1

#### 2.4 M-Matrix Property
- Definition: negative diagonal, non-negative off-diagonals
- Sufficient condition for monotonicity and positivity
- StandardCentral violation when |Pe| > 1
- ExponentialFitting guarantee proof sketch
- Fallback policy implementation

#### 2.5 Time Discretization
- Crank-Nicolson (θ = 0.5)
- CN-equivalence in 1D: Douglas and CraigSneyd reduce to CN because apply_mixed() = 0
- DampingSteps break CN-equivalence (prepend Implicit Euler steps)
- Theorem 3.1 (CN positivity for σ² > r) and Theorem 3.2 (MT stability: Δt < 1/(rM))

### 3. Implementation in QuantLib

#### 3.1 Architecture Overview
- FdmBlackScholesSpatialDesc configuration struct
- FdmBlackScholesOp operator assembly
- FdmBlackScholesSolver with CN-equivalence gate
- Fallback mechanism (FallbackToExponentialFitting)

#### 3.2 Discrete Barrier Framework
- FdmDiscreteBarrierStepCondition with tolerant time matching
- FdmBlackScholesMesher with concentrated grid points at K, L, U
- Barrier engine: European-only guard, t=0 monitoring
- Truncated call payoff on extended domain

#### 3.3 Numerical Stability
- xCothx three-regime evaluation
- M-matrix diagnostics and fallback
- LazyObject pattern for deferred calculation

### 4. Numerical Experiments

All 9 experiments from generate_data.cpp, with figures and tables:

#### 4.1 Experiment 1: Truncated Call — CN Spurious Oscillations (Fig 1, Fig 2)
- Parameters: r=0.05, q=0, σ=0.001, K=50, U=70, T=5/12
- Grid: 2801×2801 log-space uniform mesh
- Demonstrate: SC produces negative prices near U=70 (M-matrix violation)
- Three-scheme comparison showing EF and MT remain smooth and positive
- Fine-grid ExponentialFitting reference (8× refinement)

#### 4.2 Experiment 2: Moderate-Volatility Discrete Double Barrier (Fig 3)
- Parameters: K=100, σ=0.25, r=0.05, q=0, L=95, U=110, T=0.5, 5 monitoring dates
- Grid: 4000-node concentrated mesh
- All schemes agree (σ² > r satisfies Theorem 3.1)
- Monte Carlo validation: 10M paths, 95% CI error bars
- Paper Table 1 reconstruction at 9 spots: S ∈ {95, 95.0001, 95.5, 99.5, 100, 100.5, 109.5, 109.9999, 110}

#### 4.3 Experiment 3: Low-Volatility Discrete Double Barrier (Fig 4)
- Parameters: K=100, σ=0.001, r=0.05, q=0, L=95, U=110, T=1.0, 5 monitoring dates
- Grid: 800-node uniform mesh (auto-domain excludes barriers at low vol)
- SC produces 37 negative grid nodes
- EF maintains positivity; MT shows visible undershoot near barriers (large artificial diffusion)
- Monte Carlo validation: 5M paths

#### 4.4 Experiment 4: Grid Convergence (Fig 5)
- European call: S=100, K=100, r=0.05, q=0.02, σ=0.20, T=1.0
- 7 refinement levels: N_x ∈ {25, 50, ..., 1600}, N_t = 4·N_x
- Reference: analytical Black-Scholes (V_BS ≈ 9.227)
- All schemes converge O(h²); MT shows ~2.5× coarse-grid error
- Convergence rate table with least-squares fitted α

#### 4.5 Experiment 5: Effective Diffusion σ-Sweep (Fig 6)
- 50 log-spaced σ from 0.001 to 0.5
- Shows a_eff vs σ for all three schemes
- Regime boundary σ_* ≈ 0.0186 where Pe = 1
- Low σ: MT >> EF >> SC; High σ: all converge to σ²/2

#### 4.6 Experiment 6: M-Matrix Off-Diagonal σ-Sweep (Fig 7)
- Same σ range as Experiment 5
- Lower and upper off-diagonal elements vs σ
- SC crosses zero at σ ≈ 0.02
- EF and MT remain non-negative throughout

#### 4.7 Experiment 7: Performance Benchmark (Fig 8)
- European call, 6 grid levels
- Cost (N_x·N_t) vs accuracy and wall-clock time
- SC ≈ 4ms, EF ≈ 7ms at N_x=200; tridiagonal solve dominates
- Cost-accuracy tradeoff analysis

#### 4.8 Experiment 8: xCothx / Péclet Number Regimes (Fig 9)
- Pure function evaluation: ρ(Pe) from Pe = -80 to +80
- Three-regime implementation visualization
- Transition boundaries at |Pe| = 10⁻⁶ and |Pe| = 50

### 5. Results and Discussion

#### 5.1 When to Use Each Scheme
- StandardCentral: σ² > r, smooth payoffs, maximum accuracy
- ExponentialFitting: universal robustness, moderate overhead
- MilevTaglianiCN: low-vol regime with CN time scheme, niche applications

#### 5.2 Key Findings
- M-matrix critical volatility σ_crit ≈ 0.02 for test parameters
- Regime boundary σ_* ≈ 0.0186 (Pe = 1)
- MT drift correction omission bounded by grid error (audited)
- Fallback mechanism eliminates configuration risk

#### 5.3 Comparison with Literature
- Consistency with Milev-Tagliani paper examples
- Extension to discrete barrier products
- QuantLib integration benefits (reuse of existing mesher/solver infrastructure)

### 6. Conclusion
- Summary of contributions
- Practical recommendations for practitioners
- Future work: local vol, multi-asset extensions, American options

### 7. References
- Milev, M. and Tagliani, A. (2010) "Nonstandard Finite Difference Schemes with Application to Finance: Option Pricing." Serdica Mathematical Journal, 36:75-88.
- Milev, M. and Tagliani, A. (2010) "Low Volatility Options and Numerical Diffusion of Finite Difference Schemes." Serdica Mathematical Journal, 36:223-236.
- Duffy, D.J. (2004) "A Critique of the Crank-Nicolson Scheme." Wilmott Magazine.
- Duffy, D.J. (2006) "Finite Difference Methods in Financial Engineering." Wiley.

## Figures to Include

All 9 existing PDF figures from results/figures/:
1. fig1_cn_oscillations.pdf — CN spurious oscillations on truncated call
2. fig2_three_scheme_truncated.pdf — Three-scheme comparison (2-panel: full + zoom)
3. fig3_barrier_moderate.pdf — Moderate-vol barrier with MC validation
4. fig4_barrier_lowvol.pdf — Low-vol barrier with negative price highlighting
5. fig5_convergence.pdf — Grid convergence (log-log)
6. fig6_effective_diffusion.pdf — Effective diffusion vs σ
7. fig7_mmatrix.pdf — M-matrix off-diagonals vs σ (2-panel)
8. fig8_benchmark.pdf — Performance benchmark (2-panel)
9. fig9_xcothx.pdf — xCothx regime visualization

## Tables to Include

1. **Table 1: Paper Table Reconstruction** — Discrete barrier prices at 9 spots (S ∈ {95, 95.0001, ..., 110}) for all 3 schemes + MC reference
2. **Table 2: Grid Convergence** — 7 grid levels × 3 schemes, absolute/relative errors, fitted convergence rates
3. **Table 3: Performance Benchmark** — Grid size, cost, error, wall time for each scheme
4. **Table 4: Scheme Comparison Summary** — Accuracy order, M-matrix guarantee, artificial diffusion, overhead, recommended use case

## LaTeX Requirements

- Single .tex file (or main.tex with possible .bib file)
- Use standard packages: amsmath, amssymb, amsthm, graphicx, booktabs, algorithm2e or algorithmic
- Numbered equations throughout
- Theorem/Definition/Proposition environments for key results
- \includegraphics for the 9 existing PDF figures
- Professional table formatting with booktabs
- BibTeX bibliography
- Proper cross-referencing (\ref, \eqref, \cite)
- Target length: ~20-30 pages

## Data Sources

All numerical data is pre-generated in results/data/ (28 CSV files). The paper should reference and present this data, not regenerate it. The figures are pre-generated in results/figures/ (9 PDF files). The paper should include these via \includegraphics.

## Existing Analysis

The file results/results_analysis.md contains a comprehensive analysis that should inform the paper's discussion sections. The mathematical details are implemented in the QuantLib source code under ql/methods/finitedifferences/.
