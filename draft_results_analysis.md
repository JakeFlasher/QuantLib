# Draft: Professional Results Analysis for Nonstandard FD Schemes in QuantLib

## Context

We have implemented three spatial discretization schemes for the 1D Black-Scholes FDM operator in QuantLib:

1. **StandardCentral** (baseline Crank-Nicolson centered differences)
2. **ExponentialFitting** (Duffy 2004 — fitting factor ρ = (μh/2)·coth(μh/(2σ_diff)))
3. **MilevTaglianiCNEffectiveDiffusion** (Milev-Tagliani 2010 — CN variant with artificial diffusion r²h²/(8σ²))

These are implemented across 22 source files in `/ql/methods/finitedifferences/` and `/ql/pricingengines/`, with 41+ test cases across 5 test suites. The implementation operates in log-space (x = ln(S)).

## Goal

Generate a comprehensive, publication-quality results analysis suitable for submission to a top-tier computational finance journal. This requires:

### Part A: Professional Figure Generation

Create a Python script that links against the QuantLib build (or re-implements the pricing logic using QuantLib's Python bindings or direct C++ executable output) to produce matplotlib/pgf figures for:

1. **Figure 1 — Spurious Oscillations Demonstration**: Reproduce the paper's key result showing CN scheme oscillations near barriers for a truncated call option (r=0.05, σ=0.001, K=50, U=70, T=5/12, Smax=140, ΔS=0.05, Δt=0.01). Plot V(S) for StandardCentral alongside the analytical solution, clearly showing the oscillations near S=70.

2. **Figure 2 — Three-Scheme Solution Comparison (Truncated Call)**: Same parameters as Fig 1 but comparing all three schemes on the same axes. StandardCentral with oscillations, ExponentialFitting smooth, MilevTagliani smooth. Include an inset/zoom near the discontinuity at U=70.

3. **Figure 3 — Discrete Double Barrier Knock-Out Comparison**: Reproduce Example 4.1 from the Milev-Tagliani (2010) Serdica paper. Parameters: K=100, σ=0.25, T=0.5, r=0.05, L=95, U=110, 5 monitoring dates. Compare all schemes against Table 1 values from the paper.

4. **Figure 4 — Low Volatility Double Barrier**: Same as Fig 3 but with σ=0.001, r=0.05 (the σ² ≪ r regime). Show that StandardCentral produces negative prices while ExponentialFitting and MT remain positive.

5. **Figure 5 — Grid Convergence Study**: For a European call option (where analytical BS price is known), plot |V_numerical - V_analytical| vs grid size (xGrid) on log-log scale for all three schemes. Show convergence rates. Use parameters: S=100, K=100, r=0.05, q=0.02, σ=0.20, T=1.0.

6. **Figure 6 — Artificial Diffusion Comparison**: Plot the effective diffusion coefficient a_eff(S) across the mesh for Duffy vs MT at σ=0.001. Show how each scheme inflates the diffusion to maintain M-matrix property.

7. **Figure 7 — M-Matrix Property Visualization**: For each scheme at σ=0.001, plot the off-diagonal entries of the operator matrix vs node index. StandardCentral shows negative entries; Duffy and MT show all non-negative.

8. **Figure 8 — Performance Benchmark**: Runtime vs accuracy trade-off. For a fixed problem (European call), sweep grid sizes and plot (runtime, |error|) for each scheme. Show that non-standard schemes achieve positivity without significant runtime overhead.

9. **Figure 9 — Péclét Number Dependence**: Plot the Duffy fitting factor ρ = xCothx(Pe) vs Pe, and show how it transitions from ρ≈1 (low Pe, standard central) to ρ≈|Pe| (high Pe, upwind).

### Part B: Comprehensive Results Analysis Document

Write a `results_analysis.md` that includes:

1. **Introduction**: Problem statement — why standard CN fails for low volatility / discontinuous payoffs in financial PDE pricing.

2. **Mathematical Framework**: Brief summary of the Black-Scholes PDE, the three discretization approaches, and their theoretical properties (M-matrix, positivity, convergence order).

3. **Implementation Notes**: How the schemes are adapted to log-space in QuantLib. Discussion of the MT drift correction and why the diffusion-only approach is sufficient.

4. **Numerical Experiments**:
   - Experiment 1: Truncated call option (paper Example 4.1 from low_volatility.pdf)
   - Experiment 2: Discrete double barrier knock-out (paper Example 4.1 from 2010-075-088.pdf)
   - Experiment 3: Low-volatility stress test (σ²≪r regime)
   - Experiment 4: Grid convergence against analytical Black-Scholes
   - Experiment 5: Performance benchmarks

5. **Results and Discussion**: Reference all figures, compare with paper tables, discuss implications.

6. **Conclusion**: Summary of contributions and when to use each scheme.

## Technical Approach

The figure generation should use a standalone C++ test executable that outputs CSV data for each experiment, combined with a Python plotting script using matplotlib with publication-quality settings (LaTeX fonts, appropriate sizing for journal columns, colorblind-friendly palettes).

The C++ data-generation program should reuse the existing QuantLib infrastructure (FdmBlackScholesOp, FdmBlackScholesSolver, etc.) to ensure the figures reflect the actual implementation.

## Files to Create/Modify

- `results/generate_data.cpp` — C++ program to run all numerical experiments and output CSV
- `results/plot_figures.py` — Python script to generate all figures from CSV data
- `results/results_analysis.md` — The comprehensive analysis document
- Potentially modify `CMakeLists.txt` to add the data generation target

## In-Scope Papers

- Milev, Tagliani (2010). "Low Volatility Options and Numerical Diffusion of Finite Difference Schemes." Serdica Math J 36: 223-236.
- Duffy (2004). "A Critique of the Crank Nicolson Scheme Strengths and Weaknesses for Financial Instrument Pricing." Wilmott Magazine.
- Milev, Tagliani (2010). "Nonstandard Finite Difference Schemes with Application to Finance: Option Pricing." Serdica Math J 36: 75-88.
