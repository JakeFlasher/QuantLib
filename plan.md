# Professional Results Analysis for Nonstandard FD Schemes in QuantLib

## Goal Description

Create a comprehensive, publication-quality results analysis for three spatial discretization schemes (StandardCentral, ExponentialFitting, MilevTaglianiCNEffectiveDiffusion) implemented in the QuantLib 1D Black-Scholes FDM operator. This involves:

1. A **C++ data generator** (`results/generate_data.cpp`) that runs 9 numerical experiments using the existing QuantLib infrastructure and outputs CSV data files.
2. A **Python plotting script** (`results/plot_figures.py`) that produces 9 publication-quality vector PDF figures from the CSV data.
3. A **comprehensive analysis document** (`results/results_analysis.md`) with mathematical framework, experiment descriptions, figure references, and discussion suitable for a top-tier computational finance journal.

The data generator links against the existing QuantLib build at `build/` via CMake `find_package(QuantLib)`. No modifications to existing QuantLib source files (`/ql/*`, `/test-suite/*`) are required.

## Acceptance Criteria

Following TDD philosophy, each criterion includes positive and negative tests for deterministic verification.

- AC-1: C++ data generator compiles, links against the existing QuantLib build, and produces CSV output without any modifications to files in `/ql/` or `/test-suite/`
  - Positive Tests (expected to PASS):
    - `cmake --build results/build` succeeds with exit code 0
    - Running `results/build/generate_data` produces CSV files in `results/data/`
    - All existing r6 tests (`ctest --test-dir build`) still pass unchanged
  - Negative Tests (expected to FAIL):
    - `git diff ql/ test-suite/` shows zero changes (any modification fails this criterion)
    - Running with a missing QuantLib build directory produces a clear error message, not a crash

- AC-2: All CSV output files include full parameter metadata
  - Positive Tests (expected to PASS):
    - Each CSV file contains a header comment block with: scheme name, effective scheme (after any fallback), xGrid, tGrid, r, q, sigma, strike, maturity, mMatrixPolicy, and mesh type
    - CSV header row lists column names matching the data
  - Negative Tests (expected to FAIL):
    - A CSV file missing the scheme name in metadata is rejected by the plotting script with a clear error
    - A CSV file with mismatched column count between header and data rows causes a parse error

- AC-3: Every plotted numerical value is traceable to an analytical reference, a published paper table value, or a fine-grid monotone/MC reference
  - Positive Tests (expected to PASS):
    - European call prices match Black-Scholes analytical formula to within grid convergence tolerance
    - Barrier option prices from fine-grid (4x refinement) reference and MC reference (10^6 paths) agree within 3 standard errors
    - Paper Table 1 values (Milev-Tagliani 2010, Serdica 36:75-88) are reproduced to within the paper's reported precision for matching grid parameters
  - Negative Tests (expected to FAIL):
    - A price curve labeled "analytical" that deviates by more than 1e-4 from the BS formula at any plotted point is flagged
    - An MC reference with fewer than 10^5 paths is rejected as insufficient

- AC-4: All 9 figures are produced in vector PDF format with publication-quality formatting
  - Positive Tests (expected to PASS):
    - `python results/plot_figures.py` produces 9 PDF files in `results/figures/`
    - Figures use LaTeX-compatible fonts (Computer Modern or similar)
    - Color palette is colorblind-friendly (uses distinguishable line styles in addition to color)
    - Axis labels include proper mathematical notation (Greek letters, subscripts)
  - Negative Tests (expected to FAIL):
    - A raster-format output (PNG, JPG) is rejected; only vector PDF is acceptable
    - A figure with overlapping, unreadable legends fails visual inspection

- AC-5: Existing r6 test suite passes unchanged after all work
  - Positive Tests (expected to PASS):
    - `cd build && ctest -R FdmBlackScholes` passes all tests
    - `cd build && ctest -R FdmHyperbolicCot` passes all tests
    - `cd build && ctest -R FdmDiscreteBarrier` passes all tests
  - Negative Tests (expected to FAIL):
    - Any test failure in the above suites blocks completion

- AC-6: A single reproducibility command regenerates all data and figures
  - Positive Tests (expected to PASS):
    - `results/reproduce.sh` runs cmake, builds, executes data generator, and runs plotting script
    - Running `reproduce.sh` on a clean checkout (with existing QuantLib build) produces identical CSV and PDF output
  - Negative Tests (expected to FAIL):
    - A manual multi-step process without the wrapper script fails this criterion

- AC-7: `results_analysis.md` explicitly discusses: (a) log-space formulation, (b) MT artificial diffusion as a diffusion-only adaptation (not the paper's full log-space translation), (c) CN-equivalence requirement for MT, (d) fallback behavior from MT to ExponentialFitting, (e) tradeoffs of each scheme
  - Positive Tests (expected to PASS):
    - Document contains a section titled "Implementation Notes" covering log-space adaptation
    - Document contains a subsection discussing the MT drift correction omission and its O(h^2) bound
    - Document states that MT requires CrankNicolson time stepping
  - Negative Tests (expected to FAIL):
    - A document claiming MT is "the paper scheme verbatim" without qualifying the log-space adaptation is incorrect
    - A document omitting fallback behavior discussion is incomplete

- AC-8: Any MT-to-ExponentialFitting fallback during data generation is recorded in the CSV metadata as the effective scheme used
  - Positive Tests (expected to PASS):
    - When mMatrixPolicy=FallbackToExponentialFitting is used and fallback triggers, the CSV metadata line reads `effective_scheme: ExponentialFitting`
    - When mMatrixPolicy=None is used (no fallback), the effective scheme matches the requested scheme
  - Negative Tests (expected to FAIL):
    - A CSV file that reports `scheme: MilevTaglianiCNEffectiveDiffusion` when fallback actually occurred to ExponentialFitting is incorrect

- AC-9: Monte Carlo reference values for barrier options use at least 10^6 paths and report standard errors
  - Positive Tests (expected to PASS):
    - MC CSV output includes columns for price, standard_error, and num_paths
    - Standard errors are below 0.5% of the price for all reported values
  - Negative Tests (expected to FAIL):
    - MC reference without standard error column is incomplete
    - MC with fewer than 10^5 paths is rejected

## Path Boundaries

### Upper Bound (Maximum Acceptable Scope)
The implementation includes all 9 figures, a complete C++ data generator covering all experiments (truncated call, discrete double barrier at two volatility regimes, grid convergence, operator diagnostics, performance benchmarks, MC references), a Python plotting script with publication-quality output, and a comprehensive `results_analysis.md` document with mathematical framework, all experiment descriptions, figure references, and complete discussion of scheme tradeoffs.

### Lower Bound (Minimum Acceptable Scope)
The implementation includes all 9 figures with at least the C++ data generator producing valid CSV data for each experiment, a Python script generating PDF figures, and a `results_analysis.md` that covers all required discussion points from AC-7. MC references may use 10^5 paths minimum.

### Allowed Choices
- Can use: CMake for the out-of-tree build, matplotlib/pgf for plotting, any colorblind-friendly palette, Uniform1dMesher for operator diagnostics, FdmBlackScholesMesher for barrier experiments
- Can use: QuantLib's existing MC infrastructure or a custom MC loop for barrier references
- Cannot use: modifications to any file in `/ql/` or `/test-suite/`
- Cannot use: QuantLib Python/SWIG bindings (not built in this configuration)
- Must use: CrankNicolson time scheme for MT experiments (CN-equivalence requirement)
- Must use: `TruncatedCallPayoff` pattern (custom payoff on extended domain) for truncated call experiments
- Must use: `FdmDiscreteBarrierStepCondition` with direct rollback for discrete double barrier experiments (NOT single-barrier engine)
- Must use: `FdmBlackScholesMesher` with concentration at K, L, U for barrier experiments
- Must use: Explicit `mMatrixPolicy` setting in every experiment configuration

## Feasibility Hints and Suggestions

> **Note**: This section is for reference and understanding only.

### Conceptual Approach

**Data Generator Architecture:**
```
results/
  CMakeLists.txt          # find_package(QuantLib), add_executable(generate_data ...)
  generate_data.cpp       # Main driver with experiment functions
  reproduce.sh            # Wrapper: cmake + build + run + plot
  plot_figures.py          # matplotlib with LaTeX fonts
  data/                   # CSV output directory (generated)
  figures/                # PDF output directory (generated)
  results_analysis.md     # Analysis document
```

**Key implementation patterns from existing test code:**
- Truncated call payoff: `TruncatedCallPayoff(K, U)` returning `(s >= K && s <= U) ? s - K : 0.0`
- Solver rollback: `FdmBackwardSolver` with `FdmSchemeDesc::CrankNicolson()`
- Barrier step condition: `FdmDiscreteBarrierStepCondition` applied at monitoring dates
- Matrix extraction: `FdmBlackScholesOp::toMatrix()` returns SparseMatrix
- Effective diffusion recovery on uniform mesh: `aUsed = (mat(i,i-1) + mat(i,i+1)) * h^2 / 2`

**Convergence study design:**
Joint space-time convergence with `tGrid = 4 * xGrid` (temporal over-resolution to isolate spatial error). Reference: analytical Black-Scholes formula for European call.

**MC for discrete barriers:**
Simple path-level simulation: generate GBM paths, apply discrete barrier monitoring at each date, compute discounted payoff average. QuantLib's `PseudoRandom` + `MersenneTwisterUniformRng` for random numbers.

### Relevant References
- `test-suite/fdmblackscholespositivity.cpp` — TruncatedCallPayoff class, LowVolSetup fixture, solver patterns
- `test-suite/fdmdiscretebarrierengine.cpp` — Discrete barrier with step conditions
- `test-suite/fdmblackscholesspatialdiscretization.cpp` — Operator coefficient verification, matrix extraction
- `ql/methods/finitedifferences/operators/fdmblackscholesop.cpp` — Operator assembly with scheme dispatch
- `ql/methods/finitedifferences/operators/fdmblackscholesspatialdesc.hpp` — Scheme enum and desc struct
- `ql/methods/finitedifferences/operators/fdmhyperboliccot.hpp` — xCothx implementation
- `ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp` — Solver with spatial desc parameter
- `ql/pricingengines/vanilla/fdblackscholesvanillaengine.cpp` — Engine pattern for FD pricing
- `build/cmake/QuantLibConfig.cmake` — CMake package config for linking

## Dependencies and Sequence

### Milestones
1. **Infrastructure**: Build system and CSV output framework
   - Phase A: Create `results/CMakeLists.txt` with `find_package(QuantLib)`
   - Phase B: Implement CSV writer with metadata headers
   - Phase C: Create `reproduce.sh` wrapper script

2. **Numerical Experiments**: All 9 data generation experiments
   - Phase A: Truncated call experiments (Figs 1-2) using `TruncatedCallPayoff` + `FdmBackwardSolver`
   - Phase B: Discrete double barrier experiments (Figs 3-4) using `FdmDiscreteBarrierStepCondition`
   - Phase C: MC reference generation for barrier experiments (Fig 3-4 validation)
   - Phase D: Grid convergence study (Fig 5) with analytical BS reference
   - Phase E: Operator diagnostics — M-matrix off-diagonals (Fig 7) and effective diffusion (Fig 6)
   - Phase F: Performance benchmark (Fig 8) with median timing
   - Phase G: xCothx/Peclet data generation (Fig 9)

3. **Figure Generation**: Python plotting script
   - Phase A: matplotlib setup with LaTeX fonts, colorblind palette, journal sizing
   - Phase B: Individual figure plotting functions for all 9 figures
   - Phase C: Integration with `reproduce.sh`

4. **Analysis Document**: `results_analysis.md`
   - Phase A: Mathematical framework and scheme descriptions
   - Phase B: Experiment descriptions and figure references
   - Phase C: Results discussion, tradeoff analysis, conclusion

Milestone 2 depends on Milestone 1. Milestone 3 depends on Milestone 2. Milestone 4 depends on Milestone 3 (needs final figure paths). Within Milestone 2, Phases A-G are independent and can be parallelized.

## Task Breakdown

| Task ID | Description | Target AC | Tag | Depends On |
|---------|-------------|-----------|-----|------------|
| task1 | Create `results/CMakeLists.txt` with `find_package(QuantLib)` and CSV output infrastructure | AC-1, AC-2 | coding | - |
| task2 | Implement truncated call experiments (Figs 1-2): `TruncatedCallPayoff` on extended domain, sweep spot values for CN/ExpFit/MT schemes, output CSV | AC-1, AC-3 | coding | task1 |
| task3 | Implement discrete double barrier experiments (Figs 3-4): `FdmDiscreteBarrierStepCondition` with aligned mesh, moderate-vol (σ=0.25) and low-vol (σ=0.001), output CSV | AC-1, AC-3, AC-8 | coding | task1 |
| task4 | Implement MC reference for discrete barriers: GBM path simulation with discrete monitoring, 10^6 paths, report prices + standard errors | AC-3, AC-9 | coding | task1 |
| task5 | Implement grid convergence study (Fig 5): European call with analytical BS reference, joint space-time refinement, output CSV | AC-1, AC-3 | coding | task1 |
| task6 | Implement operator diagnostics (Figs 6-7): extract matrix via `toMatrix()`, recover effective diffusion from off-diagonals, output CSV | AC-1, AC-3 | coding | task1 |
| task7 | Implement performance benchmark (Fig 8): sweep grid sizes, median timing over multiple runs, output CSV | AC-1 | coding | task1 |
| task8 | Implement xCothx/Peclet data (Fig 9): evaluate `detail::xCothx` across Pe range, output CSV | AC-1 | coding | task1 |
| task9 | Create `reproduce.sh` wrapper script | AC-6 | coding | task2, task3, task4, task5, task6, task7, task8 |
| task10 | Create `results/plot_figures.py` with all 9 figures in publication-quality PDF | AC-4 | coding | task2, task3, task4, task5, task6, task7, task8 |
| task11 | Write `results/results_analysis.md` with full mathematical framework, experiment descriptions, figure references, and discussion | AC-7 | coding | task10 |
| task12 | Review correctness of MT formulation description and tradeoff discussion in analysis document | AC-7 | analyze | task11 |
| task13 | Cross-validate barrier results against paper Table 1 and MC references | AC-3, AC-9 | analyze | task3, task4 |
| task14 | Final integration test: run `reproduce.sh` end-to-end, verify all ACs | AC-1 through AC-9 | coding | task9, task10, task11 |

## Claude-Codex Deliberation

### Agreements
- Out-of-tree `results/` directory at repo root is the correct approach, avoiding any `/ql/` or `/test-suite/` modifications
- C++ data generator linking via CMake `find_package(QuantLib)` is feasible and the correct architecture
- Using `FdmBlackScholesOp::toMatrix()` for M-matrix diagnostics and effective diffusion recovery is sound
- Effective diffusion recovery via `aUsed = (lower + upper) * h^2 / 2` works on uniform log-meshes
- Sweeping `FdmBlackScholesSolver::valueAt(S)` for price curves is the right approach
- The shipped MT operator is a log-space diffusion-only adaptation, not the paper's full translation; the analysis must state this clearly
- MT experiments must use `FdmSchemeDesc::CrankNicolson()` per the CN-equivalence requirement
- `mMatrixPolicy` must be explicit in every experiment configuration

### Resolved Disagreements
- **Truncated call implementation**: Claude initially proposed using a boundary at U. Codex correctly identified this changes the product. Resolution: use `TruncatedCallPayoff` on a domain extending past U, matching existing test code pattern.
- **Double barrier engine**: Claude initially proposed using `FdBlackScholesBarrierEngine`. Codex correctly identified this is a single-barrier engine. Resolution: use `FdmDiscreteBarrierStepCondition` with direct rollback, matching `test-suite/fdmdiscretebarrierengine.cpp` pattern.
- **Mesh strategy**: Claude initially proposed uniform log-mesh for everything. Codex noted barriers need mesh alignment. Resolution: use `Uniform1dMesher` for operator diagnostics, `FdmBlackScholesMesher` with concentration at K/L/U for barrier experiments.
- **Convergence study**: Claude initially stated "tGrid = xGrid for spatial convergence isolation." Codex noted this is joint refinement, not true isolation. Resolution: renamed to "joint space-time convergence" with `tGrid = 4 * xGrid` to over-resolve temporally.
- **Reproducibility command**: Claude initially proposed `make && python plot_figures.py`. Codex noted there's no root Makefile. Resolution: use a `reproduce.sh` script wrapping cmake + build + run + plot.

### Convergence Status
- Final Status: `converged` (2 rounds, all REQUIRED_CHANGES addressed, no material DISAGREE remaining)

## Pending User Decisions

- DEC-1: File location for new work
  - Claude Position: `results/` directory at repo root
  - Codex Position: Agreed, but noted constraint boundary
  - Tradeoff Summary: In-repo keeps version control; out-of-repo avoids any repo footprint
  - Decision Status: `results/ in repo root` (user decided)

- DEC-2: Barrier reference methodology
  - Claude Position: Fine-grid monotone reference only
  - Codex Position: Paper MC table and/or fine-grid
  - Tradeoff Summary: Fine-grid is self-contained; MC adds rigor but complexity; paper tables limited to specific spots
  - Decision Status: `Fine-grid + Monte Carlo` (user decided)

- DEC-3: Figure scope
  - Claude Position: 7 core + 2 optional for faster delivery
  - Codex Position: Reduce to 5 strongest for signal
  - Tradeoff Summary: All 9 is comprehensive but slower; 5 is focused but may miss important comparisons
  - Decision Status: `All 9 figures mandatory` (user decided)

## Implementation Notes

### Code Style Requirements
- Implementation code and comments must NOT contain plan-specific terminology such as "AC-", "Milestone", "Step", "Phase", or similar workflow markers
- These terms are for plan documentation only, not for the resulting codebase
- Use descriptive, domain-appropriate naming in code instead

### Experiment Parameter Summary

| Experiment | r | q | σ | K | L | U | T | Mesh | Time Scheme | mMatrixPolicy |
|-----------|---|---|---|---|---|---|---|------|-------------|---------------|
| Truncated call (Figs 1-2) | 0.05 | 0.0 | 0.001 | 50 | - | 70 | 5/12 | Uniform1dMesher, xGrid=2801 | CrankNicolson | None (raw comparison) |
| Barrier moderate-vol (Fig 3) | 0.05 | 0.0 | 0.25 | 100 | 95 | 110 | 0.5 | FdmBlackScholesMesher aligned at K,L,U | CrankNicolson | None |
| Barrier low-vol (Fig 4) | 0.05 | 0.0 | 0.001 | 100 | 95 | 110 | 1.0 | FdmBlackScholesMesher aligned at K,L,U | CrankNicolson | None |
| European convergence (Fig 5) | 0.05 | 0.02 | 0.20 | 100 | - | - | 1.0 | FdmBlackScholesMesher | CrankNicolson | None |
| Operator diagnostics (Figs 6-7) | 0.05 | 0.0 | 0.001 | 100 | - | - | 1.0 | Uniform1dMesher, xGrid=200 | N/A (single time step) | None |
| Benchmark (Fig 8) | 0.05 | 0.02 | 0.20 | 100 | - | - | 1.0 | FdmBlackScholesMesher | CrankNicolson | FallbackToExponentialFitting |
| xCothx (Fig 9) | N/A | N/A | N/A | N/A | - | - | N/A | N/A (pure function) | N/A | N/A |

### Key Mathematical Notes
- The shipped MT operator uses `aUsed = σ²/2 + r²h²/(8σ²)` in log-space, omitting the drift correction `-r²h²/(8σ²)` that the paper's full S-space→log-space translation would include. Test `testMilevTaglianiDriftCorrectionAudit` in `fdmblackscholesspatialdiscretization.cpp` proves this omission is bounded by grid convergence error.
- ExponentialFitting uses `aUsed = (σ²/2) · xCothx(Pe)` where `Pe = μh/σ²` and `μ = r - q - σ²/2`. This is algebraically equivalent to Duffy's `ρ = (μh/2)·coth(μh/(2σ_diff))`.
- For barrier experiments with `mMatrixPolicy=None`, no fallback occurs. This gives raw scheme comparison without the safety net.

--- Original Design Draft Start ---

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

--- Original Design Draft End ---
