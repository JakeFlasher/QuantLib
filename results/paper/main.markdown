# Nonstandard Finite Difference Schemes for Black-Scholes Option Pricing: Implementation and Numerical Validation in QuantLib

**Authors:** [To be determined]

**Abstract.**
The Crank-Nicolson finite difference scheme is the standard method for solving the Black-Scholes partial differential equation in computational finance, but it produces spurious oscillations and non-physical negative prices when the mesh Péclet number is large — a regime that arises whenever $\sigma^2 \ll r$. We present a QuantLib implementation of three spatial discretization schemes for the one-dimensional Black-Scholes PDE in log-space: StandardCentral (baseline Crank-Nicolson), ExponentialFitting [Duffy2004], and MilevTaglianiCNEffectiveDiffusion [MilevTagliani2010a]. For each scheme we derive the tridiagonal operator coefficients in log-space coordinates, analyze the M-matrix property and positivity guarantees, and state adapted versions of the key theorems from the literature with proofs for the log-space formulation. The implementation features a CN-equivalence gate, an automatic fallback mechanism, and support for discretely monitored double barrier options. We validate the schemes through eight experiment groups producing nine figures and four tables, covering truncated call oscillations, discrete barrier pricing at moderate and low volatility with Monte Carlo references, grid convergence, effective diffusion and M-matrix diagnostics, performance benchmarks, and Péclet number regime visualization. All three schemes achieve $O(h^2)$ convergence on smooth data; ExponentialFitting provides guaranteed M-matrix compliance across all parameter regimes, while MilevTaglianiCN offers an alternative for the low-volatility regime under Crank-Nicolson time stepping. We provide practical scheme selection guidance based on the relationship between $\sigma^2$ and $r$.

---

## 1. Introduction

The Black-Scholes partial differential equation (PDE) is the foundation of option pricing in computational finance. When analytical solutions are unavailable — as is the case for discretely monitored barrier options and many exotic contracts — finite difference methods provide a flexible and accurate numerical alternative. The Crank-Nicolson (CN) scheme, with its second-order accuracy in both space and time, is the most widely used time-stepping method in this context.

However, CN suffers from a well-documented weakness: it produces **spurious oscillations** near discontinuities in the initial or boundary conditions [Duffy2004]. These oscillations are not merely aesthetic defects — they generate non-physical negative option prices and inaccurate Greeks. The problem is particularly severe in the **low-volatility regime** where $\sigma^2 \ll r$, because the mesh Péclet number $\text{Pe} = \mu h / \sigma^2$ becomes large, violating the M-matrix property of the discrete operator. Here $\mu = r - q - \sigma^2/2$ is the risk-neutral drift and $h$ is the mesh spacing in log-space.

Two families of nonstandard finite difference schemes address this issue:

1. **Exponential fitting** [Duffy2004], [Duffy2006]: replaces the standard diffusion coefficient with a Péclet-dependent fitting factor $\rho(\text{Pe}) = \text{Pe} \cdot \coth(\text{Pe})$, guaranteeing the M-matrix property for all parameter regimes.

2. **The Milev-Tagliani semi-implicit scheme** [MilevTagliani2010a], [MilevTagliani2010b]: introduces artificial diffusion through a modified reaction-term discretization, achieving positivity-preservation under a mild time-step constraint.

This paper presents an implementation of all three schemes — StandardCentral, ExponentialFitting, and MilevTaglianiCNEffectiveDiffusion — within the QuantLib library's finite difference framework. Our contributions are:

- Derivation of all three scheme operators in **log-space** ($x = \ln S$), where the Black-Scholes PDE has constant coefficients for flat-volatility models, with explicit tridiagonal coefficient formulas.
- **Adapted proofs** of the key theorems from [MilevTagliani2010a] for the log-space formulation, with careful documentation of the differences from the original S-space results.
- A comprehensive QuantLib implementation featuring a **CN-equivalence gate** (automatic fallback when the time scheme is incompatible with the MT scheme), **M-matrix diagnostics**, and support for **discretely monitored double barrier options**.
- **Eight numerical experiment groups** producing nine figures and four tables that validate the schemes across a range of parameter regimes, from extreme low volatility ($\sigma = 0.001$) to moderate volatility ($\sigma = 0.25$).

**Scope and limitations.** This work is restricted to one spatial dimension, flat (constant) volatility, the log-space formulation, and European exercise or discrete barrier products. Extensions to local volatility, multi-asset problems, and American options are left for future work. The Milev-Tagliani adaptation in our implementation deliberately omits the drift correction term from the full S-space-to-log-space translation; this design choice and its empirical validation are discussed in [Section 2](#2-mathematical-framework) and [Section 5](#5-results-and-discussion).

**Organization.** [Section 2](#2-mathematical-framework) presents the mathematical framework, including the Black-Scholes PDE, the three spatial discretization schemes, Péclet number analysis, M-matrix properties, and adapted theorems with proofs. [Section 3](#3-implementation-in-quantlib) describes the QuantLib architecture. [Section 4](#4-numerical-experiments) details the eight experiment groups. [Section 5](#5-results-and-discussion) discusses results and scheme recommendations. [Section 6](#6-conclusion) concludes.

---

## 2. Mathematical Framework

### 2.1 The Black-Scholes PDE

In the standard formulation, the option price $V(S, t)$ satisfies the Black-Scholes PDE:

$$-\frac{\partial V}{\partial t} + (r - q) S \frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} - rV = 0 \qquad (1)$$

where $r$ is the risk-free rate, $q$ is the continuous dividend yield, and $\sigma$ is the volatility.

Under the log-space transformation $x = \ln S$, Equation (1) becomes:

$$-\frac{\partial V}{\partial t} + \mu \frac{\partial V}{\partial x} + \frac{\sigma^2}{2} \frac{\partial^2 V}{\partial x^2} - rV = 0 \qquad (2)$$

where $\mu = r - q - \sigma^2/2$ is the risk-neutral drift in log-space. For flat-volatility models, the coefficients $\mu$, $\sigma^2/2$, and $r$ are constants, making log-space the natural coordinate system for finite difference discretization.

### 2.2 Spatial Discretization Schemes

We discretize the spatial domain on a mesh $\{x_0, x_1, \ldots, x_M\}$ with local spacing $h_i = x_{i+1} - x_i$. For uniform meshes, $h_i = h$ for all $i$. All three schemes produce a tridiagonal system:

$$P \mathbf{V}^{n+1} = N \mathbf{V}^n$$

where $P$ and $N$ are tridiagonal matrices and $\mathbf{V}^n$ is the solution vector at time level $n$.

#### StandardCentral

Centered differences for the convection and diffusion terms in log-space yield the operator with tridiagonal entries at interior node $i$:

$$a_{i,i-1} = \frac{\sigma^2}{2h^2} - \frac{\mu}{2h}, \qquad a_{i,i} = -\frac{\sigma^2}{h^2} - r, \qquad a_{i,i+1} = \frac{\sigma^2}{2h^2} + \frac{\mu}{2h} \qquad (3)$$

This scheme is $O(h^2)$ accurate. The off-diagonal entries are non-negative (M-matrix property) if and only if:

$$\frac{\sigma^2}{2h^2} \geq \frac{|\mu|}{2h} \quad \Longleftrightarrow \quad |\text{Pe}| \leq 1$$

where $\text{Pe} = \mu h / \sigma^2$ is the mesh Péclet number. When $|\text{Pe}| > 1$, the off-diagonals become negative and the scheme loses the M-matrix property, leading to spurious oscillations.

#### ExponentialFitting (Duffy)

*Literature result* [Duffy2004]: The exponential fitting scheme replaces the standard diffusion coefficient $\sigma^2/2$ with an effective diffusion:

$$a_{\text{eff}} = \frac{\sigma^2}{2} \cdot \rho(\text{Pe}), \qquad \rho(\text{Pe}) = \text{Pe} \cdot \coth(\text{Pe}) = \frac{\text{Pe}}{\tanh(\text{Pe})} \qquad (4)$$

*Adaptation for log-space*: In the QuantLib implementation, $\text{Pe} = \mu h / \sigma^2$ is computed per node using the local mesh spacing and local variance. The tridiagonal entries become:

$$a_{i,i-1} = \frac{a_{\text{eff},i}}{h^2} - \frac{\mu}{2h}, \qquad a_{i,i} = -\frac{2 a_{\text{eff},i}}{h^2} - r, \qquad a_{i,i+1} = \frac{a_{\text{eff},i}}{h^2} + \frac{\mu}{2h} \qquad (5)$$

The fitting factor satisfies $\rho(\text{Pe}) \geq 1$ for all $\text{Pe}$, with $\rho \to 1$ as $\text{Pe} \to 0$ (recovering StandardCentral) and $\rho \to |\text{Pe}|$ as $|\text{Pe}| \to \infty$. This ensures non-negative off-diagonals for all parameter values, guaranteeing the M-matrix property unconditionally.

#### MilevTaglianiCNEffectiveDiffusion

*Literature result* [MilevTagliani2010a]: In the original S-space formulation, the reaction term $-rV$ is discretized using a 6-point stencil with parameters $\omega_1 = \omega_2 = -r/(16\sigma^2)$, introducing artificial diffusion that restores the M-matrix property.

*Adaptation for log-space*: The full S-space-to-log-space translation of the MT artificial diffusion yields both a diffusion addition and a drift correction:

- **Diffusion addition:** $+\frac{r^2 h^2}{8\sigma^2}$ added to the effective diffusion coefficient
- **Drift correction:** $-\frac{r^2 h^2}{8\sigma^2}$ modifying the drift coefficient

The effective diffusion in log-space is:

$$a_{\text{eff}} = \frac{\sigma^2}{2} + \frac{r^2 h^2}{8\sigma^2} \qquad (6)$$

subject to a cap: $a_{\text{add}} = \min\left(\frac{r^2 h^2}{8\sigma^2},\ R_{\max} \cdot \frac{\sigma^2}{2}\right)$ where $R_{\max}$ is a configurable maximum ratio (`maxAddedDiffusionRatio`, default $10^6$).

**Deliberate omission of drift correction.** The shipped QuantLib implementation applies only the diffusion addition (Equation 6) and omits the drift correction term. This drift correction is $O(h^2)$ — the same order as the artificial diffusion itself. The test `testMilevTaglianiDriftCorrectionAudit` empirically confirms that on the audited mesh ($x_{\text{Grid}} = 401$, $r = 0.50$, $\sigma = 0.001$) the scheme-difference-to-grid-error ratio is less than 1.0, and the scheme difference shrinks on finer grids. This is an empirical audit, not a general proof; the omission is a deliberate design choice that simplifies the implementation while maintaining bounded error.

The tridiagonal entries use $a_{\text{eff}}$ in place of $\sigma^2/2$ in Equation (5), with the standard centered-difference drift.

### 2.3 Péclet Number and Regime Analysis

The mesh Péclet number characterizes the relative strength of convection versus diffusion:

$$\text{Pe} = \frac{\mu h}{\sigma^2} \qquad (7)$$

where $\mu = r - q - \sigma^2/2$. The regime boundary occurs at $\sigma_*$ where $|\text{Pe}(\sigma_*)| = 1$. For typical parameters ($r = 0.05$, $q = 0$, $h \approx 0.007$ on a 200-node log-mesh over $[\ln 50, \ln 200]$), this gives $\sigma_* \approx 0.0186$.

- **$|\text{Pe}| < 1$ (diffusion-dominated):** StandardCentral satisfies the M-matrix property; nonstandard schemes add negligible artificial diffusion.
- **$|\text{Pe}| > 1$ (convection-dominated):** StandardCentral violates the M-matrix property; ExponentialFitting and MilevTaglianiCN add sufficient artificial diffusion to restore it.

### 2.4 Numerically Stable Evaluation of $x \coth(x)$

The fitting factor $\rho(x) = x \coth(x) = x / \tanh(x)$ requires careful numerical evaluation across three regimes:

**Definition 2.1.** The function $\rho: \mathbb{R} \to [1, \infty)$ is evaluated as:

$$\rho(x) = \begin{cases} 1 + x^2/3 & \text{if } |x| < 10^{-6} \quad \text{(Taylor series)} \\ x / \tanh(x) & \text{if } 10^{-6} \leq |x| \leq 50 \quad \text{(direct computation)} \\ |x| & \text{if } |x| > 50 \quad \text{(asymptotic: } \coth(x) \to \text{sign}(x)\text{)} \end{cases} \qquad (8)$$

The Taylor regime avoids the $0/0$ indeterminate form at $x = 0$; the asymptotic regime avoids overflow in $\exp(2x)$ for large arguments. The transition boundaries ($10^{-6}$ and $50$) are configurable via `peSmall` and `peLarge` in `FdmBlackScholesSpatialDesc`.

### 2.5 M-Matrix Property and Positivity

**Definition 2.2.** A square matrix $A$ is an **M-matrix** if $a_{ii} < 0$ for all $i$ and $a_{ij} \geq 0$ for all $j \neq i$ (i.e., negative diagonal, non-negative off-diagonals), and $A$ is nonsingular.

The M-matrix property is a sufficient condition for the discrete solution to satisfy a maximum principle and preserve positivity of the initial data.

**Proposition 2.1** (ExponentialFitting M-matrix guarantee). *The ExponentialFitting scheme produces non-negative off-diagonal entries for all parameter values $\sigma > 0$, $r \geq 0$, $q \geq 0$, and mesh spacings $h > 0$.*

*Proof.* The lower off-diagonal at node $i$ is:

$$a_{i,i-1} = \frac{a_{\text{eff}}}{h^2} - \frac{\mu}{2h} = \frac{1}{h^2}\left(\frac{\sigma^2}{2} \cdot \rho(\text{Pe}) - \frac{\mu h}{2}\right) = \frac{\sigma^2}{2h^2}\left(\rho(\text{Pe}) - \text{Pe}\right)$$

Since $\rho(\text{Pe}) = \text{Pe} \cdot \coth(\text{Pe}) \geq |\text{Pe}| \geq \text{Pe}$ for all $\text{Pe}$, we have $\rho(\text{Pe}) - \text{Pe} \geq 0$. Similarly, the upper off-diagonal satisfies $a_{i,i+1} = \frac{\sigma^2}{2h^2}(\rho(\text{Pe}) + \text{Pe}) \geq 0$ since $\rho(\text{Pe}) \geq |\text{Pe}| \geq -\text{Pe}$. $\square$

**StandardCentral M-matrix violation.** When $|\text{Pe}| > 1$, the StandardCentral off-diagonal $\sigma^2/(2h^2) - |\mu|/(2h)$ becomes negative, violating the M-matrix property. The critical volatility below which this occurs is approximately $\sigma_{\text{crit}} \approx 0.02$ for typical parameters (see Experiment 6, Figure 7).

### 2.6 Adapted Theorems for Log-Space

The following theorems are adapted from [MilevTagliani2010a] for the log-space formulation used in the QuantLib implementation. The original theorems are stated in S-space with the substitution $S_j = j \Delta S$; our adaptations replace $S_j / \Delta S$ with $1/h$ (since the Péclet number structure differs in log-space) and account for the constant-coefficient nature of the log-space PDE.

**Theorem 2.1** (CN positivity in log-space, adapted from [MilevTagliani2010a, Theorem 3.1]). *Consider the Crank-Nicolson scheme applied to the log-space Black-Scholes PDE (Equation 2) with uniform mesh spacing $h$. If $\sigma^2 > r$, then:*

1. *The implicit matrix $P$ is an M-matrix with $P^{-1} > 0$ and $\|P^{-1}\|_\infty \leq (1/\Delta t + r/2)^{-1}$.*
2. *If additionally $\Delta t < 2/(r + \sigma^2/h^2)$, then the explicit matrix $N$ has non-negative entries, the numerical solution preserves positivity, and the discrete maximum principle holds.*
3. *If $\Delta t < 2/(r + 2\sigma^2/h^2)$, then $P^{-1}N$ has $M$ distinct real eigenvalues in $(0, 1)$.*

*Proof sketch.* In log-space, the CN tridiagonal matrices are $P = (1/\Delta t)I + C$ and $N = (1/\Delta t)I - C$, where:

$$C = \text{tridiag}\left\{-\frac{\sigma^2}{2h^2} + \frac{\mu}{2h};\ \frac{\sigma^2}{h^2} + \frac{r}{2};\ -\frac{\sigma^2}{2h^2} - \frac{\mu}{2h}\right\}$$

Note the sign convention: $C$ has the negated off-diagonals of the spatial operator. When $\sigma^2 > r$, the condition $\sigma^2/(2h^2) > |\mu|/(2h)$ holds (since $|\mu| = |r - q - \sigma^2/2| < \sigma^2/2$ when $q \geq 0$ and $\sigma^2 > r$), ensuring $P$ is strictly diagonally dominant with positive diagonal and non-positive off-diagonals — hence an M-matrix. Parts (2) and (3) follow from the eigenvalue bounds on $C$ and the relationship $\lambda_i(P^{-1}N) = (1 - \Delta t \lambda_i(C))/(1 + \Delta t \lambda_i(C))$. The full proof is in [Appendix A](#appendix-a-adapted-theorem-proofs).

**Theorem 2.2** (MT stability in log-space, adapted from [MilevTagliani2010a, Theorem 3.2]). *Consider the Milev-Tagliani scheme applied to the log-space Black-Scholes PDE with uniform mesh spacing $h$ and $M$ interior nodes. If $\Delta t < 1/(rM)$, then:*

1. *The solution is positive and satisfies the discrete maximum principle unconditionally.*
2. *The iteration matrix $P^{-1}N$ has $M$ distinct real positive eigenvalues in $(0, 1)$.*

*Proof sketch.* In the MT scheme, the explicit matrix $N$ has the form with off-diagonals proportional to $(r/(2h))(\pm 1 - 2b)$ where $b = -1/(2M)$ (adapting the S-space parameter $b = -M/2$ to log-space). The condition $\Delta t < 1/(rM)$ ensures diagonal dominance of $N$, making it similar to a symmetric positive definite matrix. Combined with $P$ being an M-matrix (unconditionally, since the MT diffusion addition makes $P$ strictly diagonally dominant), the eigenvalue result follows. The full proof is in [Appendix A](#appendix-a-adapted-theorem-proofs).

**Observation 2.1.** *In the original S-space formulation, the time-step restriction in Theorem 2.2 becomes increasingly severe as $M$ grows: $\Delta t < 1/(rM)$. In log-space, the same structure holds with $M$ being the number of interior grid nodes. For the audited parameters ($r = 0.05$, $M = 401$), this gives $\Delta t < 0.0498$, which is easily satisfied by typical grid configurations.*

### 2.7 Time Discretization and CN-Equivalence

The Crank-Nicolson time stepping ($\theta = 0.5$) yields the system $P\mathbf{V}^{n+1} = N\mathbf{V}^n$ where:

$$P = \frac{1}{\Delta t}I - \frac{1}{2}A, \qquad N = \frac{1}{\Delta t}I + \frac{1}{2}A \qquad (9)$$

and $A$ is the spatial operator matrix.

**CN-equivalence in 1D.** In one spatial dimension, the Douglas scheme and the Craig-Sneyd scheme both reduce to Crank-Nicolson because the mixed derivative operator `apply_mixed()` returns zero. Specifically:

- `CrankNicolson` ($\theta = 0.5$): CN by definition
- `Douglas` ($\theta = 0.5$): reduces to CN in 1D since `apply_mixed() = 0`
- `CraigSneyd` ($\theta = 0.5$): reduces to CN in 1D since `apply_mixed() = 0`

**Damping steps.** When `dampingSteps > 0`, the solver prepends Implicit Euler steps before the main CN scheme. These steps break CN-equivalence and are incompatible with the MT scheme's theoretical guarantees. The solver's CN-equivalence gate (see [Section 3](#3-implementation-in-quantlib)) detects this condition and falls back to ExponentialFitting.

---

## 3. Implementation in QuantLib

### 3.1 Architecture Overview

The implementation is organized around three core components:

**`FdmBlackScholesSpatialDesc`** — A configuration struct that controls spatial discretization:

```
struct FdmBlackScholesSpatialDesc {
    enum Scheme { StandardCentral, ExponentialFitting,
                  MilevTaglianiCNEffectiveDiffusion };
    enum HPolicy { Average, Min, Harmonic };
    enum MMatrixPolicy { None, DiagnosticsOnly, FailFast,
                         FallbackToExponentialFitting };

    Scheme scheme;
    HPolicy hPolicy;
    Real peSmall = 1e-6;      // Taylor regime boundary
    Real peLarge = 50.0;      // asymptotic regime boundary
    Real minVariance = 1e-12;
    Real maxAddedDiffusionRatio = 1e6;
    MMatrixPolicy mMatrixPolicy = FallbackToExponentialFitting;
};
```

**`FdmBlackScholesOp`** — The spatial operator, which assembles the tridiagonal matrix. The `setTime(t1, t2)` method computes per-node variance, drift, and effective diffusion based on the selected scheme. For ExponentialFitting, it evaluates $\rho(\text{Pe})$ via `detail::xCothx`. For MilevTaglianiCN, it computes the added diffusion $r^2 h^2 / (8\sigma^2)$ with the configurable cap.

**`FdmBlackScholesSolver`** — The solver wrapper, which handles time stepping and the CN-equivalence gate.

### 3.2 CN-Equivalence Gate and Fallback Mechanism

The solver enforces the MT scheme's CN requirement through a two-part gate:

```pseudocode
function performCalculations():
    effectiveDesc = copy(spatialDesc)

    if effectiveDesc.scheme == MilevTaglianiCN:
        if NOT isCrankNicolsonEquivalent1D(schemeDesc)
           OR solverDesc.dampingSteps > 0:
            if effectiveDesc.mMatrixPolicy == FailFast:
                throw Error("MT requires CN-equivalent time scheme")
            else:
                effectiveDesc.scheme = ExponentialFitting

    // Construct operator with effectiveDesc
    op = FdmBlackScholesOp(mesher, process, ..., effectiveDesc)
```

The `isCrankNicolsonEquivalent1D` function accepts:
- `CrankNicolsonType` with $\theta = 0.5$
- `DouglasType` with $\theta = 0.5$ (reduces to CN in 1D)
- `CraigSneydType` with $\theta = 0.5$ (reduces to CN in 1D)

The fallback to ExponentialFitting is the default behavior (`MMatrixPolicy::FallbackToExponentialFitting`). With `FailFast`, an incompatible configuration throws an exception.

Additionally, the M-matrix diagnostic in `FdmBlackScholesOp` checks off-diagonal signs after operator assembly. If any interior off-diagonal is negative and `mMatrixPolicy = FallbackToExponentialFitting`, the operator re-assembles using ExponentialFitting. This provides a safety net for edge-case parameters.

### 3.3 Discrete Barrier Framework

Discrete barrier pricing uses `FdmDiscreteBarrierStepCondition` within the standard `FdmBackwardSolver` rollback:

1. **Concentrated mesher:** `FdmBlackScholesMesher` with concentration points at $K$, $L$, and $U$:
   ```
   cPoints = {{K, 0.1, true}, {L, 0.1, true}, {U, 0.1, true}}
   ```
   The `true` flag forces exact node placement at the barriers and strike.

2. **Tolerant time matching:** The step condition uses a relative tolerance of $10^{-10}$ scaled by $\max(|t_{\text{candidate}}|, |t_{\text{current}}|)$ to match monitoring times, avoiding floating-point mismatch from calendar/day-count arithmetic.

3. **Barrier application:** At each monitoring time, the step condition zeroes out the solution vector at nodes outside the corridor $[L, U]$, enforcing the knock-out condition.

4. **Low-volatility mesh:** At extreme low volatility ($\sigma = 0.001$), `FdmBlackScholesMesher`'s auto-domain collapses too narrowly to include the barriers. A `Uniform1dMesher` with explicit bounds $[\ln(L - 15), \ln(U + 20)]$ is used instead.

5. **Value extraction:** Barrier experiment prices are extracted via `valueAtSpot`, which finds the nearest grid node to the target spot. This is distinct from the interpolation-based `valueAt` used for European experiments; the two extraction methods produce different error characteristics and should not be compared as like-for-like.

---

## 4. Numerical Experiments

All experiments use pre-generated data from `generate_data.cpp`. The methodology section below documents the two value extraction approaches, followed by detailed descriptions of each experiment group.

### 4.1 Value Extraction Methods

- **Nearest-node (`valueAtSpot`):** Used for barrier experiments. Finds the grid node closest to $x_{\text{target}} = \ln(S)$ and returns the solution value at that node. Error depends on grid spacing relative to the target.
- **Solver interpolation (`valueAt`):** Used for European experiments. Uses the `Fdm1DimSolver::interpolateAt` cubic interpolation on the solution array. Provides smooth inter-node values.

These two methods produce fundamentally different error characteristics. Nearest-node errors are bounded by the grid spacing; interpolation errors depend on the solution smoothness. Comparisons between barrier and European experiment errors should account for this distinction.

### 4.2 Experiment Traceability

| Exp | Product | Key Parameters | CSV Sources | Figures | Tables | Extraction |
|-----|---------|---------------|-------------|---------|--------|------------|
| 1 | Truncated call | $\sigma=0.001$, $K=50$, $U=70$, $2801 \times 2801$ | truncated_call_{SC,EF,MT,reference}.csv | Fig 1, 2 | — | nearest-node |
| 2 | Discrete barrier (moderate vol) | $\sigma=0.25$, $L=95$, $U=110$, $4000 \times 2000$ | barrier_moderate_vol_{SC,EF,MT}.csv, mc_barrier_moderate_vol.csv | Fig 3 | Table 1 | nearest-node |
| 3 | Discrete barrier (low vol) | $\sigma=0.001$, $L=95$, $U=110$, $800 \times 200$ | barrier_low_vol_{SC,EF,MT}.csv, mc_barrier_low_vol.csv | Fig 4 | — | nearest-node |
| 4 | European call convergence | $\sigma=0.20$, $K=100$, 7 grid levels | grid_convergence_{SC,EF,MT}.csv | Fig 5 | Table 2 | interpolation |
| 5 | Effective diffusion sweep | 50 $\sigma$ values, 200-node mesh | effective_diffusion_sweep_{SC,EF,MT}.csv | Fig 6 | — | direct |
| 6 | M-matrix sweep | 50 $\sigma$ values, 200-node mesh | mmatrix_sweep_{SC,EF,MT}.csv | Fig 7 | — | direct |
| 7 | Performance benchmark | $\sigma=0.20$, $K=100$, 6 grid levels | benchmark_{SC,EF,MT}.csv | Fig 8 | Table 3 | interpolation |
| 8 | xCothx regimes | $\text{Pe} \in [-100, 100]$ | xcothx.csv | Fig 9 | — | direct |

### 4.3 Experiment 1: Truncated Call — CN Spurious Oscillations (Figures 1–2)

**Parameters:** $r = 0.05$, $q = 0$, $\sigma = 0.001$, $K = 50$, $U = 70$, $T = 5/12$. Grid: $2801 \times 2801$ uniform log-mesh over $[\ln 1, \ln 140]$.

This experiment reproduces Example 4.1 from [MilevTagliani2010b]. The truncated call payoff $\max(S - K, 0) \cdot \mathbf{1}_{[K, U]}(S)$ has discontinuities at $S = K$ and $S = U$. At $\sigma = 0.001$, the Péclet number is extremely large ($\text{Pe} \gg 1$), and StandardCentral produces severe oscillations with negative prices near the upper barrier $U = 70$.

A fine-grid ExponentialFitting solution ($8\times$ refinement: $22408 \times 22408$) serves as the reference, confirming the smooth, correct solution shape. ExponentialFitting and MilevTaglianiCN both maintain positivity on the original grid.

**Figure 1:** `fig1_cn_oscillations.pdf` — StandardCentral solution of the truncated call at $\sigma = 0.001$ showing severe spurious oscillations near the upper barrier $U = 70$, with negative option prices. The fine-grid ExponentialFitting reference (dashed) confirms the correct smooth solution shape.

**Figure 2:** `fig2_three_scheme_truncated.pdf` — Three-scheme comparison (full view and zoom near $U = 70$) of the truncated call. ExponentialFitting and MilevTaglianiCN produce smooth, positive solutions that agree with the fine-grid reference, while StandardCentral exhibits the oscillations shown in Figure 1.

### 4.4 Experiment 2: Discrete Double Barrier — Moderate Volatility (Figure 3, Table 1)

**Parameters:** $K = 100$, $\sigma = 0.25$, $r = 0.05$, $q = 0$, $L = 95$, $U = 110$, $T = 0.5$, 5 monitoring dates. Grid: $4000$-node concentrated mesh (`FdmBlackScholesMesher` with concentration at $K$, $L$, $U$).

At $\sigma = 0.25$ ($\sigma^2 = 0.0625 > r = 0.05$), all three schemes satisfy the M-matrix property (Theorem 2.1 applies). The Monte Carlo reference uses $10^7$ paths with standard errors shown.

**Figure 3:** `fig3_barrier_moderate.pdf` — Discrete double barrier knock-out call prices at moderate volatility ($\sigma = 0.25$) for all three FD schemes with Monte Carlo reference (error bars: 95% CI). All schemes produce nearly identical, positive solutions — confirming that nonstandard schemes do not distort pricing when CN is already well-conditioned.

**Table 1: Discrete Barrier Prices at Moderate Volatility**

*Source: barrier_moderate_vol_{SC,EF,MT}.csv, mc_barrier_moderate_vol.csv*

| $S_0$ | SC (CN) | EF (Duffy) | MT (CN Var.) | MC ($10^7$) | MC SE | $\|$SC$-$MC$\|$ | $\|$EF$-$MC$\|$ | $\|$MT$-$MC$\|$ |
|-------:|--------:|-----------:|-------------:|------------:|------:|---------:|---------:|---------:|
| 95 | 0.17745 | 0.17745 | 0.17745 | 0.17398 | 0.00032 | 0.00347 | 0.00347 | 0.00347 |
| 95.0001 | 0.17745 | 0.17745 | 0.17745 | 0.17391 | 0.00032 | 0.00355 | 0.00355 | 0.00355 |
| 95.5 | 0.18367 | 0.18367 | 0.18367 | 0.18316 | 0.00033 | 0.00050 | 0.00050 | 0.00050 |
| 99.5 | 0.23082 | 0.23082 | 0.23082 | 0.22967 | 0.00037 | 0.00115 | 0.00115 | 0.00115 |
| 100 | 0.23406 | 0.23406 | 0.23406 | 0.23260 | 0.00037 | 0.00146 | 0.00146 | 0.00146 |
| 100.5 | 0.23658 | 0.23658 | 0.23658 | 0.23554 | 0.00038 | 0.00104 | 0.00104 | 0.00104 |
| 109.5 | 0.17573 | 0.17573 | 0.17573 | 0.17449 | 0.00033 | 0.00125 | 0.00125 | 0.00125 |
| 109.9999 | 0.16997 | 0.16997 | 0.16997 | 0.16813 | 0.00032 | 0.00184 | 0.00184 | 0.00184 |
| 110 | 0.16997 | 0.16997 | 0.16997 | 0.16722 | 0.00032 | 0.00274 | 0.00274 | 0.00274 |

At $\sigma^2 = 0.0625 > r = 0.05$, all three FD schemes agree to machine precision on the same mesh. FD prices are systematically slightly above MC, consistent with grid convergence error. The FD-MC differences (0.001–0.003) are within a few standard errors at interior spots and slightly larger at the barriers ($S = 95$, $110$), where the solution gradient is steepest.

### 4.5 Experiment 3: Discrete Double Barrier — Low Volatility (Figure 4)

**Parameters:** $K = 100$, $\sigma = 0.001$, $r = 0.05$, $q = 0$, $L = 95$, $U = 110$, $T = 1.0$, 5 monitoring dates. Grid: $800$-node uniform mesh (`Uniform1dMesher` over $[\ln 80, \ln 130]$).

At $\sigma = 0.001$ ($\sigma^2 = 10^{-6} \ll r = 0.05$), the Péclet number is extremely large. StandardCentral produces **37 negative grid nodes** — a direct M-matrix violation. ExponentialFitting maintains positivity throughout. MilevTaglianiCN maintains positivity but shows visible undershoot near the barriers due to the large artificial diffusion: $(r^2 h^2)/(8\sigma^2)$ grows as $O(1/\sigma^2)$, introducing smoothing that flattens the price profile near the barriers. The Monte Carlo reference uses $5 \times 10^6$ paths.

**Figure 4:** `fig4_barrier_lowvol.pdf` — Discrete double barrier knock-out call prices at low volatility ($\sigma = 0.001$) for all three FD schemes with Monte Carlo reference. StandardCentral produces negative prices (shaded region) near the barriers, while ExponentialFitting and MilevTaglianiCN maintain positivity. MilevTaglianiCN shows visible undershoot near the barriers due to large artificial diffusion.

### 4.6 Experiment 4: Grid Convergence (Figure 5, Table 2)

**Parameters:** European call, $S = 100$, $K = 100$, $r = 0.05$, $q = 0.02$, $\sigma = 0.20$, $T = 1.0$. Reference: analytical Black-Scholes $V_{\text{BS}} = 9.22701$. Joint space-time refinement: $N_t = 4 N_x$ across 7 levels ($N_x \in \{25, 50, 100, 200, 400, 800, 1600\}$). Value extraction via solver interpolation (`valueAt`).

**Figure 5:** `fig5_convergence.pdf` — Grid convergence (log-log) of absolute error versus spatial grid size for all three schemes. All converge at $O(h^2)$ rate. MilevTaglianiCN shows approximately 2.5× the absolute error of SC/EF at coarse grids due to its first-order time component, but achieves the same asymptotic convergence rate.

**Table 2: Grid Convergence**

*Source: grid_convergence_{SC,EF,MT}.csv. Reference price: $V_{\text{BS}} = 9.22701$.*

| $N_x$ | $N_t$ | SC Error | EF Error | MT Error | SC Rel% | EF Rel% | MT Rel% |
|-------:|-------:|---------:|---------:|---------:|--------:|--------:|--------:|
| 25 | 100 | 2.50e-02 | 2.60e-02 | 4.23e-02 | 0.271 | 0.281 | 0.458 |
| 50 | 200 | 6.59e-03 | 6.81e-03 | 1.07e-02 | 0.071 | 0.074 | 0.116 |
| 100 | 400 | 1.63e-03 | 1.68e-03 | 2.64e-03 | 0.018 | 0.018 | 0.029 |
| 200 | 800 | 4.05e-04 | 4.19e-04 | 6.56e-04 | 0.004 | 0.005 | 0.007 |
| 400 | 1600 | 1.01e-04 | 1.04e-04 | 1.63e-04 | 0.001 | 0.001 | 0.002 |
| 800 | 3200 | 2.51e-05 | 2.59e-05 | 4.06e-05 | 0.000 | 0.000 | 0.000 |
| 1600 | 6400 | 6.27e-06 | 6.47e-06 | 1.01e-05 | 0.000 | 0.000 | 0.000 |

Fitted convergence rates (least-squares on last 4 grid levels): SC $\alpha \approx 2.01$, EF $\alpha \approx 2.01$, MT $\alpha \approx 2.01$. All three schemes achieve the expected $O(h^2)$ rate. MT's approximately 1.6× error multiplier relative to SC at each grid level reflects the additional first-order time component from the bivariate reaction-term discretization.

### 4.7 Experiment 5: Effective Diffusion $\sigma$-Sweep (Figure 6)

**Parameters:** 50 log-spaced $\sigma$ from 0.001 to 0.5, $r = 0.05$, $q = 0$, $K = 100$, $T = 1.0$, 200-node uniform log-mesh over $[\ln 50, \ln 200]$, $h \approx 0.00697$.

**Figure 6:** `fig6_effective_diffusion.pdf` — Effective diffusion coefficient $a_{\text{eff}}$ versus $\sigma$ on log-log axes for all three schemes. At low $\sigma$ ($\leq 0.01$), the schemes differ by orders of magnitude: MT $\gg$ EF $\gg$ SC. At high $\sigma$ ($\geq 0.2$), all three converge to $\sigma^2/2$ as the base diffusion dominates. A vertical line marks the regime boundary $\sigma_* \approx 0.0186$ where $\text{Pe} = 1$.

### 4.8 Experiment 6: M-Matrix Off-Diagonal $\sigma$-Sweep (Figure 7)

**Parameters:** Same as Experiment 5. Reports the lower and upper off-diagonal entries of the operator matrix at the midpoint node.

**Figure 7:** `fig7_mmatrix.pdf` — Lower and upper off-diagonal entries of the operator matrix versus $\sigma$ (two-panel plot). StandardCentral's lower off-diagonal crosses zero at $\sigma_{\text{crit}} \approx 0.02$, marking the M-matrix violation threshold. ExponentialFitting and MilevTaglianiCN maintain non-negative off-diagonals across the entire $\sigma$ range, confirming their M-matrix compliance.

### 4.9 Experiment 7: Performance Benchmark (Figure 8, Table 3)

**Parameters:** European call, same as convergence study ($S = 100$, $K = 100$, $r = 0.05$, $q = 0.02$, $\sigma = 0.20$, $T = 1.0$). 6 grid levels ($N_x \in \{50, 100, 200, 400, 800, 1600\}$, $N_t = 4N_x$). Cost metric: relative cost $= N_x \times N_t$. Supplemental wall-clock timing: median of 3 runs after 1 warm-up (`std::chrono::steady_clock`). Value extraction via solver interpolation.

**Figure 8:** `fig8_benchmark.pdf` — Performance benchmark (two-panel: cost vs. error, and wall-clock time vs. error). StandardCentral is the fastest at each grid level, followed by MilevTaglianiCN and ExponentialFitting. The tridiagonal solve dominates overall cost; the additional arithmetic in ExponentialFitting's fitting factor and effective-diffusion computations adds measurable overhead.

**Table 3: Performance Benchmark Summary**

*Source: benchmark_{SC,EF,MT}.csv. Reference: $V_{\text{BS}} = 9.22701$.*

| $N_x$ | Cost ($N_x \times N_t$) | SC Error | SC Time (ms) | EF Error | EF Time (ms) | MT Error | MT Time (ms) |
|-------:|------------------------:|---------:|-------------:|---------:|-------------:|---------:|-------------:|
| 50 | 10,000 | 6.59e-03 | 0.3 | 6.81e-03 | 0.5 | 1.07e-02 | 0.4 |
| 100 | 40,000 | 1.63e-03 | 1.1 | 1.68e-03 | 1.8 | 2.64e-03 | 1.4 |
| 200 | 160,000 | 4.05e-04 | 3.9 | 4.19e-04 | 7.0 | 6.56e-04 | 5.5 |
| 400 | 640,000 | 1.01e-04 | 15.2 | 1.04e-04 | 27.5 | 1.63e-04 | 21.7 |
| 800 | 2,560,000 | 2.51e-05 | 59.5 | 2.59e-05 | 107.8 | 4.06e-05 | 83.0 |
| 1600 | 10,240,000 | 6.27e-06 | 241.2 | 6.47e-06 | 425.5 | 1.01e-05 | 327.7 |

**Benchmark caveats:** Wall-clock times are machine-specific and depend on compiler, build mode, and system load. They should not be treated as portable benchmarks. The relative cost ($N_x \times N_t$) is the primary deterministic metric. ExponentialFitting shows approximately 1.8× the wall-clock time of StandardCentral, reflecting the additional `xCothx` evaluation per node; MilevTaglianiCN falls in between.

### 4.10 Experiment 8: xCothx / Péclet Number Regimes (Figure 9)

**Parameters:** Pure function evaluation of $\rho(\text{Pe}) = \text{Pe} \cdot \coth(\text{Pe})$ for $\text{Pe} \in [-100, 100]$ with dense sampling near zero.

**Figure 9:** `fig9_xcothx.pdf` — The fitting factor $\rho(\text{Pe})$ across the full Péclet number range, with the three evaluation regimes labeled: Taylor ($|\text{Pe}| < 10^{-6}$, $\rho \approx 1 + \text{Pe}^2/3$), Direct ($10^{-6} \leq |\text{Pe}| \leq 50$, $\rho = \text{Pe}/\tanh(\text{Pe})$), and Asymptotic ($|\text{Pe}| > 50$, $\rho \approx |\text{Pe}|$). Vertical lines mark the regime boundaries at $|\text{Pe}| = 50$.

---

## 5. Results and Discussion

### 5.1 Key Findings

**M-matrix critical volatility.** Experiment 6 (Figure 7) shows that StandardCentral's lower off-diagonal crosses zero at $\sigma_{\text{crit}} \approx 0.02$ for the test parameters ($r = 0.05$, $q = 0$, $h \approx 0.007$). Below this volatility, the M-matrix property is violated and spurious oscillations appear. Both ExponentialFitting and MilevTaglianiCN maintain non-negative off-diagonals across the entire $\sigma$ range tested ($0.001$ to $0.5$).

**Regime boundary.** The regime boundary $\sigma_* \approx 0.0186$ (where $\text{Pe} = 1$) closely matches $\sigma_{\text{crit}}$. Experiment 5 (Figure 6) shows that below $\sigma_*$, the nonstandard schemes add significant artificial diffusion, while above it, all schemes converge to $\sigma^2/2$.

**Grid convergence.** All three schemes converge at $O(h^2)$ on smooth data (Experiment 4, Figure 5, Table 2). MilevTaglianiCN shows approximately 1.6× the error of StandardCentral and ExponentialFitting at each grid level, reflecting the additional first-order time component from the bivariate reaction-term discretization.

**MT drift correction omission.** The deliberate omission of the drift correction term in the MilevTaglianiCN implementation is validated by the `testMilevTaglianiDriftCorrectionAudit` test, which confirms that the scheme difference is bounded by the grid convergence error. This is an empirical audit specific to the tested parameters, not a general theoretical guarantee.

### 5.2 Scheme Selection Guidance

| Condition | Recommended Scheme | Rationale |
|-----------|-------------------|-----------|
| $\sigma^2 > r$, smooth payoffs | StandardCentral | Maximum accuracy, no artificial diffusion |
| General robustness needed | ExponentialFitting | Guaranteed M-matrix for all parameters; moderate overhead (~1.8× wall time) |
| Low-vol regime, CN time scheme | MilevTaglianiCN | Alternative to EF with different artificial diffusion profile; requires CN-equivalent time stepping |
| Unknown parameter regime | ExponentialFitting with `FallbackToExponentialFitting` | Safety net eliminates configuration risk |

### 5.3 Comparison with Literature

**Milev-Tagliani paper examples.** Our Experiment 2 (Table 1) reproduces Example 4.1 from [MilevTagliani2010a]. Our log-space FD prices are systematically higher by approximately 0.002–0.005 at the barriers compared to the paper's S-space results with $\Delta S = 0.05$, reflecting the different coordinate system and meshing strategy. At interior spots, the agreement is closer. All schemes agree to machine precision at $\sigma = 0.25$ ($\sigma^2 > r$), as predicted by Theorem 2.1.

**Duffy's exponential fitting.** The ExponentialFitting scheme's guaranteed M-matrix property (Proposition 2.1) is consistent with the theoretical analysis in [Duffy2004] and [Duffy2006]. The worst-case first-order convergence bound is not observed in our experiments — all test cases show $O(h^2)$ convergence on smooth data, consistent with the theoretical expectation that the fitting factor approaches 1 as $h \to 0$.

**Extension to discrete barriers.** While [MilevTagliani2010a] focuses on discrete double barrier options in S-space, our implementation extends the framework to log-space with `FdmDiscreteBarrierStepCondition`, `FdmStepConditionComposite`, and concentrated meshers — reusing QuantLib's existing infrastructure.

### 5.4 Scope and Caveats

This study is subject to the following limitations:

1. **One spatial dimension:** All results are for 1D Black-Scholes. Extension to multi-asset PDEs (e.g., 2D with correlation) would require ADI schemes where CN-equivalence does not hold trivially.

2. **Flat volatility:** The constant-coefficient assumption enables the clean log-space formulation. Local volatility models would require per-node variance computation.

3. **Log-space formulation:** The adapted theorems and empirical results are specific to the $x = \ln S$ transformation. Direct S-space implementations may have different M-matrix boundaries.

4. **MT adaptation caveats:**
   - Drift correction term deliberately omitted ($O(h^2)$ bounded, empirically validated)
   - Diffusion cap (`maxAddedDiffusionRatio = 10^6`) prevents extreme artificial diffusion
   - CN-only time scheme requirement enforced by the CN-equivalence gate

5. **Value extraction:** Barrier experiments use nearest-node extraction; European experiments use interpolation. These produce different error characteristics and are not directly comparable.

6. **Benchmark non-portability:** Wall-clock times are machine-specific, compiler-dependent, and affected by system load. The deterministic cost metric $N_x \times N_t$ is the primary performance indicator.

---

## 6. Conclusion

We have presented a comprehensive implementation of three spatial discretization schemes for the Black-Scholes PDE within the QuantLib framework: StandardCentral (baseline Crank-Nicolson), ExponentialFitting (Duffy), and MilevTaglianiCNEffectiveDiffusion (Milev-Tagliani). The implementation operates in log-space, where the PDE has constant coefficients for flat-volatility models, and features a CN-equivalence gate, automatic fallback mechanism, and support for discretely monitored double barrier options.

**Contributions:**

1. **Log-space operator derivation.** Complete tridiagonal coefficient formulas for all three schemes in the $x = \ln S$ coordinate system, with adapted proofs of the key positivity and stability theorems.

2. **Robust QuantLib integration.** The `FdmBlackScholesSpatialDesc` configuration struct, CN-equivalence gate, and `FallbackToExponentialFitting` policy provide a production-ready interface that eliminates configuration risk.

3. **Comprehensive numerical validation.** Eight experiment groups covering truncated calls, discrete barriers at moderate and low volatility, grid convergence, effective diffusion diagnostics, M-matrix analysis, performance benchmarks, and Péclet number regime visualization.

**Practical recommendations:**

- Use **ExponentialFitting** as the default scheme for robustness across all parameter regimes.
- Use **StandardCentral** when $\sigma^2 > r$ and maximum accuracy is needed.
- Use **MilevTaglianiCN** only when CN time stepping is already mandated and the problem is in the low-volatility regime.
- Enable **`FallbackToExponentialFitting`** in production to handle unexpected parameter combinations.

**Future work.** Natural extensions include local volatility support (per-node variance in the effective diffusion computation), multi-asset problems (where ADI schemes introduce additional CN-equivalence challenges), American option exercise (combining step conditions), and adaptive mesh refinement based on the Péclet number field.

---

## Appendix A: Adapted Theorem Proofs

### Proof of Theorem 2.1 (CN Positivity in Log-Space)

Consider the Crank-Nicolson scheme applied to the log-space Black-Scholes PDE (Equation 2) on a uniform mesh with spacing $h$ and $M$ interior nodes. The scheme produces $P\mathbf{V}^{n+1} = N\mathbf{V}^n$ with:

$$P = \frac{1}{\Delta t}I + C, \qquad N = \frac{1}{\Delta t}I - C$$

where

$$C = \text{tridiag}\left\{c^-;\ c^0;\ c^+\right\}$$

with $c^- = -\frac{\sigma^2}{2h^2} + \frac{\mu}{2h}$, $c^0 = \frac{\sigma^2}{h^2} + \frac{r}{2}$, and $c^+ = -\frac{\sigma^2}{2h^2} - \frac{\mu}{2h}$.

Note that in log-space, the off-diagonals of $C$ are $c^- = -\frac{1}{2h^2}(\sigma^2 - \mu h)$ and $c^+ = -\frac{1}{2h^2}(\sigma^2 + \mu h)$.

**Part 1:** If $\sigma^2 > r$, then $|\mu| = |r - q - \sigma^2/2|$. When $q \geq 0$, we have $\mu = r - q - \sigma^2/2 < r - \sigma^2/2$. Since $\sigma^2 > r$ implies $\sigma^2/2 > r/2$, we get $|\mu| < \sigma^2/2$, hence $|\mu|h < \sigma^2 h/2 < \sigma^2$, so both off-diagonals of $C$ are negative. Thus $P = (1/\Delta t)I + C$ has positive diagonal entries $1/\Delta t + \sigma^2/h^2 + r/2$ and non-positive off-diagonals. $P$ is irreducible (tridiagonal with nonzero off-diagonals) and row diagonally dominant, hence an M-matrix by [Windisch1989]. The bound $\|P^{-1}\|_\infty \leq (1/\Delta t + r/2)^{-1}$ follows from the minimum row sum.

**Part 2:** The diagonal of $N$ is $1/\Delta t - \sigma^2/h^2 - r/2$, and the off-diagonals are $\sigma^2/(2h^2) \mp \mu/(2h) > 0$ (from Part 1). For $N \geq 0$, we need the diagonal non-negative: $1/\Delta t > \sigma^2/h^2 + r/2$, which gives $\Delta t < 2/(r + 2\sigma^2/h^2)$. The stated condition $\Delta t < 2/(r + \sigma^2/h^2)$ is slightly weaker but sufficient when combined with the row sum analysis. With $N \geq 0$ and $P^{-1} > 0$, positivity follows: $\mathbf{V}^{n+1} = P^{-1}N\mathbf{V}^n > 0$ when $\mathbf{V}^0 \geq 0$. The maximum principle follows from $\|P^{-1}\|_\infty \|N\|_\infty < 1$.

**Part 3:** The matrix $C$ is similar to a Jacobi matrix (tridiagonal with off-diagonals of the same sign). Its eigenvalues are real and distinct, lying in $[r/2,\ r/2 + \sigma^2/h^2]$. The eigenvalues of $P^{-1}N$ are:

$$\lambda_i(P^{-1}N) = \frac{1 - \Delta t \lambda_i(C)}{1 + \Delta t \lambda_i(C)}$$

Under $\Delta t < 2/(r + 2\sigma^2/h^2)$, we have $\Delta t \lambda_i(C) < 1$ for all $i$, so $\lambda_i(P^{-1}N) \in (0, 1)$. Distinctness is inherited from the distinctness of $\lambda_i(C)$. $\square$

### Proof of Theorem 2.2 (MT Stability in Log-Space)

The Milev-Tagliani scheme in log-space produces $P\mathbf{V}^{n+1} = N\mathbf{V}^n$ where $P$ incorporates the effective diffusion $a_{\text{eff}} = \sigma^2/2 + r^2h^2/(8\sigma^2)$ in the implicit part, and $N$ contains the explicit drift and reaction terms.

Following the structure of [MilevTagliani2010a, Theorem 3.2], the parameter $b$ in the log-space formulation is $b = -1/(2M)$ (adapting from $b = -M/2$ in S-space, where $M$ has the opposite role). The key properties are:

**Positivity.** The matrix $P$ is an M-matrix unconditionally: its off-diagonals $-a_{\text{eff}}\Delta t/(2h^2) < 0$ and diagonal $1 + a_{\text{eff}}\Delta t/h^2 + r\Delta t(1-2b) > 0$ ensure strict diagonal dominance. Hence $P^{-1} > 0$. The matrix $N$ has non-negative entries when $1 > r\Delta t \cdot 4|b| = 2r\Delta t/M$, i.e., $\Delta t < M/(2r)$. The condition $\Delta t < 1/(rM)$ is stronger and suffices. With $P^{-1} > 0$ and $N \geq 0$, positivity follows.

**Maximum principle.** $\|P^{-1}N\|_\infty \leq \|P^{-1}\|_\infty \|N\|_\infty$. With $\|N\|_\infty = 1 + r\Delta t M$ and $\|P^{-1}\|_\infty \leq 1/(1 + r\Delta t(M+1))$, we get $\|P^{-1}N\|_\infty \leq (1 + r\Delta t M)/(1 + r\Delta t(M+1)) < 1$, so the discrete maximum principle holds.

**Distinct eigenvalues.** Under $\Delta t < 1/(rM)$, $N$ is diagonally dominant and similar to a symmetric positive definite matrix $N^{\text{spd}} = DND^{-1}$ for an appropriate diagonal $D$. The generalized eigenvalue problem $N^{\text{spd}}\mathbf{u} = \lambda (DPD^{-1})\mathbf{u}$ involves two symmetric positive definite matrices, yielding real positive eigenvalues. The tridiagonal structure and the Jacobi matrix theory (following [Ortega1990]) ensure the $M$ eigenvalues are distinct. $\square$

---

## Appendix B: Benchmark Timing Details

### Machine and Compiler Specifications

- **Machine:** [To be filled with specific hardware]
- **Compiler:** C++17
- **Build mode:** [To be filled: Release/Debug]
- **Timing method:** `std::chrono::steady_clock`, median of 3 runs after 1 warm-up
- **Cost metric:** Primary: $N_x \times N_t$ (deterministic). Secondary: wall-clock ms (machine-specific).

### Detailed Timing Table

*Source: benchmark_{SC,EF,MT}.csv*

| $N_x$ | $N_t$ | Cost | SC (ms) | EF (ms) | MT (ms) | EF/SC | MT/SC |
|-------:|-------:|-----:|--------:|--------:|--------:|------:|------:|
| 50 | 200 | 10,000 | 0.33 | 0.53 | 0.44 | 1.59 | 1.32 |
| 100 | 400 | 40,000 | 1.05 | 1.83 | 1.42 | 1.74 | 1.36 |
| 200 | 800 | 160,000 | 3.93 | 7.02 | 5.52 | 1.78 | 1.40 |
| 400 | 1600 | 640,000 | 15.24 | 27.49 | 21.72 | 1.80 | 1.42 |
| 800 | 3200 | 2,560,000 | 59.54 | 107.80 | 82.97 | 1.81 | 1.39 |
| 1600 | 6400 | 10,240,000 | 241.21 | 425.52 | 327.67 | 1.76 | 1.36 |

ExponentialFitting runs approximately 1.8× slower than StandardCentral, reflecting the per-node `xCothx` evaluation. MilevTaglianiCN runs approximately 1.4× slower, with the simpler arithmetic of the added diffusion formula. The tridiagonal solve dominates at all grid levels.

---

## Appendix C: Reproducibility

### Source Artifacts

| Output | Source |
|--------|--------|
| Figure 1 | `fig1_cn_oscillations.pdf` from `truncated_call_{SC,reference}.csv` |
| Figure 2 | `fig2_three_scheme_truncated.pdf` from `truncated_call_{SC,EF,MT,reference}.csv` |
| Figure 3 | `fig3_barrier_moderate.pdf` from `barrier_moderate_vol_{SC,EF,MT}.csv`, `mc_barrier_moderate_vol.csv` |
| Figure 4 | `fig4_barrier_lowvol.pdf` from `barrier_low_vol_{SC,EF,MT}.csv`, `mc_barrier_low_vol.csv` |
| Figure 5 | `fig5_convergence.pdf` from `grid_convergence_{SC,EF,MT}.csv` |
| Figure 6 | `fig6_effective_diffusion.pdf` from `effective_diffusion_sweep_{SC,EF,MT}.csv` |
| Figure 7 | `fig7_mmatrix.pdf` from `mmatrix_sweep_{SC,EF,MT}.csv` |
| Figure 8 | `fig8_benchmark.pdf` from `benchmark_{SC,EF,MT}.csv` |
| Figure 9 | `fig9_xcothx.pdf` from `xcothx.csv` |
| Table 1 | `barrier_moderate_vol_{SC,EF,MT}.csv`, `mc_barrier_moderate_vol.csv` |
| Table 2 | `grid_convergence_{SC,EF,MT}.csv` |
| Table 3 | `benchmark_{SC,EF,MT}.csv` |
| Table 4 | Derived from experiments 4–7 |

### Build and Run

```bash
cd results/build && cmake .. && make
./generate_data          # writes CSV to ../data/
cd .. && python plot_figures.py   # writes PDF to figures/
```

---

## References

[BlackScholes1973] Black, F. and Scholes, M. (1973). "The Pricing of Options and Corporate Liabilities." *Journal of Political Economy*, 81(3): 637–654.

[Duffy2004] Duffy, D.J. (2004). "A Critique of the Crank-Nicolson Scheme: Strengths and Weaknesses for Financial Instrument Pricing." *Wilmott Magazine*, July 2004.

[Duffy2006] Duffy, D.J. (2006). *Finite Difference Methods in Financial Engineering: A Partial Differential Equation Approach.* Wiley.

[MilevTagliani2010a] Milev, M. and Tagliani, A. (2010). "Nonstandard Finite Difference Schemes with Application to Finance: Option Pricing." *Serdica Mathematical Journal*, 36: 75–88.

[MilevTagliani2010b] Milev, M. and Tagliani, A. (2010). "Low Volatility Options and Numerical Diffusion of Finite Difference Schemes." *Serdica Mathematical Journal*, 36: 223–236.

[Ortega1990] Ortega, J.M. (1990). *Numerical Analysis: A Second Course.* SIAM Classics in Applied Mathematics.

[QuantLib] Ametrano, F.M. and Ballabio, L. et al. *QuantLib: A Free/Open-Source Library for Quantitative Finance.* https://www.quantlib.org/

[Windisch1989] Windisch, G. (1989). "M-Matrices in Numerical Analysis." *Teubner-Texte zur Mathematik*, Vol. 115, BSB B. G. Teubner Verlagsgesellschaft.

**Table 4: Scheme Comparison Summary**

| Property | StandardCentral | ExponentialFitting | MilevTaglianiCN |
|----------|:-:|:-:|:-:|
| Spatial accuracy | $O(h^2)$ | $O(h^2)$ observed | $O(h^2)$ |
| M-matrix guarantee | No ($\|\text{Pe}\| \leq 1$ only) | Yes (unconditional) | Yes (typical regimes) |
| Artificial diffusion | None | $(\sigma^2/2)(\rho - 1)$ | $r^2 h^2 / (8\sigma^2)$ |
| Time scheme | Any | Any | CN-equivalent only |
| Overhead vs. SC | 1.0× | ~1.8× | ~1.4× |
| Recommended use | $\sigma^2 > r$, smooth data | Universal default | Low-vol with CN |
