# Black-Scholes as Quantum Evolution: A Physics-Informed Analysis of Nonstandard Finite Difference Schemes

*A conceptual essay mapping the Black-Scholes finite difference framework — as implemented in the QuantLib Huatai fork — to paradigms from quantum mechanics and statistical physics.*

---

## Chapter 1: The Quantum Nature of Finance — Black-Scholes as Imaginary-Time Schrödinger

### 1.1 The Hidden Symmetry

The Black-Scholes equation describes how option price $V(S, t)$ evolves with the underlying asset price $S$ and calendar time $t$:

$$\frac{\partial V}{\partial t} + (r - q) S \frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} - rV = 0$$

where $r$ is the risk-free rate, $q$ the continuous dividend yield, and $\sigma$ the volatility. This equation looks ungainly — variable coefficients $S$ and $S^2$ multiplying every derivative. But there is a hidden symmetry.

A stock dropping from 100 to 50 inflicts the same percentage loss as one dropping from 10 to 5. This **scale invariance** tells us that the natural coordinate is not the price $S$ but its logarithm $x = \ln S$. And because options have a terminal payoff at maturity $T$, their value is determined by backward induction — so the natural time variable is time-to-maturity $\tau = T - t$, pointing the arrow of time from future payoff toward present value.

### 1.2 The Coordinate Transformation

Under $x = \ln S$ and $\tau = T - t$, with $U(x, \tau) = V(S, t)$, the chain rule gives **(exact)**:

$$\frac{\partial U}{\partial \tau} = \frac{\sigma^2}{2} \frac{\partial^2 U}{\partial x^2} + \mu \frac{\partial U}{\partial x} - rU$$

where $\mu = r - q - \sigma^2/2$ is the risk-neutral drift in log-space. The variable coefficients have vanished. What remains is a constant-coefficient convection-diffusion-reaction equation — already far more amenable to physical interpretation.

The derivation is straightforward: $S \partial_S V = \partial_x U$ from the chain rule, and $S^2 \partial_{SS} V = \partial_{xx} U - \partial_x U$ from differentiating again. Substituting and reversing the time arrow via $\partial_t V = -\partial_\tau U$ yields the result.

In the language of physics, this is a **diffusion equation** with three terms: pure diffusion $(\sigma^2/2) U_{xx}$, advection bias $\mu U_x$ (the drift of the risk-neutral random walk), and amplitude decay $-rU$ (the time value of money).

### 1.3 The Gauge Transformation — Stripping to Pure Diffusion

To reveal the deepest structure, we perform a **gauge transformation** **(exact)** — the physicist's tool for removing inessential complications from a wave equation. Set:

$$U(x, \tau) = e^{-\alpha x - \beta \tau} \, \Psi(x, \tau)$$

Computing derivatives and substituting into the log-space BS equation, we require the drift and decay terms to vanish. This determines the gauge parameters uniquely:

$$\alpha = \frac{\mu}{\sigma^2} = \frac{r - q}{\sigma^2} - \frac{1}{2}$$

$$\beta = r + \frac{\mu^2}{2\sigma^2} = r + \frac{(r - q - \sigma^2/2)^2}{2\sigma^2}$$

The parameter $\alpha$ removes the advection bias — it shifts to a co-moving frame traveling with the drift. The parameter $\beta$ absorbs both the discount rate and the residual kinetic energy of that drift. What survives is the purest possible evolution:

$$\frac{\partial \Psi}{\partial \tau} = \frac{\sigma^2}{2} \frac{\partial^2 \Psi}{\partial x^2}$$

This is the **heat equation** — one of the most studied objects in mathematical physics. No drift, no decay, no variable coefficients. Pure diffusion.

### 1.4 The Schrödinger Correspondence

Now consider the free-particle Schrödinger equation in one dimension:

$$i\hbar \frac{\partial \psi}{\partial t} = -\frac{\hbar^2}{2m} \frac{\partial^2 \psi}{\partial x^2}$$

Under the **Wick rotation** $t \to -i\tau$ (formally replacing real time with imaginary time), the $i$ on the left side absorbs into the time derivative, and we obtain the **imaginary-time Schrödinger equation**:

$$\frac{\partial \psi}{\partial \tau} = \frac{\hbar}{2m} \frac{\partial^2 \psi}{\partial x^2}$$

Comparing with our gauge-transformed Black-Scholes equation, the identification is immediate **(exact identity)**:

$$\frac{\sigma^2}{2} \longleftrightarrow \frac{\hbar}{2m_{\text{eff}}}$$

Setting $\hbar = 1$ (natural units — the market's intrinsic randomness plays the role of Planck's constant), we obtain:

$$m_{\text{eff}} = \frac{1}{\sigma^2}$$

**A low-volatility underlying is a heavy quantum particle.** Its probability distribution diffuses slowly — enormous "news shocks" are required to shift its state. A major equity index with $\sigma = 0.15$ has effective mass $m_{\text{eff}} \approx 44$. A volatile altcoin with $\sigma = 1.0$ has effective mass $m_{\text{eff}} = 1$ — practically massless, buffeted by every perturbation.

This is not metaphor. The gauge-transformed BS equation and the imaginary-time free-particle Schrödinger equation are **the same mathematical object**. The solutions are identical under the coefficient identification. What differs is the physical interpretation — and, critically, the absence of the imaginary unit $i$.

### 1.5 Why No $i$? — Diffusion vs. Unitary Evolution

The real-time Schrödinger equation has $i$ on the left side: $i\hbar\psi_t = H\psi$. This $i$ makes the evolution operator $e^{-iHt/\hbar}$ **unitary** — it preserves the norm of the wave function, $\|\psi\|^2 = 1$, forever. Quantum probability is conserved. The wave function oscillates but never decays.

The Black-Scholes equation has no $i$. The evolution operator $e^{-H\tau}$ (with the financial "Hamiltonian" $H = -\sigma^2/(2)\partial_{xx}$) is a **contraction semigroup** — it shrinks the solution. Option values diffuse and decay; they do not oscillate.

In quantum mechanics, probability is the squared modulus of the wave function: $\rho = |\psi|^2$. The $i$ ensures that $\psi$ evolves unitarily and $\rho$ satisfies a continuity equation. In finance, the option price $V$ is itself the probability-weighted expected value — there is no need to square anything. The risk-neutral expectation directly gives the price.

This distinction — **unitary evolution vs. diffusive decay** — is the deepest point where the analogy reaches its limit **(formal analogy)**. The mathematical structures are identical; the physical interpretations diverge at the meaning of the "wave function."

### 1.6 The Financial Hamiltonian — Before Gauge Transformation

Before the gauge strip, the log-space BS equation reads:

$$\partial_\tau U = \underbrace{\frac{\sigma^2}{2} \partial_{xx} U}_{\text{kinetic energy}} + \underbrace{\mu \, \partial_x U}_{\text{advection bias}} - \underbrace{r \, U}_{\text{potential / discount}}$$

In quantum language **(formal analogy)**:

- The **kinetic energy** operator $-(\sigma^2/2)\partial_{xx}$ governs diffusion — the quantum particle spreading through space.
- The **advection bias** $\mu \partial_x$ acts like a constant force field, giving the particle an average drift velocity. (More precisely, it is advection in the moving frame, not a conservative force — a subtlety that matters for self-adjointness; see Chapter 5.)
- The **discount term** $-rU$ is analogous to a constant scalar potential $V_0 = r$ — it does not change the shape of the probability distribution, only the overall amplitude. This is the time value of money: total option value decays exponentially at rate $r$ as time passes.

Together, these define the **financial Hamiltonian**:

$$\hat{H}_{\text{BS}} = -\frac{\sigma^2}{2} \partial_{xx} - \mu \, \partial_x + r$$

The log-space BS equation is then $\partial_\tau U = -\hat{H}_{\text{BS}} U$ — **exactly the imaginary-time Schrödinger equation with a constant-potential Hamiltonian** **(exact)** for the Euclidean interpretation.

### 1.7 Feynman-Kac: The Exact Path Integral

The deepest connection between finance and quantum physics is not an analogy at all — it is a theorem.

**Feynman-Kac Theorem (exact).** Let $X_t$ follow $dX_t = \mu \, dt + \sigma \, dW_t$ under the risk-neutral measure $\mathbb{Q}$. Then the solution to

$$V_\tau = \frac{\sigma^2}{2} V_{xx} + \mu V_x - rV, \qquad V(x, 0) = f(x)$$

is given by:

$$V(x, \tau) = \mathbb{E}^{\mathbb{Q}}_x\!\left[e^{-r\tau} f(X_\tau)\right]$$

The proof is elegant: define $Y_s = e^{-rs} V(X_s, \tau - s)$ and apply Itô's lemma. The PDE ensures that the drift of $Y$ vanishes, making $Y$ a martingale. Evaluating at $s = 0$ and $s = \tau$ gives the result.

This is the financial incarnation of the **Feynman path integral**. In Euclidean quantum mechanics, the propagator is:

$$K(x', \tau | x) = \int_{x(0) = x}^{x(\tau) = x'} \mathcal{D}x \, e^{-S_E[x]}$$

where $S_E$ is the Euclidean action. For the Black-Scholes case **(exact via Wiener measure)**:

$$S_E[x] = \int_0^\tau \left[\frac{(\dot{x} - \mu)^2}{2\sigma^2} + r\right] ds$$

The first term is the "kinetic energy" — paths that deviate from the drift $\mu$ are exponentially suppressed, weighted by the inverse effective mass $\sigma^2$. The second term is the "potential energy" — the time value of money discount applied uniformly along each path.

Every possible price trajectory contributes. Paths close to the risk-neutral drift dominate; wild trajectories are exponentially suppressed by the Gaussian weight. The option price at time zero is the sum over all these weighted outcomes — **a Euclidean path integral over the space of price trajectories**.

The Feynman-Kac theorem makes this rigorous, not as an analogy but as an exact mathematical identity. The path integral, constructed via Wiener measure and time-sliced Gaussian kernels, converges to the PDE solution. The physicist's notation $\int \mathcal{D}x \, e^{-S_E}$ is formal shorthand **(formal)** for this exact construction.

### 1.8 Risk-Neutral Measure and Wick Rotation

In quantum field theory, the **Wick rotation** $t \to -i\tau$ transforms the oscillatory Minkowski path integral $\int \mathcal{D}\phi \, e^{iS}$ into the convergent Euclidean path integral $\int \mathcal{D}\phi \, e^{-S_E}$. This is an analytic continuation in the time variable.

In finance, the **Girsanov theorem** transforms the physical probability measure $\mathbb{P}$ (where stocks have expected return $\mu_{\text{phys}}$) into the risk-neutral measure $\mathbb{Q}$ (where stocks drift at $r - q$). This is a change of probability measure — a shift in the drift of the Brownian motion.

Both operations achieve a similar structural result: they transform an ill-behaved (oscillatory or drift-dependent) integral into a well-behaved (convergent, drift-standardized) one. But the mechanisms are fundamentally different **(formal analogy)**:

- Wick rotation is an **analytic continuation** in the time variable.
- Girsanov is a **change of probability measure** — a reweighting of paths, not a rotation of time.

The resemblance is striking but not exact. Risk-neutral pricing is not Wick rotation; it is a separate mathematical operation that happens to produce analogous structures. This is one of the boundaries where the quantum-finance correspondence, however deep, must be stated with care.

---

## Chapter 2: The Discrete Hamiltonian — Spatial Schemes as Lattice Physics

### 2.1 Discretizing the Continuum

Having established that the Black-Scholes equation is an imaginary-time Schrödinger equation, we now face the practical question: how do we solve it numerically?

In quantum physics, moving from continuous space to a discrete lattice is a profound step. The continuum Hamiltonian $\hat{H} = -(\hbar^2/2m)\partial_{xx}$ becomes a matrix — a tridiagonal operator on lattice sites. The physics of the lattice approximation depends critically on how we discretize, and different discretization strategies correspond to different physical insights about the underlying dynamics.

The QuantLib Huatai fork implements three spatial discretization schemes for the log-space BS operator. Each one corresponds to a distinct strategy from computational physics.

### 2.2 StandardCentral — The Naive Lattice

The simplest discretization replaces derivatives with centered finite differences on a mesh $\{x_0, x_1, \ldots, x_M\}$ with spacing $h$:

$$\partial_{xx} U \approx \frac{U_{i+1} - 2U_i + U_{i-1}}{h^2}, \qquad \partial_x U \approx \frac{U_{i+1} - U_{i-1}}{2h}$$

The resulting tridiagonal "Hamiltonian matrix" has entries:

$$a_{i,i-1} = \frac{\sigma^2}{2h^2} - \frac{\mu}{2h}, \qquad a_{i,i} = -\frac{\sigma^2}{h^2} - r, \qquad a_{i,i+1} = \frac{\sigma^2}{2h^2} + \frac{\mu}{2h}$$

In the lattice physics analogy **(formal)**, this is the **tight-binding Hamiltonian** — each lattice site couples only to its nearest neighbors, with hopping amplitudes $\sigma^2/(2h^2) \pm \mu/(2h)$ encoding the interplay between diffusion and drift.

This is the scheme assembled in the `StandardCentral` branch of `FdmBlackScholesOp::setTime()` (`fdmblackscholesop.cpp:107–152`), where `mapT_.axpyb(drift, dxMap_, dxxMap_.mult(0.5*v), Array(1, -r))` constructs the operator $L = \mu \partial_x + (\sigma^2/2) \partial_{xx} - r$.

**The problem appears when convection dominates diffusion.** The off-diagonal $a_{i,i-1} = \sigma^2/(2h^2) - \mu/(2h)$ becomes negative when $|\mu|h > \sigma^2$ — that is, when the mesh **Péclet number**

$$\text{Pe} = \frac{\mu h}{\sigma^2}$$

exceeds 1 in absolute value. In this regime, the lattice Hamiltonian has the wrong sign structure: the "hopping amplitude" in one direction becomes negative, which is unphysical — it means the discrete operator no longer preserves positivity of the solution.

For a physicist, this is immediately recognizable: it is the **lattice artifact** that arises when the grid is too coarse to resolve the convection scale. The centered difference approximation to $\partial_x$ introduces spurious oscillation modes that the physical system does not have.

### 2.3 ExponentialFitting — The Scharfetter-Gummel Flux

In semiconductor device physics, Scharfetter and Gummel (1969) encountered the same problem: the drift-diffusion equation for carrier transport,

$$J = -D \frac{dn}{dx} + v \cdot n$$

needs to be discretized on meshes where the electric-field-driven drift $v$ can dominate the thermal diffusion $D$. Their solution was to use the **exact steady-state solution** on each mesh cell to construct the numerical flux.

On a cell $[x_i, x_{i+1}]$ with constant $D$ and $v$, the steady-state current $J$ is constant, and the carrier density follows an exponential profile. Imposing this exact local solution as the basis for the numerical flux yields an effective diffusion coefficient **(exact equivalence)**:

$$D_{\text{eff}} = D \cdot \text{Pe} \cdot \coth(\text{Pe}), \qquad \text{Pe} = \frac{v h}{2D}$$

This is the **same algebraic construction** as the exponential fitting prescription used in the QuantLib implementation — equivalent up to the Péclet number convention (SG uses $\text{Pe} = vh/(2D)$ while BS uses $\text{Pe} = \mu h/\sigma^2$; the two coincide under $D = \sigma^2/2$, $v = \mu$). For the Black-Scholes operator:

$$a_{\text{eff}} = \frac{\sigma^2}{2} \cdot \rho(\text{Pe}), \qquad \rho(\text{Pe}) = \text{Pe} \cdot \coth(\text{Pe}), \qquad \text{Pe} = \frac{\mu h}{\sigma^2}$$

The connection is not merely analogous — it is the **same algebraic construction** applied to the convection-diffusion part of the operator **(exact algebraic equivalence for 1D cellwise convection-diffusion)**. The reaction term $-rU$ is separate and handled independently.

In the code (`fdmblackscholesop.cpp:51–57`):

```cpp
case FdmBlackScholesSpatialDesc::Scheme::ExponentialFitting: {
    const Real Pe = drift[i] * h[i] / vEff;
    const Real rho = detail::xCothx(Pe, desc.peSmall, desc.peLarge);
    aUsed[i] = aBase * rho;
    break;
}
```

The function `xCothx` (`fdmhyperboliccot.hpp:25–39`) evaluates $\rho(\text{Pe}) = \text{Pe}/\tanh(\text{Pe})$ with three-regime numerical stability:

- **Small $|\text{Pe}| < \text{peSmall}$**: Taylor series $1 + \text{Pe}^2/3$ (avoids 0/0).
- **Moderate $\text{peSmall} \leq |\text{Pe}| \leq \text{peLarge}$**: Direct computation via `std::tanh`.
- **Large $|\text{Pe}| > \text{peLarge}$**: Asymptotic form $|\text{Pe}|$ (avoids overflow in $e^{2\text{Pe}}$).

The transition thresholds `peSmall` ($10^{-6}$) and `peLarge` ($50$) are configurable parameters in `FdmBlackScholesSpatialDesc`, not hard-coded constants.

In physics language, the three regimes correspond to: diffusion-dominated (thermal equilibrium), intermediate (mixed transport), and convection-dominated (ballistic transport). The $x \coth(x)$ function appears throughout physics — in the Langevin function for paramagnetism, in Bose-Einstein statistics for phonon heat capacity, and in the Debye model. Its appearance here, mediating between diffusive and convective transport regimes, is the same mathematical structure in a different physical context **(formal analogy)**.

**The key property**: $x \coth(x) \geq |x|$ for all real $x$, with equality only at $|x| \to \infty$. This guarantees that the fitted off-diagonals

$$\frac{a_{\text{eff}}}{h^2} \pm \frac{\mu}{2h} = \frac{\sigma^2}{2h^2}\big(\rho(\text{Pe}) \pm \text{Pe}\big) \geq 0$$

are non-negative **unconditionally** — the M-matrix property holds for all parameter values, all mesh spacings, all volatility regimes.

The proof in `results/paper/main.markdown` (Proposition 2.1) makes this rigorous: since $\rho(\text{Pe}) \geq |\text{Pe}|$, the lower off-diagonal satisfies $A^-_i = (\sigma^2/2h^2)(\rho - \text{Pe}) \geq 0$, and similarly for the upper off-diagonal. The CN implicit matrix $P = (1/\Delta t)I - (1/2)A$ then has positive diagonal, non-positive off-diagonals, and strict diagonal dominance — a nonsingular M-matrix with $P^{-1} > 0$.

### 2.4 MilevTaglianiCNEffectiveDiffusion — Artificial Viscosity

The third scheme takes a different approach, rooted in the computational physics tradition of **artificial viscosity** **(formal analogy)** — adding controlled numerical diffusion to stabilize a scheme that would otherwise oscillate.

In the original S-space formulation (Milev-Tagliani, *Serdica* 36:75–88, 2010), the authors discretize the reaction term $-rV$ using a 6-point stencil with parameters $\omega_1 = \omega_2 = -r/(16\sigma^2)$, which introduces an artificial diffusion of magnitude $(1/8)(r \Delta S / \sigma)^2 V_{SS}$. This is their nonstandard scheme, proven to preserve positivity under a time-step restriction $\Delta t < 1/(rM)$ (Theorem 3.2).

The QuantLib log-space adaptation translates this artificial diffusion:

$$a_{\text{eff}} = \frac{\sigma^2}{2} + \frac{r^2 h^2}{8\sigma^2}$$

subject to a safety cap: $a_{\text{add}} = \min(r^2 h^2/(8\sigma^2), \, R_{\max} \cdot \sigma^2/2)$ where $R_{\max}$ defaults to $10^6$.

In the code (`fdmblackscholesop.cpp:58–66`):

```cpp
case FdmBlackScholesSpatialDesc::Scheme::MilevTaglianiCNEffectiveDiffusion: {
    Real aAdd = r * r * h[i] * h[i] / (8.0 * vEff);
    aAdd = std::min(aAdd, desc.maxAddedDiffusionRatio * aBase);
    aUsed[i] = aBase + aAdd;
    break;
}
```

**Paper vs. Code — The Dual Track.** The full S-space-to-log-space translation yields both a diffusion addition and a drift correction $-r^2 h^2/(8\sigma^2)$ on the first-derivative term. The shipped implementation **deliberately omits the drift correction**: the operator assembles `mapT_.axpyb(drift, dxMap_, dxxMap_.mult(aUsed), Array(1, -r))` where `drift` is the unmodified physical drift $\mu$, not $\mu - a_{\text{add}}$. The code comment explains that "the drift correction is bounded by $O(h^2)$ grid error and validated by `testMilevTaglianiDriftCorrectionAudit()`." This means the shipped operator is a **diffusion-only adaptation** — it solves a slightly different PDE than the paper's full S-space scheme. In low-volatility/high-rate regimes where $a_{\text{add}}$ can be comparable to $\mu$, the omitted drift correction is not negligible in absolute terms, but remains $O(h^2)$ relative to the grid discretization error. The scheme should be understood as *inspired by* rather than *identical to* the Milev-Tagliani nonstandard scheme.

**Paper result (S-space, Theorem 3.2):** Under $\sigma^2 > r$ and $\Delta t < 1/(rM)$, with the full 6-point stencil, both $P^{-1} > 0$ and $N \geq 0$, guaranteeing positivity and the discrete maximum principle. The eigenvalues of $P^{-1}N$ are positive and distinct.

**Code reality (log-space):** The M-matrix condition depends on the mesh Péclet number $|\text{Pe}| = |\mu|h/\sigma^2 \leq 1$, which is a joint condition on parameters and mesh spacing — not the simple $\sigma^2 > r$ of the S-space result. The time-step constraint $\Delta t < 1/(rM)$ from the paper is noted but **not enforced at runtime** (`fdmblackscholessolver.cpp:55–58`).

**The CN-equivalence gate.** The MT scheme requires Crank-Nicolson time stepping. The solver enforces this (`fdmblackscholessolver.cpp:59–76`): if the time scheme is not CN-equivalent (CrankNicolson, Douglas, or CraigSneyd with $\theta = 0.5$) or if `dampingSteps > 0`, the solver silently falls back to ExponentialFitting. This is a **runtime safety mechanism**, not a theoretical result — it prevents the MT scheme from being used with incompatible time steppers that would break the algebraic structure on which its positivity properties depend.

In the physics analogy, artificial viscosity is a well-understood technique: when a numerical scheme produces unphysical oscillations (analogous to Gibbs phenomenon, or to lattice artifacts in quantum chromodynamics), one adds a controlled amount of dissipation to smooth them out. The MT scheme's artificial diffusion $r^2 h^2/(8\sigma^2)$ grows quadratically with mesh spacing — it is significant only on coarse grids where the convection scale is underresolved, and vanishes as $h \to 0$, preserving the formal convergence order.

### 2.5 The M-Matrix Property — Positivity as a Physical Requirement

In quantum mechanics, the probability density $\rho = |\psi|^2$ must be non-negative everywhere. A numerical scheme that produces negative probabilities is unphysical — it signals that the lattice approximation has broken down.

In finance, option prices must be non-negative (you cannot owe money for holding a call option). The discrete analog of this requirement is the **M-matrix property** of the implicit-side matrix $P$: if $P$ is an M-matrix (positive diagonal, non-positive off-diagonals, $P^{-1} \geq 0$), then a non-negative input produces a non-negative output at each time step **(exact for the discrete system)**.

The M-matrix property is **not** unitarity — it is **positivity preservation** **(formal analogy)**. In quantum mechanics, the evolution preserves the norm ($\|\psi\|^2 = 1$); in discrete BS, the M-matrix property preserves the sign ($V_i \geq 0$ for all $i$). Both are conservation laws of the discrete system, but they conserve different things.

The runtime diagnostic in `fdmmatrixdiagnostic.hpp` performs a **partial** M-matrix check: it scans off-diagonal entries of the assembled tridiagonal operator for sign violations and non-finite values, typically excluding boundary rows. This is a necessary condition for the M-matrix property (non-positive off-diagonals), but it does not constitute a full M-matrix proof (which additionally requires the inverse to be entry-wise non-negative). When violations are detected, the `FallbackToExponentialFitting` policy (`fdmblackscholesop.cpp:247–254`) re-assembles the operator using the exponential fitting scheme — a **runtime adaptive stabilization mechanism** **(heuristic)**. This is not a phase transition or symmetry breaking; it is an engineering safety net that detects when the chosen scheme has failed to maintain positivity guarantees and switches to one that always succeeds.

---

## Chapter 3: The Cayley Propagator — Time Evolution in the Financial Hilbert Space

### 3.1 From Differential Equation to Discrete Propagator

Having discretized space (the Hamiltonian), we must now discretize time (the evolution). In quantum mechanics, the time evolution operator for a time-independent Hamiltonian is:

$$U(t) = e^{-iHt/\hbar}$$

For small time steps $\Delta t$, the **Cayley form** provides a rational approximation:

$$U_{\text{Cayley}} = \left(I + \frac{i\Delta t H}{2\hbar}\right)^{-1}\left(I - \frac{i\Delta t H}{2\hbar}\right)$$

When $H$ is self-adjoint, this approximation is **exactly unitary**: $U_{\text{Cayley}}^\dagger U_{\text{Cayley}} = I$. This makes it the preferred time stepper in computational quantum mechanics — it preserves probability exactly at each step, regardless of step size.

### 3.2 Crank-Nicolson as the Financial Cayley Propagator

The Crank-Nicolson scheme for $V_\tau = AV$ (where $A$ is the spatial operator matrix) produces:

$$PV^{n+1} = NV^n, \qquad P = \frac{1}{\Delta t}I - \frac{1}{2}A, \qquad N = \frac{1}{\Delta t}I + \frac{1}{2}A$$

The propagator is:

$$G_{\text{CN}} = P^{-1}N = \left(I - \frac{\Delta t}{2}A\right)^{-1}\left(I + \frac{\Delta t}{2}A\right)$$

This is **exactly the Cayley transform** of $(\Delta t/2)A$ **(exact algebraic identity)**. The algebraic structure is identical to the QM Cayley propagator — the same rational map applied to the generator of evolution.

In the QuantLib implementation (`cranknicolsonscheme.cpp:37–45`), CN is decomposed into an explicit half-step and an implicit half-step:

```cpp
void CrankNicolsonScheme::step(array_type& a, Time t) {
    if (theta_ != 1.0)
        explicit_->step(a, t, 1.0-theta_);
    if (theta_ != 0.0)
        implicit_->step(a, t, theta_);
}
```

With $\theta = 0.5$ (standard CN), this implements $(I - (\Delta t/2)A)^{-1}(I + (\Delta t/2)A)$ — the explicit step applies $(I + (\Delta t/2)A)$ and the implicit step solves $(I - (\Delta t/2)A)$.

### 3.3 Contractivity, Not Unitarity

Here the analogy must be stated carefully **(formal analogy only)**.

In QM, the generator $A_{\text{QM}} = -(i/\hbar)H$ is **skew-adjoint** when $H$ is self-adjoint. This is what makes the Cayley transform exactly unitary: $G^\dagger G = I$.

The Black-Scholes spatial operator $A$ is **not** skew-adjoint. It has:
- A symmetric diffusion part: $(\sigma^2/2)\partial_{xx}$
- An antisymmetric convection part: $\mu \partial_x$
- A scalar dissipation: $-r$

The convection term breaks the symmetry, and the dissipation shifts all eigenvalues into the left half-plane. As a result, the CN propagator for BS is **contractive**, not unitary:

$$|\lambda_i(G_{\text{CN}})| \leq 1 \quad \text{for all } i$$

with **strict** contraction $|\lambda_i| < 1$ when $r > 0$. The spectral radius satisfies:

$$|\lambda| = \left|\frac{1 + z\lambda_A}{1 - z\lambda_A}\right|, \qquad z = \frac{\Delta t}{2}$$

where $\lambda_A = a + ib$ is an eigenvalue of $A$. When $a \leq 0$ (left half-plane), $|\lambda| \leq 1$.

So the CN scheme for BS **damps** all modes — consistent with the physics of diffusion. Option values decay toward their long-time steady state. There is no conservation of "probability" in the quantum sense; instead, there is monotone contraction — each time step reduces the supremum norm of the solution when the M-matrix condition holds.

### 3.4 Eigenvalue Pathology — Parasitic Checkerboard Modes

The Milev-Tagliani paper (Theorems 3.1–3.2) reveals a subtle pathology in the CN eigenvalue structure. Writing $P = (1/\Delta t)I + C$ and $N = (1/\Delta t)I - C$ with $C = -A/2$, the eigenvalues of $P^{-1}N$ are:

$$\lambda_i(P^{-1}N) = \frac{1 - \Delta t \lambda_i(C)}{1 + \Delta t \lambda_i(C)}$$

Under the conditions of Theorem 3.1 ($|\text{Pe}| < 1$ in log-space, with appropriate $\Delta t$ bounds), these eigenvalues lie in $(0, 1)$ — all positive, all less than 1. The scheme is stable, contractive, and oscillation-free.

But when the M-matrix condition is violated ($|\text{Pe}| > 1$), the explicit-side matrix $N$ develops negative entries. The eigenvalues $\lambda_i(C)$ can become large enough that

$$\lambda_i(P^{-1}N) \to -1$$

These are **parasitic high-frequency checkerboard modes** — the discrete solution alternates sign at adjacent grid points, producing the spurious oscillations visible in Figure 1 of the Milev-Tagliani paper. The eigenvalues remain inside the unit disk (the scheme is formally stable), but the near-$(-1)$ modes decay so slowly that they persist for many time steps, contaminating the numerical solution with visible oscillations.

**Scoping the pathology.** This occurs specifically in the **low-volatility / convection-dominated regime** where $\sigma^2$ is small relative to $r$ (or more precisely, where $|\mu|h > \sigma^2$). It affects:

- **StandardCentral + CN**: always, when $|\text{Pe}| > 1$
- **MT + CN under paper assumptions**: the paper's nonstandard scheme is designed to prevent this, but only under its stated conditions ($\sigma^2 > r$ in S-space, $\Delta t < 1/(rM)$)
- **ExponentialFitting + CN**: never, because the fitting factor ensures non-negative off-diagonals unconditionally

In the QuantLib implementation, the CN-equivalence gate (`fdmblackscholessolver.cpp:33–45`) identifies whether the time scheme reduces to CN in 1D:

```cpp
bool FdmBlackScholesSolver::isCrankNicolsonEquivalent1D(
        const FdmSchemeDesc& schemeDesc, Real tol) {
    switch (schemeDesc.type) {
      case FdmSchemeDesc::CrankNicolsonType:
      case FdmSchemeDesc::DouglasType:
      case FdmSchemeDesc::CraigSneydType:
        return std::fabs(schemeDesc.theta - 0.5) <= tol;
      default:
        return false;
    }
}
```

Douglas and Craig-Sneyd reduce to CN in 1D because the mixed-derivative operator `apply_mixed()` returns zero. When the MT scheme is requested with a non-CN-equivalent time scheme or with `dampingSteps > 0`, the solver falls back to ExponentialFitting — ensuring that the eigenvalue pathology is never triggered by an incompatible time-space combination.

---

## Chapter 4: Barriers as Quantum Measurements — Discrete Monitoring as Projection

### 4.1 The Knock-Out Event

A discretely monitored double barrier knock-out option works as follows: at each monitoring date $t_i$, the underlying price $S$ is observed. If $S$ lies outside the corridor $[L, U]$, the option is immediately terminated and a rebate is paid. Between monitoring dates, the option value evolves according to the Black-Scholes PDE on the full domain.

In the PDE framework, this is implemented by **resetting the solution at monitoring dates**: at time $t_i$, the value function is multiplied by the indicator function $\mathbf{1}_{[L,U]}$, setting the value to the rebate outside the corridor.

$$V(S, t_i) = V(S, t_i^-) \cdot \mathbf{1}_{[\ln L, \ln U]}(x) + \text{rebate} \cdot \mathbf{1}_{\mathbb{R} \setminus [\ln L, \ln U]}(x)$$

### 4.2 The Code Mechanism

In the QuantLib implementation (`fdmdiscretebarrierstepcondition.cpp:62–96`), the barrier condition is applied through a step condition that fires at monitoring times:

```cpp
void FdmDiscreteBarrierStepCondition::applyTo(Array& a, Time t) const {
    // ...tolerant time matching...
    for (Size idx : knockOutIndices_)
        a[idx] = rebate_;
}
```

The `knockOutIndices_` are precomputed in the constructor (`fdmdiscretebarrierstepcondition.cpp:55–59`): every mesh node with $x < \ln L$ or $x > \ln U$ is marked for knockout. At each monitoring time, the solution array at these nodes is **overwritten** with the rebate value. This is not a smooth modification — it is a sharp, discontinuous reset.

The time matching uses a tolerant comparison (`fdmdiscretebarrierstepcondition.cpp:71–87`) with relative tolerance $10^{-10}$, robust against floating-point representation differences in monitoring times. The mesh is constructed with concentration points at the barriers: `cPoints = {{K,0.1,true},{L,0.1,true},{U,0.1,true}}` ensures that the strike and both barriers land exactly on mesh nodes, providing maximum resolution at the discontinuities.

### 4.3 Absorption and Killing — The Classical Interpretation

The standard numerical analysis interpretation is **absorption** or **killing**: the barrier is an absorbing boundary. When the particle (price) reaches the barrier, it is removed from the system and replaced with a fixed value (the rebate). This is equivalent to solving the PDE on a domain with Dirichlet boundary conditions at the barriers, except that the boundaries are imposed discretely in time rather than continuously.

In the language of stochastic processes, this is a **killed diffusion**: the process $X_t$ is terminated at the first monitoring time where $X_{t_i} \notin [\ln L, \ln U]$, and the terminal payoff is replaced by the rebate.

### 4.4 Projective Measurement — The Quantum Interpretation

In quantum mechanics, a **projective measurement** collapses the wave function. If we measure whether a particle is inside or outside a region $[L, U]$, the **von Neumann projection postulate** says **(formal analogy)**:

$$|\psi\rangle \to \frac{P_{[L,U]} |\psi\rangle}{\| P_{[L,U]} |\psi\rangle \|}$$

where $P_{[L,U]}$ is the projection operator onto states localized in $[L, U]$.

The discrete barrier condition performs precisely this operation (without the normalization, since in finance we track value, not probability):

1. **Before monitoring**: The solution $V(x, t_i^-)$ exists on the entire domain — the "wave function" is spread across all possible prices.
2. **At monitoring**: The solution is projected: $V(x, t_i) = V(x, t_i^-) \cdot \mathbf{1}_{[\ln L, \ln U]}(x) + \text{rebate} \cdot \mathbf{1}_{\text{outside}}(x)$.
3. **After monitoring**: Evolution continues from the projected state.

This is **repeated quantum measurement**: at each monitoring date, we "observe" whether the price is inside the corridor. If it is not, the option "collapses" to the rebate value. The information gained by the measurement (the price is outside the corridor) irreversibly alters the state.

The analogy extends to the **quantum Zeno effect**: as monitoring frequency increases (continuous monitoring limit), the barrier becomes effectively continuous — the option is killed the instant the price exits the corridor. In the Zeno limit, the absorption is instantaneous and complete, and the discrete barrier option converges to the continuously monitored barrier option.

**Important caveat** **(formal analogy only)**: The financial "measurement" is not quantum measurement in any physical sense. There is no superposition of being "knocked out" and "alive" — the barrier condition is a classical observation. The projection analogy captures the mathematical structure (sharp state reset at discrete times) but not the physical content (wave function collapse, Born rule, entanglement). The absence of interference effects is the clearest marker: in quantum mechanics, the projection interacts with the phase of the wave function; in finance, there is no phase — only non-negative values.

### 4.5 Why Not Tunneling?

A tempting but incorrect analogy is **quantum tunneling** — a particle passing through a potential barrier that classically forbids it. But the barrier option mechanism is the opposite:

- **Tunneling**: the particle penetrates the barrier and continues on the other side. The barrier is an energy obstacle.
- **Knock-out**: the option is **destroyed** when the price reaches the barrier. The barrier is a killing surface, not a potential wall.

The correct physics analog is **absorption at a boundary** or, equivalently, **projective measurement that collapses the state**. The barrier does not repel the price; it annihilates the option.

---

## Chapter 5: Where the Analogy Breaks — Honest Limits of the Correspondence

### 5.1 Three Layers of Correspondence

Throughout this essay, we have annotated each analogy with its strength. The correspondence between Black-Scholes and quantum mechanics operates at three levels:

1. **Exact identities**: Mathematical equivalences that hold without qualification.
2. **Formal analogies**: Same algebraic/structural form, different physical content.
3. **Heuristic connections**: Suggestive parallels that illuminate but do not prove.

### 5.2 Where the Analogy Breaks

#### Break 1: Risk-Neutral Measure $\neq$ Wick Rotation

The Girsanov change of measure $\mathbb{P} \to \mathbb{Q}$ and the Wick rotation $t \to -i\tau$ achieve similar structural results — both transform an intractable integral into a tractable one. But they are **different mathematical operations**:

- Wick rotation is an analytic continuation in the complex time plane.
- Girsanov is a reweighting of probability paths via a Radon-Nikodym derivative.

Conflating the two leads to false conclusions. For instance, Wick rotation preserves the action functional (just evaluated at imaginary time); Girsanov changes the drift of the stochastic process, altering the action itself. The Euclidean path integral obtained by Wick rotation has a different measure structure than the risk-neutral path integral obtained by Girsanov.

**Strength: formal analogy at best.**

#### Break 2: The BS Operator Is Not Self-Adjoint

In quantum mechanics, the Hamiltonian $H$ is self-adjoint (Hermitian): $H = H^\dagger$. This ensures real eigenvalues, orthogonal eigenstates, and unitary time evolution.

The Black-Scholes spatial operator $\hat{H}_{\text{BS}} = -(\sigma^2/2)\partial_{xx} - \mu\partial_x + r$ is **not self-adjoint** in the standard $L^2$ inner product. The drift term $\mu\partial_x$ is antisymmetric: $(\partial_x)^* = -\partial_x$. This means:

- Eigenvalues can be complex (though for our tridiagonal operator on a finite grid, they are real).
- The left and right eigenvectors are different.
- The evolution is contractive, not unitary.

In the QM narrative, this means the "financial Hamiltonian" generates a **dissipative** evolution, not a conservative one. The analogy between the Cayley-CN propagator and the unitary QM propagator is algebraic, not geometric.

**Strength: the self-adjoint structure is a formal analogy only.**

#### Break 3: CN Is Contractive, Not Unitary

As detailed in Chapter 3, the CN propagator $G = P^{-1}N$ satisfies $\|G\| \leq 1$ (contraction), not $\|G\| = 1$ (unitarity). The spectral radius is strictly less than 1 when $r > 0$. This means:

- Option values decay monotonically in the supremum norm (when the M-matrix holds).
- There is no conserved "probability" — the total value decreases at rate $r$ per time step.
- The Cayley form preserves contraction, not norm conservation.

**Strength: formal analogy.**

#### Break 4: Barriers Are Absorption, Not Tunneling

As discussed in Chapter 4, the barrier mechanism is **killing** (the option is destroyed), not **tunneling** (the particle penetrates). There is no wavefunction on the other side of the barrier — the value is simply set to the rebate. This is projective measurement / absorption, not potential scattering.

**Strength: formal analogy for projective measurement; heuristic at best for tunneling (which is actually incorrect).**

#### Break 5: No Interference, No Superposition

In quantum mechanics, the wave function $\psi$ is complex, and different paths through the Feynman integral interfere: $|A_1 + A_2|^2 \neq |A_1|^2 + |A_2|^2$. This interference is the hallmark of quantum behavior.

In finance, the "wave function" $V$ is real and non-negative. Path contributions to the Feynman-Kac integral add as positive reals:

$$V = \mathbb{E}^{\mathbb{Q}}[e^{-r\tau} f(X_\tau)] = \int K(x', \tau | x) f(x') dx'$$

There is no cancellation between paths, no destructive interference, no double-slit experiment for option prices. This is because the Euclidean (imaginary-time) path integral suppresses oscillation: $e^{-S_E}$ is positive and monotone, unlike $e^{iS}$ which oscillates.

**Strength: this is where the analogy breaks most completely.**

### 5.3 The Analogy Report Card

| Claim | Strength | Basis |
|-------|----------|-------|
| Gauge-transformed BS = heat equation | **Exact identity** | Change of variables + similarity transform |
| Heat equation = imaginary-time free Schrödinger | **Exact identity** | Same PDE under coefficient matching |
| $\sigma^2/2 \leftrightarrow \hbar/(2m)$, $m_{\text{eff}} = 1/\sigma^2$ | **Exact** | Coefficient identification |
| Feynman-Kac theorem | **Exact theorem** | Itô calculus proof |
| Path integral via Wiener measure | **Exact** | Rigorous measure-theoretic construction |
| $\int\mathcal{D}x\,e^{-S_E}$ notation | **Formal shorthand** | Physicist's notation for the exact construction |
| Risk-neutral measure $\leftrightarrow$ Wick rotation | **Formal analogy** | Different operations, similar structural outcomes |
| BS operator as Hamiltonian with potential | **Formal analogy** | Same PDE structure, operator not self-adjoint |
| CN = Cayley transform | **Exact algebra** | Same rational map |
| CN propagator $\leftrightarrow$ unitary QM propagator | **Formal analogy** | Same algebra, different geometric property (contraction vs. unitarity) |
| EF $\leftrightarrow$ Scharfetter-Gummel | **Exact equivalence** | Same algebraic construction for 1D cellwise convection-diffusion (up to Péclet convention) |
| MT artificial diffusion $\leftrightarrow$ numerical viscosity | **Formal analogy** | Same stabilization principle, different specific mechanisms |
| M-matrix $\leftrightarrow$ positivity preservation | **Exact** (for discrete system) | M-matrix implies non-negative solution |
| M-matrix $\leftrightarrow$ unitarity | **Incorrect** | M-matrix preserves sign, not norm |
| Barrier $\leftrightarrow$ projective measurement | **Formal analogy** | Same mathematical structure (sharp state reset), no quantum phase |
| Barrier $\leftrightarrow$ tunneling | **Incorrect** | Barriers kill; they don't transmit |
| FallbackToExponentialFitting $\leftrightarrow$ phase transition | **Heuristic** | It is adaptive stabilization, not a thermodynamic transition |
| $x\coth(x)$ in BS $\leftrightarrow$ $x\coth(x)$ in Bose-Einstein / Debye | **Formal analogy** | Same function, different physical context |

### 5.4 Finance-PDE-Physics Dictionary

| # | Finance Term | PDE Term | Physics Term | Strength |
|---|-------------|----------|-------------|----------|
| 1 | Option price $V(S,t)$ | Solution $U(x,\tau)$ | Wave function $\Psi(x,\tau)$ | **Exact** (after gauge transform) |
| 2 | Underlying price $S$ | Physical coordinate | Position | **Exact** (coordinate) |
| 3 | Log-price $x = \ln S$ | Log-coordinate | Position in natural units | **Exact** |
| 4 | Time to maturity $\tau = T-t$ | Forward time | Imaginary time | **Exact** |
| 5 | Volatility $\sigma$ | Diffusion coefficient $\sqrt{2D}$ | $\sqrt{\hbar/m_{\text{eff}}}$ | **Exact** |
| 6 | Effective mass $1/\sigma^2$ | Inverse diffusion | Particle mass $m_{\text{eff}}$ | **Exact** |
| 7 | Risk-free rate $r$ | Reaction/decay coefficient | Constant potential $V_0$ | **Formal** |
| 8 | Risk-neutral drift $\mu$ | Convection velocity | Advection bias / force field | **Formal** |
| 9 | Discount factor $e^{-r\tau}$ | Amplitude decay | Ground-state energy decay | **Formal** |
| 10 | Terminal payoff $f(S_T)$ | Initial condition | Initial wave packet | **Exact** |
| 11 | BS PDE operator | Spatial operator $A$ | Hamiltonian $\hat{H}$ | **Formal** |
| 12 | CN time step | Cayley transform | Cayley propagator | **Exact** (algebra) |
| 13 | EF fitting factor $\rho$ | Modified diffusion | Scharfetter-Gummel flux | **Exact** (1D cellwise) |
| 14 | MT added diffusion | Artificial viscosity | Numerical viscosity | **Formal** |
| 15 | M-matrix property | Non-negative inverse | Positivity preservation | **Exact** (discrete) |
| 16 | Barrier knock-out | Dirichlet reset | Absorbing boundary / projection | **Formal** |
| 17 | Monitoring date | Reset time | Measurement time | **Formal** |
| 18 | Rebate | Boundary value | Post-measurement state | **Formal** |
| 19 | Risk-neutral expectation | PDE solution representation | Path integral | **Exact** (Feynman-Kac) |
| 20 | Girsanov change of measure | Drift transformation | (Not Wick rotation) | **Formal** |
| 21 | Mesh Péclet number | Convection/diffusion ratio | Transport regime parameter | **Exact** |
| 22 | Grid convergence | Mesh refinement | Continuum limit | **Exact** |

### 5.5 Code Traceability

| Component | Source File | Paper Reference | Physics Analog |
|-----------|-----------|----------------|----------------|
| BS spatial operator assembly | `fdmblackscholesop.cpp:103–265` | — | Hamiltonian matrix construction |
| StandardCentral scheme | `fdmblackscholesop.cpp:107–152` | — | Naive tight-binding lattice |
| ExponentialFitting scheme | `fdmblackscholesop.cpp:51–57` | Duffy (2004), §3: exponentially fitted schemes | Scharfetter-Gummel flux |
| MT effective diffusion | `fdmblackscholesop.cpp:58–66` | Milev-Tagliani (2010a), Eq. 8–10; (2010b), §3 | Artificial viscosity |
| Drift correction omission | `fdmblackscholesop.cpp:59–62` (comment) | Milev-Tagliani (2010a), full 6-pt stencil | Deliberate simplification |
| $x\coth(x)$ evaluation | `fdmhyperboliccot.hpp:25–39` | Duffy (2004), fitting factor formula | Three-regime transport function |
| Spatial descriptor / scheme enum | `fdmblackscholesspatialdesc.hpp:12–67` | — | Configuration of lattice physics |
| M-matrix diagnostic | `fdmmatrixdiagnostic.hpp:52–101` | Milev-Tagliani (2010a), Theorem 3.1 conditions | Positivity check |
| Fallback to EF | `fdmblackscholesop.cpp:247–254` | — | Adaptive stabilization |
| CN time propagator | `cranknicolsonscheme.cpp:37–45` | Milev-Tagliani (2010a), §3.1; Duffy (2004), §2 | Cayley propagator |
| CN-equivalence gate | `fdmblackscholessolver.cpp:33–76` | Milev-Tagliani (2010a), Theorem 3.2 (CN requirement) | Propagator compatibility check |
| $\Delta t < 1/(rM)$ not enforced | `fdmblackscholessolver.cpp:55–58` (comment) | Milev-Tagliani (2010a), Theorem 3.2 | Time-step constraint (paper only) |
| Discrete barrier step condition | `fdmdiscretebarrierstepcondition.cpp:62–96` | Milev-Tagliani (2010a), Eq. 5 (indicator reset) | Projective measurement / absorption |
| Barrier knock-out indices | `fdmdiscretebarrierstepcondition.cpp:55–59` | — | Measurement operator support |
| Tolerant time matching | `fdmdiscretebarrierstepcondition.cpp:71–87` | — | Robust monitoring detection |
| Solver with spatial descriptor | `fdmblackscholessolver.cpp:47–89` | — | Full quantum system assembly |
| Spatial scheme tests | `fdmblackscholesspatialdiscretization.cpp` | All three papers | Lattice physics validation |
| Barrier tests | `fdmdiscretebarrier.cpp` | Milev-Tagliani (2010a), §4 examples | Measurement protocol validation |

### 5.6 Paper References

1. **Milev, M. and Tagliani, A.** "Nonstandard Finite Difference Schemes with Application to Finance: Option Pricing." *Serdica Mathematical Journal* 36:75–88, 2010.
   - Theorem 3.1: CN positivity under $\sigma^2 > r$ (S-space); adapted as Theorem 2.1 in log-space with $|\text{Pe}| < 1$.
   - Theorem 3.2: Nonstandard scheme with 6-point stencil; eigenvalue analysis showing positive distinct eigenvalues under time-step constraint.
   - Section 4: Numerical examples for discrete double barrier knock-out options ($K=100$, $L=90$, $U=110$).

2. **Milev, M. and Tagliani, A.** "Low Volatility Options and Numerical Diffusion of Finite Difference Schemes." *Serdica Mathematical Journal* 36:223–236, 2010.
   - Fitting factor formula for the nonstandard scheme.
   - Artificial diffusion comparison: $(1/8)(r\Delta S/\sigma)^2 V_{SS}$.
   - Low-volatility regime analysis where $\sigma^2 < r$.

3. **Duffy, D.J.** "A Critique of the Crank-Nicolson Scheme Strengths and Weaknesses for Financial Instrument Pricing." *Wilmott Magazine*, July 2004.
   - Exponentially fitted schemes for singular perturbation problems.
   - Scharfetter-Gummel connection for drift-diffusion equations.
   - Analysis of CN spurious oscillations near discontinuities.

---

## Epilogue: The Defensibility Test

This essay has pushed the quantum analogy aggressively — presenting Black-Scholes as a quantum system in imaginary time, the discrete operators as lattice Hamiltonians, the time stepper as a Cayley propagator, and the barrier conditions as projective measurements. But every claim has been annotated with its strength, and the "Where the Analogy Breaks" section catalogs the limits.

The test of a good analogy is not whether it is literally true, but whether it **survives removal of the metaphor**. Strip away the quantum language, and what remains?

- An exact change of variables reducing BS to the heat equation.
- An exact coefficient matching with imaginary-time Schrödinger.
- An exact theorem (Feynman-Kac) connecting PDE solutions to path integrals.
- A rigorous algebraic identity between CN and the Cayley transform.
- An exact equivalence between exponential fitting and semiconductor flux discretization.
- A precise analysis of when discrete positivity holds and when it fails.

The quantum framing illuminates these structures — it suggests why certain discretization strategies work (they respect the transport physics of the underlying equation) and why others fail (they introduce lattice artifacts inconsistent with the physical dynamics). But the results stand on their own mathematical merits.

The Black-Scholes equation is not "like" an imaginary-time Schrödinger equation. After the gauge transformation, it **is** one. The rest is a matter of how far the analogy extends into the discrete world of numerical computation — and the answer, as we have seen, is remarkably far, but not without limit.
