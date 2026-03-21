#!/usr/bin/env python3
"""plot_figures.py — Generate 9 publication-quality PDF figures from CSV data.

Reads per-scheme CSV files from results/data/ and produces PDF figures in
results/figures/.  Validates CSV metadata before plotting.

Requires: matplotlib, numpy.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Publication-quality plot configuration
# ---------------------------------------------------------------------------

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 8,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.figsize': (6.5, 4.0),
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'lines.linewidth': 1.2,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'pdf.fonttype': 42,       # TrueType fonts (no embedded timestamp)
    'ps.fonttype': 42,
})

# Colorblind-friendly palette with distinct line styles
COLORS = {
    'StandardCentral': '#0072B2',
    'ExponentialFitting': '#D55E00',
    'MilevTaglianiCN': '#009E73',
    'reference': '#000000',
    'MC': '#CC79A7',
}
STYLES = {
    'StandardCentral': '-',
    'ExponentialFitting': '--',
    'MilevTaglianiCN': '-.',
    'reference': ':',
    'MC': 'o',
}
LABELS = {
    'StandardCentral': 'Standard Central (CN)',
    'ExponentialFitting': 'Exponential Fitting (Duffy)',
    'MilevTaglianiCN': r'Milev-Tagliani CN Variant',
    'reference': 'Fine-Grid Reference',
    'MC': r'MC ($10^7$ paths)',
}

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / 'data'
FIG_DIR = SCRIPT_DIR / 'figures'

REQUIRED_META = {'scheme', 'effective_scheme', 'xGrid', 'tGrid',
                 'r', 'q', 'sigma', 'strike', 'maturity',
                 'mMatrixPolicy', 'mesh'}

SCHEME_ORDER = ['StandardCentral', 'ExponentialFitting', 'MilevTaglianiCN']

# Annotation text box style
ANNOT_BOX = dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8)


def read_csv(filename):
    """Read CSV with metadata comments.  Returns (meta, header, data_rows)."""
    path = DATA_DIR / filename
    meta = {}
    header = None
    data = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('# '):
                parts = line[2:].split(': ', 1)
                if len(parts) == 2:
                    meta[parts[0].strip()] = parts[1].strip()
            elif header is None:
                header = [c.strip() for c in line.split(',')]
            else:
                vals = line.split(',')
                if len(vals) != len(header):
                    raise ValueError(
                        f"{filename}: row has {len(vals)} cols, "
                        f"header has {len(header)}: {line!r}")
                data.append(vals)
    return meta, header, data


def validate_meta(filename, meta):
    """Validate that all required metadata fields are present."""
    missing = REQUIRED_META - set(meta.keys())
    if missing:
        raise ValueError(
            f"{filename}: missing required metadata fields: {missing}")


def to_float_arrays(header, data):
    """Convert rows to dict of numpy float arrays keyed by column name."""
    arrays = {}
    for j, col in enumerate(header):
        try:
            arrays[col] = np.array([float(row[j]) for row in data])
        except ValueError:
            arrays[col] = [row[j] for row in data]
    return arrays


def load_csv(filename, validate=True):
    """Read CSV, optionally validate metadata, return (meta, arrays)."""
    meta, header, data = read_csv(filename)
    if validate:
        validate_meta(filename, meta)
    return meta, to_float_arrays(header, data)


def save_pdf(fig, name):
    """Save figure as PDF without timestamp metadata."""
    path = FIG_DIR / name
    fig.savefig(path, metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 1: CN spurious oscillations (truncated call)
# ---------------------------------------------------------------------------

def fig1():
    sc_meta, sc = load_csv('truncated_call_StandardCentral.csv')
    _, ref = load_csv('truncated_call_reference.csv')

    fig, ax = plt.subplots()
    ax.plot(sc['S'], sc['price'], color=COLORS['StandardCentral'],
            ls=STYLES['StandardCentral'], label=LABELS['StandardCentral'])
    ax.plot(ref['S'], ref['price'], color=COLORS['reference'],
            ls=STYLES['reference'], lw=0.8, label=LABELS['reference'])

    K, U = float(sc_meta['strike']), float(sc_meta['upper_barrier'])
    mask = (sc['S'] > K) & (sc['S'] < U + 5)
    ymax = max(sc['price'][mask]) * 1.3 if mask.any() else 25
    ax.set_xlim(K - 5, U + 10)
    ax.set_ylim(-5, ymax)
    ax.set_xlabel(r'Stock Price $S$')
    ax.set_ylabel(r'Option Value $V(S)$')
    ax.set_title(r'Spurious Oscillations in Crank-Nicolson '
                 r'($\sigma=0.001$, $r=0.05$)')
    ax.legend(loc='upper left')
    ax.axvline(U, color='gray', ls=':', lw=0.5, alpha=0.5)

    # Annotation: explain physical meaning of oscillations
    ax.text(0.98, 0.55,
            'Negative prices near $U$:\n'
            r'M-matrix violated ($\sigma^2 \ll r$).'
            '\nEigenvalues of $P^{-1}N \\to -1$.',
            transform=ax.transAxes, fontsize=7,
            verticalalignment='top', horizontalalignment='right',
            bbox=ANNOT_BOX)

    save_pdf(fig, 'fig1_cn_oscillations.pdf')
    print('  Fig 1: CN oscillations')


# ---------------------------------------------------------------------------
# Figure 2: Three-scheme comparison (truncated call)
# ---------------------------------------------------------------------------

def fig2():
    schemes = {}
    for s in SCHEME_ORDER:
        _, schemes[s] = load_csv(f'truncated_call_{s}.csv')

    ref_meta, ref = load_csv('truncated_call_reference.csv')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4),
                                     gridspec_kw={'width_ratios': [2, 1]})

    K = float(ref_meta['strike'])
    U = float(ref_meta['upper_barrier'])

    for s in SCHEME_ORDER:
        ax1.plot(schemes[s]['S'], schemes[s]['price'],
                 color=COLORS[s], ls=STYLES[s], label=LABELS[s])
    ax1.plot(ref['S'], ref['price'], color=COLORS['reference'],
             ls=STYLES['reference'], lw=0.8, label=LABELS['reference'])

    ax1.set_xlim(K - 5, U + 10)
    ax1.set_ylim(-5, 25)
    ax1.set_xlabel(r'Stock Price $S$')
    ax1.set_ylabel(r'Option Value $V(S)$')
    ax1.set_title('Three-Scheme Comparison: Truncated Call')
    ax1.legend(loc='upper left', fontsize=7)
    ax1.axvline(U, color='gray', ls=':', lw=0.5, alpha=0.5)

    # Zoom near U
    for s in SCHEME_ORDER:
        mask = (schemes[s]['S'] > U - 3) & (schemes[s]['S'] < U + 3)
        ax2.plot(schemes[s]['S'][mask], schemes[s]['price'][mask],
                 color=COLORS[s], ls=STYLES[s])
    mask_r = (ref['S'] > U - 3) & (ref['S'] < U + 3)
    ax2.plot(ref['S'][mask_r], ref['price'][mask_r],
             color=COLORS['reference'], ls=STYLES['reference'], lw=0.8)
    ax2.set_xlabel(r'$S$')
    ax2.set_ylabel(r'$V(S)$')
    ax2.set_title(f'Zoom near $U={int(U)}$')
    ax2.axvline(U, color='gray', ls=':', lw=0.5, alpha=0.5)

    # Annotation on zoom panel
    ax2.text(0.95, 0.95, 'EF and MT:\nsmooth, positive',
             transform=ax2.transAxes, fontsize=7,
             verticalalignment='top', horizontalalignment='right',
             bbox=ANNOT_BOX)

    fig.tight_layout()
    save_pdf(fig, 'fig2_three_scheme_truncated.pdf')
    print('  Fig 2: Three-scheme truncated call')


# ---------------------------------------------------------------------------
# Figure 3: Moderate-vol discrete barrier
# ---------------------------------------------------------------------------

def fig3():
    fig, ax = plt.subplots()
    first_meta = None
    for s in SCHEME_ORDER:
        m, arr = load_csv(f'barrier_moderate_vol_{s}.csv')
        if first_meta is None:
            first_meta = m
            L = float(m['lower_barrier'])
            U = float(m['upper_barrier'])
        mask = (arr['S'] >= L - 5) & (arr['S'] <= U + 5)
        ax.plot(arr['S'][mask], arr['price'][mask],
                color=COLORS[s], ls=STYLES[s], label=LABELS[s])

    # Fine-grid reference — thicker, distinct style for visibility
    try:
        _, ref = load_csv(
            'barrier_moderate_vol_reference_ExponentialFitting.csv')
        mask_r = (ref['S'] >= L - 5) & (ref['S'] <= U + 5)
        ax.plot(ref['S'][mask_r], ref['price'][mask_r],
                color=COLORS['reference'], ls=(0, (3, 1, 1, 1)), lw=2.0,
                label='Fine-Grid Reference (16k)', zorder=0)
    except FileNotFoundError:
        pass

    # MC reference
    try:
        mc_meta, mc = load_csv('mc_barrier_moderate_vol.csv', validate=False)
        n_paths = mc_meta.get('num_paths', '10000000')
        mc_label = f'MC ({int(n_paths)/1e6:.0f}M paths, 95% CI)'
        mc_mask = (mc['S0'] >= L - 5) & (mc['S0'] <= U + 5)
        ax.errorbar(mc['S0'][mc_mask], mc['price'][mc_mask],
                     yerr=2*mc['standard_error'][mc_mask],
                     fmt='o', color=COLORS['MC'], ms=2.5, lw=0.6,
                     capsize=1.5, label=mc_label)
    except FileNotFoundError:
        pass

    ax.set_xlim(L - 5, U + 5)
    ax.set_xlabel(r'Stock Price $S$')
    ax.set_ylabel(r'Option Value $V(S)$')
    ax.set_title(r'Discrete Double Barrier '
                 r'($\sigma=0.25$, $K=100$, $L=95$, $U=110$)')
    ax.legend(fontsize=6.5, loc='upper left')
    ax.axvline(L, color='gray', ls='--', lw=0.8, alpha=0.6)
    ax.axvline(U, color='gray', ls='--', lw=0.8, alpha=0.6)
    ax.text(L, ax.get_ylim()[1]*0.95, f' $L={int(L)}$', fontsize=7,
            va='top', color='gray')
    ax.text(U, ax.get_ylim()[1]*0.95, f' $U={int(U)}$', fontsize=7,
            va='top', color='gray')

    # Annotation: MC metadata
    ax.text(0.98, 0.45, f'10M paths, 95% CI\n5 monitoring dates\n$T=0.5$',
            transform=ax.transAxes, fontsize=7,
            verticalalignment='top', horizontalalignment='right',
            bbox=ANNOT_BOX)

    save_pdf(fig, 'fig3_barrier_moderate.pdf')
    print('  Fig 3: Moderate-vol barrier')


# ---------------------------------------------------------------------------
# Figure 4: Low-vol discrete barrier
# ---------------------------------------------------------------------------

def fig4():
    fig, ax = plt.subplots()
    first_meta = None
    for s in SCHEME_ORDER:
        m, arr = load_csv(f'barrier_low_vol_{s}.csv')
        if first_meta is None:
            first_meta = m
            L = float(m['lower_barrier'])
            U = float(m['upper_barrier'])
        mask = (arr['S'] >= L - 5) & (arr['S'] <= U + 5)
        ax.plot(arr['S'][mask], arr['price'][mask],
                color=COLORS[s], ls=STYLES[s], label=LABELS[s])

    # MC reference (low-vol)
    try:
        mc_meta, mc = load_csv('mc_barrier_low_vol.csv', validate=False)
        n_paths = mc_meta.get('num_paths', '5000000')
        mc_label = f'MC ({int(n_paths)/1e6:.0f}M paths, 95% CI)'
        mc_mask = (mc['S0'] >= L - 5) & (mc['S0'] <= U + 5)
        ax.errorbar(mc['S0'][mc_mask], mc['price'][mc_mask],
                     yerr=2*mc['standard_error'][mc_mask],
                     fmt='o', color=COLORS['MC'], ms=2.5, lw=0.6,
                     capsize=1.5, label=mc_label)
    except FileNotFoundError:
        pass

    ax.axhline(0, color='red', ls=':', lw=0.5, alpha=0.7)
    ax.set_xlim(L - 5, U + 5)
    ax.set_xlabel(r'Stock Price $S$')
    ax.set_ylabel(r'Option Value $V(S)$')
    ax.set_title(r'Low-Vol Discrete Barrier '
                 r'($\sigma=0.001$, $\sigma^2 \ll r$)')
    ax.legend(fontsize=7)
    ax.axvline(L, color='gray', ls='--', lw=0.8, alpha=0.6)
    ax.axvline(U, color='gray', ls='--', lw=0.8, alpha=0.6)
    ax.text(L, ax.get_ylim()[1]*0.95, f' $L={int(L)}$', fontsize=7,
            va='top', color='gray')
    ax.text(U, ax.get_ylim()[1]*0.95, f' $U={int(U)}$', fontsize=7,
            va='top', color='gray')

    # Shade the V < 0 region to highlight SC's non-physical negative prices
    ylim = ax.get_ylim()
    if ylim[0] < 0:
        ax.axhspan(ylim[0], 0, color='red', alpha=0.08, zorder=0)
        ax.text(0.02, 0.02,
                'Shaded: SC negative prices\n(M-matrix violation)',
                transform=ax.transAxes, fontsize=6.5,
                verticalalignment='bottom', horizontalalignment='left',
                color='red', bbox=dict(boxstyle='round,pad=0.2',
                                       facecolor='mistyrose', alpha=0.8))

    # Annotation: MT diffusion limitation
    ax.text(0.98, 0.55,
            r'MT undershoot at $\sigma^2 \ll r$:' '\n'
            r'artificial diffusion $\frac{1}{8}(\frac{r}{\sigma}\Delta S)^2$'
            '\nmay exceed true diffusion.',
            transform=ax.transAxes, fontsize=6.5,
            verticalalignment='top', horizontalalignment='right',
            bbox=ANNOT_BOX)

    save_pdf(fig, 'fig4_barrier_lowvol.pdf')
    print('  Fig 4: Low-vol barrier')


# ---------------------------------------------------------------------------
# Figure 5: Grid convergence
# ---------------------------------------------------------------------------

def fig5():
    fig, ax = plt.subplots()
    scheme_data = {}
    for s in SCHEME_ORDER:
        _, arr = load_csv(f'grid_convergence_{s}.csv')
        scheme_data[s] = arr
        ax.loglog(arr['xGrid'], arr['error'],
                  color=COLORS[s], ls=STYLES[s], marker='s', ms=4,
                  label=LABELS[s])

    # O(h^2) reference slope (use last scheme; all share the same xGrid)
    last = scheme_data[SCHEME_ORDER[-1]]
    x = last['xGrid']
    y0 = last['error'][0] if last['error'][0] > 0 else 1e-2
    ax.loglog(x, y0 * (x[0]/x)**2, 'k:', lw=0.5, alpha=0.5,
              label=r'$O(h^2)$ slope')

    # Convergence rate annotations via least-squares fit on last 4 points
    for i, s in enumerate(SCHEME_ORDER):
        d = scheme_data[s]
        log_x = np.log(d['xGrid'][-4:])
        log_e = np.log(d['error'][-4:])
        slope, _ = np.polyfit(log_x, log_e, 1)
        rate = -slope
        # Place annotation near the last point of each curve
        last_x = d['xGrid'][-1]
        last_y = d['error'][-1]
        ax.annotate(
            f'$\\alpha = {rate:.2f}$',
            xy=(last_x, last_y),
            xytext=(15, 5 + i*12),
            textcoords='offset points',
            fontsize=7, color=COLORS[s],
            arrowprops=dict(arrowstyle='->', color=COLORS[s], lw=0.5))

    # Label annotation explaining refinement type
    ax.text(0.02, 0.02,
            r'$\alpha$: joint (spatial+temporal) rate, last 4 points',
            transform=ax.transAxes, fontsize=6.5,
            verticalalignment='bottom', horizontalalignment='left',
            bbox=ANNOT_BOX)

    ax.set_xlabel(r'Spatial Grid Size ($N_x$)')
    ax.set_ylabel(r'$|V_{\mathrm{num}} - V_{\mathrm{BS}}|$')
    ax.set_title(r'Grid Convergence: European Call ($\sigma=0.20$)')
    ax.legend(fontsize=7)
    save_pdf(fig, 'fig5_convergence.pdf')
    print('  Fig 5: Grid convergence')


# ---------------------------------------------------------------------------
# Figure 6: Effective diffusion comparison
# ---------------------------------------------------------------------------

def fig6():
    fig, ax = plt.subplots()
    sweep_data = {}
    sweep_meta = {}
    for s in SCHEME_ORDER:
        m, arr = load_csv(f'effective_diffusion_sweep_{s}.csv')
        sweep_meta[s] = m
        sweep_data[s] = arr
        ax.loglog(arr['sigma'], arr['a_eff'], color=COLORS[s], ls=STYLES[s],
                  label=LABELS[s])

    ax.set_xlabel(r'Volatility $\sigma$')
    ax.set_ylabel(r'Effective Diffusion $a_{\mathrm{eff}}$')
    ax.set_title(r'Effective Diffusion vs Volatility ($\sigma$-Sweep)')
    ax.legend()

    # Regime transition marker: σ_* where Pe(σ) = 1 on this mesh.
    # Pe = μh/σ² where μ = r-q-σ²/2.  Solve numerically from sweep metadata.
    sc_meta = sweep_meta['StandardCentral']
    h_sweep = float(sc_meta.get('h', '0.006966'))
    r_sweep = float(sc_meta.get('r', '0.05'))
    q_sweep = float(sc_meta.get('q', '0.0'))
    sc_sigma = sweep_data['StandardCentral']['sigma']
    pe_vals = np.array([(r_sweep - q_sweep - s**2/2) * h_sweep / s**2
                        for s in sc_sigma])
    # Interpolate to find σ where Pe crosses 1
    cross_idx = np.where(np.diff(np.sign(pe_vals - 1.0)))[0]
    if len(cross_idx) > 0:
        i = cross_idx[0]
        # Linear interpolation between samples
        f = (1.0 - pe_vals[i]) / (pe_vals[i+1] - pe_vals[i])
        sigma_star = sc_sigma[i] + f * (sc_sigma[i+1] - sc_sigma[i])
        ax.axvline(sigma_star, color='gray', ls=':', lw=0.8, alpha=0.7)
        ax.text(sigma_star * 1.5, ax.get_ylim()[0] * 3,
                f'$\\sigma_* \\approx {sigma_star:.4f}$\n'
                '$Pe(\\sigma_*) = 1$',
                fontsize=6.5, color='#555555', bbox=ANNOT_BOX)

    # Annotation: regime summary
    ax.text(0.03, 0.97,
            r'Low $\sigma$: MT $\gg$ EF $\gg$ SC' '\n'
            r'High $\sigma$: all converge to $\sigma^2/2$',
            transform=ax.transAxes, fontsize=7,
            verticalalignment='top', horizontalalignment='left',
            bbox=ANNOT_BOX)

    save_pdf(fig, 'fig6_effective_diffusion.pdf')
    print('  Fig 6: Effective diffusion')


# ---------------------------------------------------------------------------
# Figure 7: M-matrix off-diagonal visualization
# ---------------------------------------------------------------------------

def fig7():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    sc_data = None
    for s in SCHEME_ORDER:
        _, arr = load_csv(f'mmatrix_sweep_{s}.csv')
        if s == 'StandardCentral':
            sc_data = arr
        ax1.semilogx(arr['sigma'], arr['lower'], color=COLORS[s],
                     ls=STYLES[s], label=LABELS[s])
        ax2.loglog(arr['sigma'], arr['upper'], color=COLORS[s],
                   ls=STYLES[s], label=LABELS[s])

    ax1.axhline(0, color='red', ls='--', lw=1.0, alpha=0.8)
    ax1.set_yscale('symlog', linthresh=1e-2)
    ax1.set_xlabel(r'Volatility $\sigma$')
    ax1.set_ylabel(r'Lower Off-Diagonal $a_{i,i-1}$')
    ax1.set_title('Lower Off-Diagonal')
    ax1.legend(fontsize=7)

    # Find and annotate critical sigma where SC crosses zero
    if sc_data is not None:
        lower = sc_data['lower']
        sigma = sc_data['sigma']
        sign_changes = np.where(np.diff(np.sign(lower)))[0]
        if len(sign_changes) > 0:
            idx = sign_changes[0]
            sig_crit = sigma[idx]
            ax1.axvline(sig_crit, color='red', ls=':', lw=0.6, alpha=0.6)
            ax1.text(sig_crit * 1.3, ax1.get_ylim()[1] * 0.3,
                     f'SC: M-matrix\nviolated below\n'
                     f'$\\sigma \\approx {sig_crit:.3f}$',
                     fontsize=6.5, color='red', bbox=ANNOT_BOX)

    ax2.set_xlabel(r'Volatility $\sigma$')
    ax2.set_ylabel(r'Upper Off-Diagonal $a_{i,i+1}$')
    ax2.set_title('Upper Off-Diagonal')
    ax2.legend(fontsize=7)

    fig.suptitle(r'M-Matrix Off-Diagonals vs Volatility ($\sigma$-Sweep)',
                 y=1.02)
    fig.tight_layout()
    save_pdf(fig, 'fig7_mmatrix.pdf')
    print('  Fig 7: M-matrix off-diagonals')


# ---------------------------------------------------------------------------
# Figure 8: Performance benchmark
# ---------------------------------------------------------------------------

def fig8():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    SHORT_NAMES = {'StandardCentral': 'SC', 'ExponentialFitting': 'EF',
                   'MilevTaglianiCN': 'MT'}

    scheme_meta = {}
    scheme_data = {}
    for s in SCHEME_ORDER:
        m, arr = load_csv(f'benchmark_{s}.csv')
        scheme_meta[s] = m
        scheme_data[s] = arr
        ax1.loglog(arr['xGrid'], arr['relative_cost'],
                   color=COLORS[s], ls=STYLES[s], marker='s', ms=4,
                   label=LABELS[s])
        ax2.loglog(arr['relative_cost'], arr['error'],
                   color=COLORS[s], ls=STYLES[s], marker='s', ms=4,
                   label=LABELS[s])

    ax1.set_xlabel(r'Grid Size ($N_x$)')
    ax1.set_ylabel(r'Work Units ($N_x \times N_t$)')
    ax1.set_title('Cost vs Grid Size')
    ax1.legend(fontsize=7)

    ax2.set_xlabel(r'Work Units ($N_x \times N_t$)')
    ax2.set_ylabel(r'$|V_{\mathrm{num}} - V_{\mathrm{BS}}|$')
    ax2.set_title('Cost-Accuracy Tradeoff')
    ax2.legend(fontsize=7)

    # Data-grounded accuracy annotation from each scheme at N_x=200
    lines = []
    for s in SCHEME_ORDER:
        s_arr = scheme_data[s]
        if 'wall_ms' in s_arr:
            idx_200 = np.argmin(np.abs(s_arr['xGrid'] - 200))
            err_200 = s_arr['error'][idx_200]
            ref_price = float(scheme_meta[s].get('reference_price', '9.227'))
            rel_err = err_200 / ref_price * 100
            wall_200 = s_arr['wall_ms'][idx_200]
            lines.append(f'{SHORT_NAMES[s]}: {rel_err:.3f}%, {wall_200:.1f} ms')
    if lines:
        ax2.text(0.98, 0.95,
                 f'At $N_x=200$:\n' + '\n'.join(lines),
                 transform=ax2.transAxes, fontsize=6.5,
                 verticalalignment='top', horizontalalignment='right',
                 bbox=ANNOT_BOX)

    fig.suptitle('Performance Benchmark: European Call', y=1.02)
    fig.tight_layout()
    save_pdf(fig, 'fig8_benchmark.pdf')
    print('  Fig 8: Benchmark')


# ---------------------------------------------------------------------------
# Figure 9: xCothx / Peclet number
# ---------------------------------------------------------------------------

def fig9():
    _, arr = load_csv('xcothx.csv', validate=False)

    fig, ax = plt.subplots()
    ax.plot(arr['Pe'], arr['xCothx'], color=COLORS['ExponentialFitting'],
            ls='-', lw=1.5, label=r'$x \coth(x) = \rho(Pe)$')
    ax.plot(arr['Pe'], arr['abs_Pe'], color='gray', ls=':', lw=0.8,
            label=r'$|Pe|$ (upwind limit)')
    ax.axhline(1.0, color='gray', ls='--', lw=0.5, alpha=0.5,
               label=r'$\rho=1$ (central limit)')

    ax.set_xlabel('P\u00e9clet Number $Pe$')
    ax.set_ylabel(r'Fitting Factor $\rho = x\coth(x)$')
    ax.set_title('Duffy Fitting Factor: Implementation Regimes')
    ax.legend()
    ax.set_xlim(-80, 80)
    ax.set_ylim(0, 82)

    # Implementation regime boundaries: |x| < 1e-6 (Taylor), 1e-6 <= |x| <= 50 (Direct), |x| > 50 (Asymptotic)
    ax.axvline(50, color='purple', ls='--', lw=0.6, alpha=0.5)
    ax.axvline(-50, color='purple', ls='--', lw=0.6, alpha=0.5)

    # Regime labels placed in their actual implementation ranges
    ax.text(0, 5, 'Taylor\n$\\rho \\approx 1 + x^2/3$\n($|x| < 10^{-6}$)',
            fontsize=7.5, ha='center', color='#555555',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow',
                      alpha=0.7))
    ax.text(25, 35, 'Direct\n$\\rho = x\\coth(x)$\n($10^{-6} \\leq |x| \\leq 50$)',
            fontsize=7.5, ha='center', color='#555555',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow',
                      alpha=0.7))
    ax.text(65, 68, 'Asymptotic\n$\\rho \\approx |x|$\n($|x| > 50$)',
            fontsize=7.5, ha='center', color='#555555',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow',
                      alpha=0.7))

    save_pdf(fig, 'fig9_xcothx.pdf')
    print('  Fig 9: xCothx/Peclet')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    if not DATA_DIR.exists():
        print(f"Error: Data directory {DATA_DIR} not found. "
              "Run generate_data first.", file=sys.stderr)
        sys.exit(1)

    print("Generating figures...")
    fig9()
    fig7()
    fig6()
    fig5()
    fig1()
    fig2()
    fig3()
    fig4()
    fig8()
    print(f"All 9 figures saved to {FIG_DIR}/")


if __name__ == '__main__':
    main()
