#!/usr/bin/env python3
"""plot_figures.py — Generate 9 publication-quality PDF figures from CSV data.

Reads per-scheme CSV files from results/data/ and produces PDF figures in
results/figures/.  Validates CSV metadata before plotting.

Requires: matplotlib, numpy.
"""

import csv
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
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
    'MC': r'MC ($5 \times 10^6$ paths)',
}

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / 'data'
FIG_DIR = SCRIPT_DIR / 'figures'

REQUIRED_META = {'scheme', 'effective_scheme', 'xGrid', 'tGrid',
                 'r', 'q', 'sigma', 'strike', 'maturity',
                 'mMatrixPolicy', 'mesh'}

SCHEME_ORDER = ['StandardCentral', 'ExponentialFitting', 'MilevTaglianiCN']


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


def save_pdf(fig, name):
    """Save figure as PDF without timestamp metadata."""
    path = FIG_DIR / name
    fig.savefig(path, metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 1: CN spurious oscillations (truncated call)
# ---------------------------------------------------------------------------

def fig1():
    sc_meta, sc_h, sc_d = read_csv('truncated_call_StandardCentral.csv')
    validate_meta('truncated_call_StandardCentral.csv', sc_meta)
    sc = to_float_arrays(sc_h, sc_d)

    ref_meta, ref_h, ref_d = read_csv('truncated_call_reference.csv')
    validate_meta('truncated_call_reference.csv', ref_meta)
    ref = to_float_arrays(ref_h, ref_d)

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
    save_pdf(fig, 'fig1_cn_oscillations.pdf')
    print('  Fig 1: CN oscillations')


# ---------------------------------------------------------------------------
# Figure 2: Three-scheme comparison (truncated call)
# ---------------------------------------------------------------------------

def fig2():
    schemes = {}
    for s in SCHEME_ORDER:
        m, h, d = read_csv(f'truncated_call_{s}.csv')
        validate_meta(f'truncated_call_{s}.csv', m)
        schemes[s] = to_float_arrays(h, d)

    ref_meta, ref_h, ref_d = read_csv('truncated_call_reference.csv')
    ref = to_float_arrays(ref_h, ref_d)

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
        m, h, d = read_csv(f'barrier_moderate_vol_{s}.csv')
        validate_meta(f'barrier_moderate_vol_{s}.csv', m)
        arr = to_float_arrays(h, d)
        if first_meta is None:
            first_meta = m
        L = float(m['lower_barrier'])
        U = float(m['upper_barrier'])
        mask = (arr['S'] >= L - 1) & (arr['S'] <= U + 1)
        ax.plot(arr['S'][mask], arr['price'][mask],
                color=COLORS[s], ls=STYLES[s], label=LABELS[s])

    # Fine-grid reference
    try:
        ref_m, ref_h, ref_d = read_csv(
            'barrier_moderate_vol_reference_ExponentialFitting.csv')
        ref = to_float_arrays(ref_h, ref_d)
        mask_r = (ref['S'] >= L - 1) & (ref['S'] <= U + 1)
        ax.plot(ref['S'][mask_r], ref['price'][mask_r],
                color=COLORS['reference'], ls=STYLES['reference'], lw=0.8,
                label='Fine-Grid Reference (16k)')
    except FileNotFoundError:
        pass

    # MC reference
    try:
        mc_m, mc_h, mc_d = read_csv('mc_barrier_moderate_vol.csv')
        mc = to_float_arrays(mc_h, mc_d)
        ax.errorbar(mc['S0'], mc['price'], yerr=2*mc['standard_error'],
                     fmt='o', color=COLORS['MC'], ms=3, lw=0.8,
                     capsize=2, label=LABELS['MC'])
    except FileNotFoundError:
        pass

    L = float(first_meta['lower_barrier'])
    U = float(first_meta['upper_barrier'])
    ax.set_xlabel(r'Stock Price $S$')
    ax.set_ylabel(r'Option Value $V(S)$')
    ax.set_title(r'Discrete Double Barrier '
                 r'($\sigma=0.25$, $K=100$, $L=95$, $U=110$)')
    ax.legend(fontsize=7)
    ax.axvline(L, color='gray', ls=':', lw=0.5, alpha=0.5)
    ax.axvline(U, color='gray', ls=':', lw=0.5, alpha=0.5)
    save_pdf(fig, 'fig3_barrier_moderate.pdf')
    print('  Fig 3: Moderate-vol barrier')


# ---------------------------------------------------------------------------
# Figure 4: Low-vol discrete barrier
# ---------------------------------------------------------------------------

def fig4():
    fig, ax = plt.subplots()
    first_meta = None
    for s in SCHEME_ORDER:
        m, h, d = read_csv(f'barrier_low_vol_{s}.csv')
        validate_meta(f'barrier_low_vol_{s}.csv', m)
        arr = to_float_arrays(h, d)
        if first_meta is None:
            first_meta = m
        L = float(m['lower_barrier'])
        U = float(m['upper_barrier'])
        mask = (arr['S'] >= L - 1) & (arr['S'] <= U + 1)
        ax.plot(arr['S'][mask], arr['price'][mask],
                color=COLORS[s], ls=STYLES[s], label=LABELS[s])

    # MC reference (low-vol)
    try:
        mc_m, mc_h, mc_d = read_csv('mc_barrier_low_vol.csv')
        mc = to_float_arrays(mc_h, mc_d)
        ax.errorbar(mc['S0'], mc['price'], yerr=2*mc['standard_error'],
                     fmt='o', color=COLORS['MC'], ms=3, lw=0.8,
                     capsize=2, label=LABELS['MC'])
    except FileNotFoundError:
        pass

    ax.axhline(0, color='red', ls=':', lw=0.5, alpha=0.7)
    ax.set_xlabel(r'Stock Price $S$')
    ax.set_ylabel(r'Option Value $V(S)$')
    ax.set_title(r'Low-Vol Discrete Barrier '
                 r'($\sigma=0.001$, $\sigma^2 \ll r$)')
    ax.legend()
    L = float(first_meta['lower_barrier'])
    U = float(first_meta['upper_barrier'])
    ax.axvline(L, color='gray', ls=':', lw=0.5, alpha=0.5)
    ax.axvline(U, color='gray', ls=':', lw=0.5, alpha=0.5)
    save_pdf(fig, 'fig4_barrier_lowvol.pdf')
    print('  Fig 4: Low-vol barrier')


# ---------------------------------------------------------------------------
# Figure 5: Grid convergence
# ---------------------------------------------------------------------------

def fig5():
    fig, ax = plt.subplots()
    for s in SCHEME_ORDER:
        m, h, d = read_csv(f'grid_convergence_{s}.csv')
        validate_meta(f'grid_convergence_{s}.csv', m)
        arr = to_float_arrays(h, d)
        ax.loglog(arr['xGrid'], arr['error'],
                  color=COLORS[s], ls=STYLES[s], marker='s', ms=4,
                  label=LABELS[s])

    # O(h^2) reference slope
    x = arr['xGrid']
    y0 = arr['error'][0] if arr['error'][0] > 0 else 1e-2
    ax.loglog(x, y0 * (x[0]/x)**2, 'k:', lw=0.5, alpha=0.5,
              label=r'$O(h^2)$ slope')

    ax.set_xlabel(r'Spatial Grid Size ($N_x$)')
    ax.set_ylabel(r'$|V_{\mathrm{num}} - V_{\mathrm{BS}}|$')
    ax.set_title(r'Grid Convergence: European Call ($\sigma=0.20$)')
    ax.legend()
    save_pdf(fig, 'fig5_convergence.pdf')
    print('  Fig 5: Grid convergence')


# ---------------------------------------------------------------------------
# Figure 6: Effective diffusion comparison
# ---------------------------------------------------------------------------

def fig6():
    fig, ax = plt.subplots()
    for s in SCHEME_ORDER:
        m, h, d = read_csv(f'effective_diffusion_{s}.csv')
        validate_meta(f'effective_diffusion_{s}.csv', m)
        arr = to_float_arrays(h, d)
        ax.plot(arr['S'], arr['aUsed'], color=COLORS[s], ls=STYLES[s],
                label=LABELS[s])

    ax.set_xlabel(r'Stock Price $S$')
    ax.set_ylabel(r'Effective Diffusion $a_{\mathrm{eff}}$')
    ax.set_title(r'Artificial Diffusion: Duffy vs MT ($\sigma=0.001$)')
    ax.legend()
    ax.set_yscale('log')
    save_pdf(fig, 'fig6_effective_diffusion.pdf')
    print('  Fig 6: Effective diffusion')


# ---------------------------------------------------------------------------
# Figure 7: M-matrix off-diagonal visualization
# ---------------------------------------------------------------------------

def fig7():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    for s in SCHEME_ORDER:
        m, h, d = read_csv(f'mmatrix_offdiag_{s}.csv')
        validate_meta(f'mmatrix_offdiag_{s}.csv', m)
        arr = to_float_arrays(h, d)
        ax1.plot(arr['S'], arr['lower'], color=COLORS[s], ls=STYLES[s],
                 label=LABELS[s])
        ax2.plot(arr['S'], arr['upper'], color=COLORS[s], ls=STYLES[s],
                 label=LABELS[s])

    for ax, title in [(ax1, r'Lower Off-Diagonal $a_{i,i-1}$'),
                      (ax2, r'Upper Off-Diagonal $a_{i,i+1}$')]:
        ax.axhline(0, color='red', ls=':', lw=0.8, alpha=0.7)
        ax.set_xlabel(r'Stock Price $S$')
        ax.set_ylabel(title)
        ax.set_title(title.split('$')[0].strip())
        ax.legend(fontsize=7)

    fig.suptitle(r'M-Matrix Property: Off-Diagonal Signs ($\sigma=0.001$)',
                 y=1.02)
    fig.tight_layout()
    save_pdf(fig, 'fig7_mmatrix.pdf')
    print('  Fig 7: M-matrix off-diagonals')


# ---------------------------------------------------------------------------
# Figure 8: Performance benchmark
# ---------------------------------------------------------------------------

def fig8():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    for s in SCHEME_ORDER:
        m, h, d = read_csv(f'benchmark_{s}.csv')
        validate_meta(f'benchmark_{s}.csv', m)
        arr = to_float_arrays(h, d)
        ax1.loglog(arr['xGrid'], arr['relative_cost'],
                   color=COLORS[s], ls=STYLES[s], marker='s', ms=4,
                   label=LABELS[s])
        ax2.loglog(arr['relative_cost'], arr['error'],
                   color=COLORS[s], ls=STYLES[s], marker='s', ms=4,
                   label=LABELS[s])

    ax1.set_xlabel(r'Grid Size ($N_x$)')
    ax1.set_ylabel(r'Relative Cost ($N_x \times N_t$)')
    ax1.set_title('Cost vs Grid Size')
    ax1.legend(fontsize=7)

    ax2.set_xlabel(r'Relative Cost ($N_x \times N_t$)')
    ax2.set_ylabel(r'$|V_{\mathrm{num}} - V_{\mathrm{BS}}|$')
    ax2.set_title('Cost-Accuracy Tradeoff')
    ax2.legend(fontsize=7)

    fig.suptitle('Performance Benchmark: European Call', y=1.02)
    fig.tight_layout()
    save_pdf(fig, 'fig8_benchmark.pdf')
    print('  Fig 8: Benchmark')


# ---------------------------------------------------------------------------
# Figure 9: xCothx / Peclet number
# ---------------------------------------------------------------------------

def fig9():
    m, h, d = read_csv('xcothx.csv')
    # xcothx has N/A for solver fields — skip strict validation
    arr = to_float_arrays(h, d)

    fig, ax = plt.subplots()
    ax.plot(arr['Pe'], arr['xCothx'], color=COLORS['ExponentialFitting'],
            ls='-', lw=1.5, label=r'$x \coth(x) = \rho(Pe)$')
    ax.plot(arr['Pe'], arr['abs_Pe'], color='gray', ls=':', lw=0.8,
            label=r'$|Pe|$ (upwind limit)')
    ax.axhline(1.0, color='gray', ls='--', lw=0.5, alpha=0.5,
               label=r'$\rho=1$ (central limit)')

    ax.set_xlabel(r'P\'eclet Number $Pe$')
    ax.set_ylabel(r'Fitting Factor $\rho = x\coth(x)$')
    ax.set_title('Duffy Fitting Factor: Regime Transition')
    ax.legend()
    ax.set_xlim(-20, 20)
    ax.set_ylim(0, 22)
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
