#!/usr/bin/env python3
"""plot_figures.py — Generate 9 publication-quality PDF figures from CSV data.

Reads CSV files from results/data/ and produces PDF figures in results/figures/.
Requires: matplotlib, numpy.
"""

import csv
import os
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
})

# Use LaTeX-style fonts if available
try:
    plt.rcParams.update({
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath}\usepackage{amssymb}',
    })
except Exception:
    pass

# Colorblind-friendly palette with distinct line styles
COLORS = {
    'StandardCentral': '#0072B2',       # Blue
    'ExponentialFitting': '#D55E00',    # Vermillion
    'MilevTaglianiCN': '#009E73',       # Green
    'Analytical': '#000000',            # Black
    'MC': '#CC79A7',                    # Pink
}
STYLES = {
    'StandardCentral': '-',
    'ExponentialFitting': '--',
    'MilevTaglianiCN': '-.',
    'Analytical': ':',
    'MC': 'o',
}
LABELS = {
    'StandardCentral': 'Standard Central (CN)',
    'ExponentialFitting': 'Exponential Fitting (Duffy)',
    'MilevTaglianiCN': r'Milev-Tagliani CN Variant',
    'Analytical': 'Analytical (Black-Scholes)',
    'MC': 'Monte Carlo Reference',
}

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / 'data'
FIG_DIR = SCRIPT_DIR / 'figures'


def read_csv(filename):
    """Read CSV with metadata comments. Returns (metadata_dict, header, data_rows)."""
    meta = {}
    header = None
    data = []
    with open(DATA_DIR / filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith('# '):
                parts = line[2:].split(': ', 1)
                if len(parts) == 2:
                    meta[parts[0]] = parts[1]
            elif header is None:
                header = line.split(',')
            else:
                data.append([float(x) if x.replace('.','',1).replace('-','',1).replace('e','',1).replace('+','',1).isdigit() else x
                             for x in line.split(',')])
    return meta, header, data


def to_arrays(header, data):
    """Convert list-of-lists to dict of numpy arrays keyed by column name."""
    arrays = {}
    for j, col in enumerate(header):
        try:
            arrays[col] = np.array([row[j] for row in data], dtype=float)
        except (ValueError, TypeError):
            arrays[col] = [row[j] for row in data]
    return arrays


# ---------------------------------------------------------------------------
# Figure 1: CN spurious oscillations (truncated call)
# ---------------------------------------------------------------------------

def fig1_cn_oscillations():
    meta, header, data = read_csv('truncated_call.csv')
    d = to_arrays(header, data)

    fig, ax = plt.subplots()
    ax.plot(d['S'], d['StandardCentral'],
            color=COLORS['StandardCentral'], ls=STYLES['StandardCentral'],
            label=LABELS['StandardCentral'])

    # Fine-grid reference as "analytical" (ExpFit on same grid is smooth enough)
    ax.plot(d['S'], d['ExponentialFitting'],
            color=COLORS['Analytical'], ls=STYLES['Analytical'], lw=0.8,
            label='Reference (Exponential Fitting)')

    K = float(meta.get('strike', 50))
    U = float(meta.get('upper_barrier', 70))
    ax.set_xlim(K - 5, U + 10)
    ax.set_ylim(-5, max(d['StandardCentral'][(d['S'] > K) & (d['S'] < U + 5)]) * 1.3)
    ax.set_xlabel(r'Stock Price $S$')
    ax.set_ylabel(r'Option Value $V(S)$')
    ax.set_title(r'Spurious Oscillations in Crank-Nicolson ($\sigma=0.001$, $r=0.05$)')
    ax.legend(loc='upper left')
    ax.axvline(U, color='gray', ls=':', lw=0.5, alpha=0.5)
    fig.savefig(FIG_DIR / 'fig1_cn_oscillations.pdf')
    plt.close(fig)
    print('  Fig 1: CN oscillations')


# ---------------------------------------------------------------------------
# Figure 2: Three-scheme comparison (truncated call)
# ---------------------------------------------------------------------------

def fig2_three_scheme_truncated():
    meta, header, data = read_csv('truncated_call.csv')
    d = to_arrays(header, data)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4),
                                     gridspec_kw={'width_ratios': [2, 1]})

    for scheme in ['StandardCentral', 'ExponentialFitting', 'MilevTaglianiCN']:
        ax1.plot(d['S'], d[scheme],
                 color=COLORS[scheme], ls=STYLES[scheme],
                 label=LABELS[scheme])

    K = float(meta.get('strike', 50))
    U = float(meta.get('upper_barrier', 70))
    ax1.set_xlim(K - 5, U + 10)
    ax1.set_ylim(-5, 25)
    ax1.set_xlabel(r'Stock Price $S$')
    ax1.set_ylabel(r'Option Value $V(S)$')
    ax1.set_title(r'Three-Scheme Comparison: Truncated Call')
    ax1.legend(loc='upper left')
    ax1.axvline(U, color='gray', ls=':', lw=0.5, alpha=0.5)

    # Inset/zoom near U
    mask = (d['S'] > U - 3) & (d['S'] < U + 3)
    for scheme in ['StandardCentral', 'ExponentialFitting', 'MilevTaglianiCN']:
        ax2.plot(d['S'][mask], d[scheme][mask],
                 color=COLORS[scheme], ls=STYLES[scheme],
                 label=LABELS[scheme])
    ax2.set_xlabel(r'$S$')
    ax2.set_ylabel(r'$V(S)$')
    ax2.set_title(r'Zoom near $U=70$')
    ax2.axvline(U, color='gray', ls=':', lw=0.5, alpha=0.5)

    fig.tight_layout()
    fig.savefig(FIG_DIR / 'fig2_three_scheme_truncated.pdf')
    plt.close(fig)
    print('  Fig 2: Three-scheme truncated call')


# ---------------------------------------------------------------------------
# Figure 3: Moderate-vol discrete barrier
# ---------------------------------------------------------------------------

def fig3_barrier_moderate():
    meta, header, data = read_csv('barrier_moderate_vol.csv')
    d = to_arrays(header, data)

    # MC reference
    try:
        mc_meta, mc_header, mc_data = read_csv('mc_barrier_moderate_vol.csv')
        mc = to_arrays(mc_header, mc_data)
        has_mc = True
    except FileNotFoundError:
        has_mc = False

    fig, ax = plt.subplots()
    L, U = float(meta.get('lower_barrier', 95)), float(meta.get('upper_barrier', 110))
    mask = (d['S'] >= L - 1) & (d['S'] <= U + 1)

    for scheme in ['StandardCentral', 'ExponentialFitting', 'MilevTaglianiCN']:
        ax.plot(d['S'][mask], d[scheme][mask],
                color=COLORS[scheme], ls=STYLES[scheme],
                label=LABELS[scheme])

    if has_mc:
        ax.errorbar(mc['S0'], mc['price'], yerr=2*mc['standard_error'],
                     fmt='o', color=COLORS['MC'], ms=3, lw=0.8,
                     capsize=2, label=r'MC ($10^6$ paths)')

    ax.set_xlabel(r'Stock Price $S$')
    ax.set_ylabel(r'Option Value $V(S)$')
    ax.set_title(r'Discrete Double Barrier ($\sigma=0.25$, $K=100$, $L=95$, $U=110$)')
    ax.legend()
    ax.axvline(L, color='gray', ls=':', lw=0.5, alpha=0.5)
    ax.axvline(U, color='gray', ls=':', lw=0.5, alpha=0.5)
    fig.savefig(FIG_DIR / 'fig3_barrier_moderate.pdf')
    plt.close(fig)
    print('  Fig 3: Moderate-vol barrier')


# ---------------------------------------------------------------------------
# Figure 4: Low-vol discrete barrier
# ---------------------------------------------------------------------------

def fig4_barrier_lowvol():
    meta, header, data = read_csv('barrier_low_vol.csv')
    d = to_arrays(header, data)

    fig, ax = plt.subplots()
    L, U = float(meta.get('lower_barrier', 95)), float(meta.get('upper_barrier', 110))
    mask = (d['S'] >= L - 1) & (d['S'] <= U + 1)

    for scheme in ['StandardCentral', 'ExponentialFitting', 'MilevTaglianiCN']:
        ax.plot(d['S'][mask], d[scheme][mask],
                color=COLORS[scheme], ls=STYLES[scheme],
                label=LABELS[scheme])

    ax.axhline(0, color='red', ls=':', lw=0.5, alpha=0.7)
    ax.set_xlabel(r'Stock Price $S$')
    ax.set_ylabel(r'Option Value $V(S)$')
    ax.set_title(r'Low-Vol Discrete Barrier ($\sigma=0.001$, $\sigma^2 \ll r$)')
    ax.legend()
    ax.axvline(L, color='gray', ls=':', lw=0.5, alpha=0.5)
    ax.axvline(U, color='gray', ls=':', lw=0.5, alpha=0.5)
    fig.savefig(FIG_DIR / 'fig4_barrier_lowvol.pdf')
    plt.close(fig)
    print('  Fig 4: Low-vol barrier')


# ---------------------------------------------------------------------------
# Figure 5: Grid convergence
# ---------------------------------------------------------------------------

def fig5_convergence():
    meta, header, data = read_csv('grid_convergence.csv')
    d = to_arrays(header, data)

    fig, ax = plt.subplots()
    for scheme in ['StandardCentral', 'ExponentialFitting', 'MilevTaglianiCN']:
        err = d[f'err_{scheme}']
        ax.loglog(d['xGrid'], err,
                  color=COLORS[scheme], ls=STYLES[scheme],
                  marker='s', ms=4,
                  label=LABELS[scheme])

    # Reference slope lines
    x = d['xGrid']
    ax.loglog(x, 0.5 * (x[0]/x)**2, 'k:', lw=0.5, alpha=0.5, label=r'$O(h^2)$ slope')

    ax.set_xlabel(r'Spatial Grid Size ($N_x$)')
    ax.set_ylabel(r'$|V_{\mathrm{num}} - V_{\mathrm{BS}}|$')
    ax.set_title(r'Grid Convergence: European Call ($\sigma=0.20$)')
    ax.legend()
    fig.savefig(FIG_DIR / 'fig5_convergence.pdf')
    plt.close(fig)
    print('  Fig 5: Grid convergence')


# ---------------------------------------------------------------------------
# Figure 6: Effective diffusion comparison
# ---------------------------------------------------------------------------

def fig6_effective_diffusion():
    meta, header, data = read_csv('effective_diffusion.csv')
    d = to_arrays(header, data)

    fig, ax = plt.subplots()
    ax.plot(d['S'], d['aUsed_SC'],
            color=COLORS['StandardCentral'], ls=STYLES['StandardCentral'],
            label=LABELS['StandardCentral'])
    ax.plot(d['S'], d['aUsed_EF'],
            color=COLORS['ExponentialFitting'], ls=STYLES['ExponentialFitting'],
            label=LABELS['ExponentialFitting'])
    ax.plot(d['S'], d['aUsed_MT'],
            color=COLORS['MilevTaglianiCN'], ls=STYLES['MilevTaglianiCN'],
            label=LABELS['MilevTaglianiCN'])

    ax.set_xlabel(r'Stock Price $S$')
    ax.set_ylabel(r'Effective Diffusion Coefficient $a_{\mathrm{eff}}$')
    ax.set_title(r'Artificial Diffusion: Duffy vs Milev-Tagliani ($\sigma=0.001$)')
    ax.legend()
    ax.set_yscale('log')
    fig.savefig(FIG_DIR / 'fig6_effective_diffusion.pdf')
    plt.close(fig)
    print('  Fig 6: Effective diffusion')


# ---------------------------------------------------------------------------
# Figure 7: M-matrix off-diagonal visualization
# ---------------------------------------------------------------------------

def fig7_mmatrix():
    meta, header, data = read_csv('mmatrix_offdiag.csv')
    d = to_arrays(header, data)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Lower off-diagonal
    for scheme, lo_col in [('StandardCentral', 'lower_SC'),
                            ('ExponentialFitting', 'lower_EF'),
                            ('MilevTaglianiCN', 'lower_MT')]:
        ax1.plot(d['S'], d[lo_col],
                 color=COLORS[scheme], ls=STYLES[scheme],
                 label=LABELS[scheme])
    ax1.axhline(0, color='red', ls=':', lw=0.8, alpha=0.7)
    ax1.set_xlabel(r'Stock Price $S$')
    ax1.set_ylabel(r'Lower Off-Diagonal $a_{i,i-1}$')
    ax1.set_title('Lower Off-Diagonal Entries')
    ax1.legend(fontsize=7)

    # Upper off-diagonal
    for scheme, up_col in [('StandardCentral', 'upper_SC'),
                            ('ExponentialFitting', 'upper_EF'),
                            ('MilevTaglianiCN', 'upper_MT')]:
        ax2.plot(d['S'], d[up_col],
                 color=COLORS[scheme], ls=STYLES[scheme],
                 label=LABELS[scheme])
    ax2.axhline(0, color='red', ls=':', lw=0.8, alpha=0.7)
    ax2.set_xlabel(r'Stock Price $S$')
    ax2.set_ylabel(r'Upper Off-Diagonal $a_{i,i+1}$')
    ax2.set_title('Upper Off-Diagonal Entries')
    ax2.legend(fontsize=7)

    fig.suptitle(r'M-Matrix Property: Off-Diagonal Signs ($\sigma=0.001$)', y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'fig7_mmatrix.pdf')
    plt.close(fig)
    print('  Fig 7: M-matrix off-diagonals')


# ---------------------------------------------------------------------------
# Figure 8: Performance benchmark
# ---------------------------------------------------------------------------

def fig8_benchmark():
    meta, header, data = read_csv('benchmark.csv')
    d = to_arrays(header, data)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    schemes_data = {}
    for i, row in enumerate(data):
        xGrid = int(row[0])
        scheme = str(row[1]).strip()
        time_ms = float(row[2])
        error = float(row[4])
        if scheme not in schemes_data:
            schemes_data[scheme] = {'grids': [], 'times': [], 'errors': []}
        schemes_data[scheme]['grids'].append(xGrid)
        schemes_data[scheme]['times'].append(time_ms)
        schemes_data[scheme]['errors'].append(error)

    for scheme, sd in schemes_data.items():
        color = COLORS.get(scheme, '#999999')
        ls = STYLES.get(scheme, '-')
        label = LABELS.get(scheme, scheme)
        ax1.loglog(sd['grids'], sd['times'],
                   color=color, ls=ls, marker='s', ms=4, label=label)
        ax2.loglog(sd['times'], sd['errors'],
                   color=color, ls=ls, marker='s', ms=4, label=label)

    ax1.set_xlabel(r'Grid Size ($N_x$)')
    ax1.set_ylabel('Median Runtime (ms)')
    ax1.set_title('Runtime vs Grid Size')
    ax1.legend(fontsize=7)

    ax2.set_xlabel('Median Runtime (ms)')
    ax2.set_ylabel(r'$|V_{\mathrm{num}} - V_{\mathrm{BS}}|$')
    ax2.set_title('Runtime-Accuracy Tradeoff')
    ax2.legend(fontsize=7)

    fig.suptitle('Performance Benchmark: European Call', y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'fig8_benchmark.pdf')
    plt.close(fig)
    print('  Fig 8: Benchmark')


# ---------------------------------------------------------------------------
# Figure 9: xCothx / Peclet number
# ---------------------------------------------------------------------------

def fig9_xcothx():
    meta, header, data = read_csv('xcothx.csv')
    d = to_arrays(header, data)

    fig, ax = plt.subplots()
    ax.plot(d['Pe'], d['xCothx'], color=COLORS['ExponentialFitting'],
            ls='-', lw=1.5, label=r'$x \coth(x) = \rho(Pe)$')
    ax.plot(d['Pe'], d['abs_Pe'], color='gray', ls=':', lw=0.8,
            label=r'$|Pe|$ (upwind limit)')
    ax.axhline(1.0, color='gray', ls='--', lw=0.5, alpha=0.5,
               label=r'$\rho=1$ (central limit)')

    ax.set_xlabel(r'P\'eclet Number $Pe$')
    ax.set_ylabel(r'Fitting Factor $\rho = x\coth(x)$')
    ax.set_title(r'Duffy Fitting Factor: Regime Transition')
    ax.legend()
    ax.set_xlim(-20, 20)
    ax.set_ylim(0, 22)
    fig.savefig(FIG_DIR / 'fig9_xcothx.pdf')
    plt.close(fig)
    print('  Fig 9: xCothx/Peclet')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    if not DATA_DIR.exists():
        print(f"Error: Data directory {DATA_DIR} not found. Run generate_data first.")
        sys.exit(1)

    print("Generating figures...")

    fig9_xcothx()
    fig7_mmatrix()
    fig6_effective_diffusion()
    fig5_convergence()
    fig1_cn_oscillations()
    fig2_three_scheme_truncated()
    fig3_barrier_moderate()
    fig4_barrier_lowvol()
    fig8_benchmark()

    print(f"All 9 figures saved to {FIG_DIR}/")


if __name__ == '__main__':
    main()
