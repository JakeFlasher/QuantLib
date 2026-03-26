#!/usr/bin/env python3
"""
RLCR Loop Optimization Walkthrough — QuantLib CN Variant Code Review
Generates a professional, academic-style visual showing the detailed
review/fix/verify optimization path across 5 rounds + finalize.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ─── Data ────────────────────────────────────────────────────────

# Acceptance criteria (7 total)
AC_LABELS = [
    "AC-1: Paper-annotated comments",
    "AC-2: Dead API removal",
    "AC-3: Fallback observability",
    "AC-4: Build integration",
    "AC-5: Review findings memo",
    "AC-6: Full test suite pass",
    "AC-7: Math verification (Codex)",
]

# AC completion status per round (0=not started, 0.5=partial, 1=complete)
# Rows: AC-1..AC-7, Columns: R0, R1, R2, R3, R4, Finalize
ac_status = np.array([
    # R0    R1    R2    R3    R4    Fin
    [0.50, 0.70, 0.85, 0.95, 1.00, 1.00],  # AC-1: comments
    [1.00, 1.00, 1.00, 1.00, 1.00, 1.00],  # AC-2: dead API
    [0.40, 0.60, 0.75, 1.00, 1.00, 1.00],  # AC-3: observability
    [1.00, 1.00, 1.00, 1.00, 1.00, 1.00],  # AC-4: build
    [0.00, 1.00, 1.00, 1.00, 1.00, 1.00],  # AC-5: memo
    [0.00, 0.50, 1.00, 1.00, 1.00, 1.00],  # AC-6: tests
    [1.00, 1.00, 1.00, 1.00, 1.00, 1.00],  # AC-7: math
])

# Aggregate "AC completion score" (out of 7)
ac_score = ac_status.sum(axis=0)

# Round labels and phases
round_labels = ["Round 0", "Round 1", "Round 2", "Round 3", "Round 4", "Finalize"]
round_x = np.arange(len(round_labels))

# Files changed per round
files_changed = [24, 13, 8, 10, 4, 4]

# Tests passing per round (cumulative new test evidence)
tests_status = [
    "Partial\n(49 CN)",
    "Partial\n(49 CN + barrier)",
    "Full suite\n(206.60s)",
    "Full suite\n(204.21s)",
    "Full suite\n(204.42s)",
    "Full suite\n(204.59s)",
]

# Codex review verdicts
codex_verdicts = [
    "Not complete",   # R0
    "Not complete",   # R1
    "Not complete",   # R2
    "Not complete",   # R3
    "Complete",       # R4
    "N/A",            # Finalize
]

# Key milestones (round_idx, short_label)
milestones = [
    (0, "Math verified\nDead API removed\n24 files commented"),
    (1, "Solver observability\n.cpp comments\nTest headers + memo"),
    (2, "All exit paths\n3 regression tests\nFull ctest pass"),
    (3, "Real aggregation\nSub-solve OR logic\n5 regression tests"),
    (4, "Discriminating test\nUmbrella citations"),
    (5, "schemeName()\nCode dedup"),
]

# Phase regions (start_round, end_round, label, color)
phases = [
    (0, 0, "Foundation &\nMath Verification", "#E3F2FD"),
    (1, 1, "Observability &\nCommenting", "#FFF3E0"),
    (2, 2, "Coverage &\nVerification", "#F3E5F5"),
    (3, 3, "Aggregation &\nCitations", "#E8F5E9"),
    (4, 4, "Discriminating\nTests", "#FFFDE7"),
    (5, 5, "Simplify &\nFinalize", "#FCE4EC"),
]

# Codex issue counts per round
issues_found     = [8, 5, 4, 2, 0, 0]     # issues found by Codex
issues_resolved  = [0, 6, 4, 3, 2, 3]     # issues resolved in that round

# ─── Figure ──────────────────────────────────────────────────────

fig, axes = plt.subplots(3, 1, figsize=(16, 13.5),
                          gridspec_kw={'height_ratios': [3.5, 2, 1.8]},
                          sharex=True)
fig.subplots_adjust(hspace=0.12, top=0.90, bottom=0.06, left=0.08, right=0.92)

# Title
fig.suptitle(
    "QuantLib CN Variant Fork: Code Review RLCR Optimization Path\n"
    "Paper-Annotated Commenting, Fallback Observability & Build Integration",
    fontsize=14, fontweight='bold', y=0.97,
)
fig.text(0.5, 0.925,
         "5 rounds  |  7 acceptance criteria  |  30 files touched  |  "
         "5 regression tests  |  Full ctest: 204s",
         ha='center', fontsize=10, color='#555555')

# ─── Panel 1: AC Completion Trajectory ───────────────────────────
ax1 = axes[0]

# Phase background bands
for (rs, re, lbl, col) in phases:
    ax1.axvspan(rs - 0.4, re + 0.4, alpha=0.25, color=col, zorder=0)
    ax1.text((rs + re) / 2, 7.35, lbl, ha='center', va='bottom',
             fontsize=7.5, fontstyle='italic', color='#333333')

# Plot AC score (main line)
ax1.plot(round_x, ac_score, 'o-', color='#1B5E20', linewidth=2.5,
         markersize=10, zorder=5, label='AC completion score (of 7)')

# Annotate each point
for i, (x, y) in enumerate(zip(round_x, ac_score)):
    verd_color = '#D32F2F' if codex_verdicts[i] == "Not complete" else '#2E7D32'
    if codex_verdicts[i] != "N/A":
        ax1.annotate(
            f"{y:.1f}/7\n{codex_verdicts[i]}",
            (x, y), textcoords="offset points",
            xytext=(0, 14), ha='center', fontsize=8,
            fontweight='bold', color=verd_color,
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=verd_color,
                      alpha=0.85),
        )
    else:
        ax1.annotate(
            f"{y:.1f}/7\nFinalized",
            (x, y), textcoords="offset points",
            xytext=(0, 14), ha='center', fontsize=8,
            fontweight='bold', color='#1565C0',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#1565C0',
                      alpha=0.85),
        )

# Reference line at 7
ax1.axhline(y=7.0, color='#2E7D32', linestyle='--', alpha=0.4, linewidth=1)
ax1.text(5.45, 7.05, "All ACs\nsatisfied", fontsize=7, color='#2E7D32',
         ha='left', va='bottom')

ax1.set_ylabel("Acceptance Criteria Score", fontsize=10)
ax1.set_ylim(2.5, 8.0)
ax1.set_xlim(-0.5, 5.7)
ax1.legend(loc='lower right', fontsize=9)
ax1.grid(axis='y', alpha=0.3)

# ─── Panel 2: Stacked AC heatmap ────────────────────────────────
ax2 = axes[1]

im = ax2.imshow(ac_status, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1,
                interpolation='nearest')

# Text annotations inside cells
for i in range(ac_status.shape[0]):
    for j in range(ac_status.shape[1]):
        val = ac_status[i, j]
        txt = {0.0: '--', 0.40: '40%', 0.50: '50%', 0.60: '60%',
               0.70: '70%', 0.75: '75%', 0.85: '85%', 0.95: '95%',
               1.0: '\u2713'}.get(val, f'{val:.0%}')
        color = 'white' if val < 0.4 or val >= 0.85 else 'black'
        ax2.text(j, i, txt, ha='center', va='center', fontsize=8,
                 fontweight='bold', color=color)

ax2.set_yticks(range(len(AC_LABELS)))
ax2.set_yticklabels(AC_LABELS, fontsize=8)
ax2.set_xticks(round_x)
ax2.set_xticklabels(round_labels, fontsize=9)
ax2.set_title("Per-AC Completion Matrix", fontsize=10, fontweight='bold', pad=6)

# Colorbar
cbar = fig.colorbar(im, ax=ax2, fraction=0.02, pad=0.015)
cbar.set_label("Completion", fontsize=8)

# ─── Panel 3: Effort + Issues bars ──────────────────────────────
ax3 = axes[2]
ax3b = ax3.twinx()

bar_w = 0.30
x_off = bar_w / 2

# Files changed (left y-axis)
bars1 = ax3.bar(round_x - x_off, files_changed, bar_w, color='#1565C0',
                alpha=0.75, label='Files changed', zorder=3)

# Issues resolved (left y-axis, stacked next to files)
bars2 = ax3.bar(round_x + x_off, issues_resolved, bar_w, color='#2E7D32',
                alpha=0.75, label='Issues resolved', zorder=3)

# Issues found by Codex (right y-axis, line)
line = ax3b.plot(round_x, issues_found, 's--', color='#D32F2F', linewidth=1.5,
                 markersize=7, label='Issues raised by Codex', zorder=4)

# Annotate files bars
for i, v in enumerate(files_changed):
    ax3.text(round_x[i] - x_off, v + 0.3, str(v), ha='center',
             fontsize=7.5, fontweight='bold', color='#1565C0')

# Annotate resolved bars
for i, v in enumerate(issues_resolved):
    if v > 0:
        ax3.text(round_x[i] + x_off, v + 0.3, str(v), ha='center',
                 fontsize=7.5, fontweight='bold', color='#2E7D32')

# Annotate issues-found markers
for i, v in enumerate(issues_found):
    if v > 0:
        ax3b.text(round_x[i], v + 0.35, str(v), ha='center',
                  fontsize=7.5, fontweight='bold', color='#D32F2F')

ax3.set_ylabel("Count", fontsize=10, color='#333')
ax3b.set_ylabel("Codex Issues Found", fontsize=10, color='#D32F2F')
ax3.set_ylim(0, 28)
ax3b.set_ylim(0, 10)
ax3.set_xticks(round_x)
ax3.set_xticklabels(round_labels, fontsize=9)
ax3.grid(axis='y', alpha=0.2)

# Combined legend
handles1, labels1 = ax3.get_legend_handles_labels()
handles2, labels2 = ax3b.get_legend_handles_labels()
ax3.legend(handles1 + handles2, labels1 + labels2, loc='upper right',
           fontsize=8, ncol=3)

ax3.set_title("Effort & Issue Convergence", fontsize=10, fontweight='bold', pad=6)

# ─── Milestone annotations along bottom ─────────────────────────
for (rx, lbl) in milestones:
    ax3.text(rx, -5.5, lbl, ha='center', va='top', fontsize=6.5,
             color='#333', fontstyle='italic',
             bbox=dict(boxstyle='round,pad=0.3', fc='#F5F5F5', ec='#BDBDBD',
                       alpha=0.9))

# ─── Info box ────────────────────────────────────────────────────
info_text = (
    "QuantLib v1.42 CN Variant Fork\n"
    "3 spatial schemes: StandardCentral, ExponentialFitting, MilevTaglianiCN\n"
    "Papers: [MT10], [Duffy04], [Ballabio20]\n"
    "5 RLCR rounds + Finalize  |  Codex reviewer  |  Full ctest: 204s"
)
fig.text(0.92, 0.06, info_text, fontsize=7, ha='right', va='bottom',
         fontfamily='monospace', color='#666',
         bbox=dict(boxstyle='round,pad=0.5', fc='white', ec='#CCC', alpha=0.9))

# ─── Save ────────────────────────────────────────────────────────
output = '/home/jakeshea/QuantLib_huatai2/rlcr_optimization_walkthrough.png'
fig.savefig(output, dpi=180, bbox_inches='tight', facecolor='white')
print(f"Saved to {output}")
