#!/usr/bin/env python3
"""
RLCR Loop Optimization Walkthrough — QuantLib CN Variant Code Review
Single-narrative, reference-grade visualization modeled after sol_001_optimization.png.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ─── Source-backed data ──────────────────────────────────────────

# ACs satisfied (categorical, from Codex review verdicts):
#   R0: AC-2,4,7 fully done; AC-1,3 partial; AC-5,6 not started  → 3/7
#   R1: +AC-5 done, AC-1,3,6 partial                              → 4/7
#   R2: +AC-6 done (full ctest), AC-1,3 still partial             → 5/7
#   R3: +AC-3 done (real aggregation), AC-1 still partial         → 6/7
#   R4: +AC-1 done (all citations), Codex verdict: "Complete"     → 7/7
#   Fin: simplification only, all ACs held                         → 7/7
acs_satisfied = np.array([3, 4, 5, 6, 7, 7])

round_labels = ["Round 0", "Round 1", "Round 2", "Round 3", "Round 4", "Finalize"]
round_x = np.arange(len(round_labels))

# Open Codex review issues (cumulative open at end of each round)
# R0: 8 raised, 0 resolved → 8 open
# R1: 5 new raised, 6 resolved from R0 backlog → 7 open
# R2: 4 new, 4 resolved → 7 open  (still had 3 from R1)
# R3: 2 new, 5 resolved → 4 open
# R4: 0 new, 4 resolved → 0 open
# Fin: 0/0 → 0
open_issues = np.array([8, 7, 7, 4, 0, 0])

# Files changed per commit (source: git diff --stat per round)
files_changed = [24, 13, 8, 10, 4, 4]

# Codex verdict per round
verdicts = ["Fail", "Fail", "Fail", "Fail", "Pass", "N/A"]

# Phase regions (start, end, label, color)
phases = [
    (-0.45, 0.45, "Foundation &\nMath Verification", "#BBDEFB"),
    (0.55,  1.45, "Observability &\nCommenting",     "#FFE0B2"),
    (1.55,  2.45, "Coverage &\nVerification",        "#E1BEE7"),
    (2.55,  3.45, "Aggregation &\nCitations",        "#C8E6C9"),
    (3.55,  4.45, "Discriminating\nTests",           "#FFF9C4"),
    (4.55,  5.45, "Simplify &\nFinalize",            "#F8BBD0"),
]

# Milestone annotations (round_idx, label, marker_type)
# marker_type: 'kept' = improvement, 'neutral' = neutral, 'pass' = final pass
milestones_kept = [
    (0, "Math verified\nDead API removed\n3/7 ACs"),
    (1, "Solver observability\nMemo + test headers\n4/7 ACs"),
    (2, "All exit paths\n3 regression tests\n5/7 ACs"),
    (3, "Real aggregation\nSub-solve OR logic\n6/7 ACs"),
    (4, "Discriminating test\nAll citations\n7/7 ACs"),
]
milestone_final = (5, "schemeName()\nCode dedup\n7/7 ACs")

# ─── Figure (wide landscape, single dominant panel) ──────────────

fig, ax = plt.subplots(figsize=(19, 8.5))
fig.subplots_adjust(top=0.82, bottom=0.22, left=0.065, right=0.925)

# ─── Subtitle line ───────────────────────────────────────────────
# ─── Title ───────────────────────────────────────────────────────
fig.suptitle(
    "QuantLib CN Variant Fork: Code Review RLCR Optimization Path\n"
    "(Spatial Discretization Scheme Review & Observability Retrofit)",
    fontsize=15, fontweight='bold', y=0.98)

fig.text(0.5, 0.905,
         "5 rounds  |  7 acceptance criteria  |  "
         "29 changed source/header files  |  5 regression tests  |  "
         "Full ctest: 204 s",
         ha='center', fontsize=10.5, color='#555555')

# ─── Phase bands ─────────────────────────────────────────────────
for (x0, x1, lbl, col) in phases:
    ax.axvspan(x0, x1, alpha=0.22, color=col, zorder=0)
    ax.text((x0 + x1) / 2, 8.1, lbl, ha='center', va='bottom',
            fontsize=9, fontstyle='italic', color='#444444',
            fontweight='medium')

# ─── Reference baseline ─────────────────────────────────────────
ax.axhline(y=7.0, color='#43A047', linestyle='--', alpha=0.45,
           linewidth=1.2, zorder=1)
ax.text(5.58, 7.08, "all 7 ACs\nsatisfied", fontsize=8.5,
        color='#2E7D32', ha='left', va='bottom', fontstyle='italic')

# ─── Best-so-far frontier (step line) ───────────────────────────
best_so_far = np.maximum.accumulate(acs_satisfied)
ax.step(round_x, best_so_far, where='post', color='#1B5E20',
        linewidth=2, alpha=0.35, zorder=2, linestyle='-')

# ─── Main trajectory ────────────────────────────────────────────
ax.plot(round_x, acs_satisfied, '-', color='#1B5E20', linewidth=2.8,
        zorder=4)

# ─── Per-round markers ──────────────────────────────────────────

# Kept (improvement) markers — Rounds 0-4
for rx, lbl in milestones_kept:
    ax.plot(rx, acs_satisfied[rx], 'o', color='#1B5E20',
            markersize=12, zorder=5, markeredgewidth=1.5,
            markeredgecolor='white')

# Neutral marker — Finalize (simplify, no AC change)
rx, lbl = milestone_final
ax.plot(rx, acs_satisfied[rx], 'o', color='#9E9E9E',
        markersize=11, zorder=5, markeredgewidth=1.5,
        markeredgecolor='white')

# ─── Annotation boxes ───────────────────────────────────────────
annot_cfg = [
    # (round, y_offset, label_short, color)
    (0, -1.6,  "Foundation\n3/7",                "#D32F2F"),
    (1, +1.1,  "Observability\n4/7",             "#D32F2F"),
    (2, -1.5,  "Full ctest\n5/7",                "#D32F2F"),
    (3, +1.1,  "Real aggregation\n6/7",          "#D32F2F"),
    (4, -1.3,  "All ACs\n7/7",                   "#2E7D32"),
    (5, +1.0,  "Simplified\n7/7",                "#1565C0"),
]

for rx, yoff, lbl, col in annot_cfg:
    y = acs_satisfied[rx]
    ax.annotate(
        lbl, (rx, y),
        textcoords="offset points",
        xytext=(0, yoff * 22),
        ha='center', fontsize=9.5, fontweight='bold', color=col,
        arrowprops=dict(arrowstyle='->', color=col, lw=1.2),
        bbox=dict(boxstyle='round,pad=0.35', fc='white', ec=col,
                  alpha=0.92),
        zorder=6,
    )

# ─── Codex verdict badges ───────────────────────────────────────
for i in range(len(round_labels)):
    if verdicts[i] == "Fail":
        ax.plot(i, 0.5, 'X', color='#D32F2F', markersize=11, zorder=5,
                markeredgewidth=1.5)
    elif verdicts[i] == "Pass":
        ax.plot(i, 0.5, '*', color='#2E7D32', markersize=14, zorder=5,
                markeredgewidth=0.8)

# ─── Open issues secondary line (subtle) ────────────────────────
ax2 = ax.twinx()
ax2.bar(round_x, open_issues, width=0.22, color='#EF9A9A', alpha=0.5,
        zorder=1, label='Open review issues')
for i, v in enumerate(open_issues):
    if v > 0:
        ax2.text(i, v + 0.25, str(v), ha='center', fontsize=8,
                 color='#C62828', fontweight='bold')
ax2.set_ylabel("Open Review Issues", fontsize=10, color='#C62828')
ax2.set_ylim(0, 12)
ax2.tick_params(axis='y', colors='#C62828')
ax2.legend(loc='upper right', fontsize=9, framealpha=0.9)

# ─── Axes formatting ────────────────────────────────────────────
ax.set_xlim(-0.55, 5.7)
ax.set_ylim(-0.2, 9.5)
ax.set_xticks(round_x)
ax.set_xticklabels(round_labels, fontsize=11, fontweight='medium')
ax.set_ylabel("Acceptance Criteria Satisfied (of 7)", fontsize=11)
ax.set_xlabel("Review Round", fontsize=11)
ax.grid(axis='y', alpha=0.25, zorder=0)

# ─── Legend (main) ───────────────────────────────────────────────
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='#1B5E20', markersize=9,
               linewidth=2.5, markeredgecolor='white',
               label='Kept (AC improvement)'),
    plt.Line2D([0], [0], marker='o', color='#9E9E9E', markersize=9,
               linewidth=0, markeredgecolor='white',
               label='Neutral (simplification)'),
    plt.Line2D([0], [0], marker='X', color='#D32F2F', markersize=9,
               linewidth=0, markeredgewidth=1.5,
               label='Codex verdict: not complete'),
    plt.Line2D([0], [0], marker='*', color='#2E7D32', markersize=11,
               linewidth=0, label='Codex verdict: complete'),
    plt.Line2D([0], [0], color='#1B5E20', linestyle='--', alpha=0.5,
               label='Best so far (frontier)'),
    plt.Line2D([0], [0], color='#43A047', linestyle='--', alpha=0.5,
               label='All ACs satisfied'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=8.5,
          framealpha=0.92, ncol=3,
          bbox_to_anchor=(0.0, 1.0))

# ─── Milestone detail boxes along bottom ─────────────────────────
detail_labels = [
    "Math verified (Codex)\nDead API removed\n24 files commented",
    "Solver op_ forwarding\n.cpp + test headers\nReview memo (6 sect.)",
    "reportSpatialScheme()\n3 regression tests\nctest: 206.60 s",
    "Sub-instrument OR\nKnock-in + t=0 tests\nctest: 204.21 s",
    "Discriminating test\nall.hpp citations\nctest: 204.42 s",
    "schemeName() dedup\nRemove payoff2\nctest: 204.59 s",
]
for i, lbl in enumerate(detail_labels):
    fig.text(0.065 + (i + 0.5) * (0.86 / 6), 0.03, lbl,
             ha='center', va='bottom', fontsize=7.5,
             color='#444', fontstyle='italic',
             bbox=dict(boxstyle='round,pad=0.3', fc='#FAFAFA',
                       ec='#BDBDBD', alpha=0.92))

# ─── Info box (bottom right, compact) ────────────────────────────
info = (
    "QuantLib v1.42 CN Variant Fork\n"
    "3 schemes: StandardCentral | ExponentialFitting | MilevTaglianiCN\n"
    "Papers: [MT10]  [Duffy04]  [Ballabio20]\n"
    "5 rounds + Finalize  |  Codex reviewer  |  Full ctest: 204 s"
)
ax.text(5.65, 2.2, info, fontsize=7.5, ha='right', va='bottom',
        fontfamily='monospace', color='#666',
        bbox=dict(boxstyle='round,pad=0.5', fc='white', ec='#BBB',
                  alpha=0.92))

# ─── Save ────────────────────────────────────────────────────────
output = '/home/jakeshea/QuantLib_huatai2/rlcr_optimization_walkthrough.png'
fig.savefig(output, dpi=180, bbox_inches='tight', facecolor='white')
print(f"Saved to {output}")
