#!/usr/bin/env python3
"""
Publication figures for blood pressure, PWV, and heart rate: one box plot
with individual trajectories per variable, plus a 3-panel PWV figure
(A = CF-PWV, B = CB-PWV, C = CR-PWV). The AIx@75 figure is made by
adhd_graphing_aix75.py.

Each panel shows the RM-ANOVA P and partial eta-squared, n, and brackets
for Holm-Bonferroni significant pairs (drawn only when the omnibus test is
significant).

Input:  NIHem_Data_Cleaned_3.xlsx (sheet tonometry_cleaned), same folder
Output: Graphs_pub/
"""

import os
from itertools import combinations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pingouin as pg
from scipy import stats

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(SCRIPT_DIR, 'NIHem_Data_Cleaned_3.xlsx')
SHEET = 'tonometry_cleaned'
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'Graphs_pub')

MISSING_CODES = [-9999, -9998, -9997]

TP_MAP = {
    'Baseline': 'Baseline',
    'Condition 1': '30min',
    'Condition 2': '1hr',
    'Condition 3': '5hr'
}
TP_ORDER = ['Baseline', '30min', '1hr', '5hr']
TP_LABELS = ['Baseline', '30 min', '1 h', '5 h']

# (column, y-axis label, output file)
MEASUREMENTS = [
    ('systolic.blood.pressure.mmhg',   'SBP (mmHg)',    'SBP.png'),
    ('diastolic.blood.pressure.mmhg',  'DBP (mmHg)',    'DBP.png'),
    ('mean.arterial.pressure.mmhg',    'MAP (mmHg)',    'MAP.png'),
    ('carotid.femoral.pwv.cm.s',       'CF-PWV (cm/s)', 'CF_PWV.png'),
    ('carotid.brachial.pwv.cm.s',      'CB-PWV (cm/s)', 'CB_PWV.png'),
    ('carotid.radial.pwv.cm.s',        'CR-PWV (cm/s)', 'CR_PWV.png'),
    ('hr_carotid',                     'HR (bpm)',      'HR.png'),
]
PWV_COLUMNS = ['carotid.femoral.pwv.cm.s', 'carotid.brachial.pwv.cm.s', 'carotid.radial.pwv.cm.s']


def sig_stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.'


def rm_anova(wide_df):
    """One-way RM-ANOVA on a participants x timepoints table.

    P is Greenhouse-Geisser corrected when Mauchly's test is significant.
    """
    n, k = wide_df.shape
    vals = wide_df.values
    grand_mean = vals.mean()

    ss_time = n * np.sum((vals.mean(axis=0) - grand_mean) ** 2)
    ss_subjects = k * np.sum((vals.mean(axis=1) - grand_mean) ** 2)
    ss_total = np.sum((vals - grand_mean) ** 2)
    ss_error = ss_total - ss_time - ss_subjects

    df_time = k - 1
    df_error = (k - 1) * (n - 1)
    f_stat = (ss_time / df_time) / (ss_error / df_error)
    p_uncorrected = stats.f.sf(f_stat, df_time, df_error)

    subject = wide_df.index.name
    long_df = wide_df.reset_index().melt(id_vars=subject, var_name='timepoint', value_name='value')
    aov = pg.rm_anova(data=long_df, dv='value', within='timepoint', subject=subject,
                      detailed=True, correction=True)
    # older pingouin versions use hyphens in column names (p-GG-corr)
    aov.columns = [c.replace('-', '_') for c in aov.columns]
    sphericity_met = bool(aov.loc[0, 'sphericity'])

    return {
        'p': p_uncorrected if sphericity_met else float(aov.loc[0, 'p_GG_corr']),
        'partial_eta_sq': ss_time / (ss_time + ss_error),
    }


def holm_bonferroni(p_values):
    """Holm-adjusted P values (step-down, monotone)."""
    m = len(p_values)
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * m
    running_max = 0.0
    for rank, (orig_idx, p) in enumerate(indexed):
        running_max = max(running_max, p * (m - rank))
        adjusted[orig_idx] = min(running_max, 1.0)
    return adjusted


def significant_pairs(wide_df):
    """(i, j, stars) for timepoint pairs significant after Holm-Bonferroni."""
    pairs = list(combinations(range(len(TP_ORDER)), 2))
    p_raw = [stats.ttest_rel(wide_df[TP_ORDER[i]].values, wide_df[TP_ORDER[j]].values)[1]
             for i, j in pairs]
    p_holm = holm_bonferroni(p_raw)
    return [(i, j, sig_stars(p)) for (i, j), p in zip(pairs, p_holm) if p < 0.05]


def analyse(col):
    wide = raw.pivot_table(index='participant', columns='timepoint_clean', values=col)
    wide = wide[TP_ORDER].dropna()
    res = rm_anova(wide)
    sig_pairs = significant_pairs(wide) if res['p'] < 0.05 else []
    return wide, res, sig_pairs


def add_bracket(ax, x1, x2, y, h, text, fs=13):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.3, color='black')
    ax.text((x1 + x2) / 2, y + h, text, ha='center', va='bottom', fontsize=fs, fontweight='bold')


def add_panel_labels(fig, axes_labels, fontsize=20, dy=0.012):
    """Bold panel letters, left-aligned with each panel's y-axis labels.
    Call after the layout is final."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for ax, label in axes_labels:
        bboxes = []
        if ax.yaxis.label.get_text():
            bboxes.append(ax.yaxis.label.get_window_extent(renderer=renderer))
        bboxes += [t.get_window_extent(renderer=renderer)
                   for t in ax.get_yticklabels() if t.get_text()]
        x0_px = min(b.x0 for b in bboxes) if bboxes else ax.get_window_extent(renderer=renderer).x0
        x0_fig = inv.transform((x0_px, 0))[0]
        fig.text(x0_fig, ax.get_position().y1 + dy, label,
                 fontsize=fontsize, fontweight='bold', va='bottom', ha='left')


def plot_panel(ax, wide, res, sig_pairs, ylabel, fontsize_anova=10):
    bp = ax.boxplot(
        [wide[tp].values for tp in TP_ORDER],
        positions=[0, 1, 2, 3], widths=0.55, patch_artist=True,
        showmeans=False, showfliers=False, whis=(0, 100),
        medianprops=dict(color='black', linewidth=2.2),
        whiskerprops=dict(color='black', linewidth=1.3),
        capprops=dict(color='black', linewidth=0)
    )
    for patch in bp['boxes']:
        patch.set_facecolor('none')
        patch.set_edgecolor('black')
        patch.set_linewidth(1.4)
    for cap in bp['caps']:
        cap.set_visible(False)

    # dashed mean line
    for tp_idx, tp in enumerate(TP_ORDER):
        ax.plot([tp_idx - 0.27, tp_idx + 0.27], [wide[tp].mean()] * 2,
                color='black', linestyle='--', linewidth=1.8, zorder=5)

    for participant in wide.index:
        ax.plot([0, 1, 2, 3], wide.loc[participant].values,
                color='gray', alpha=0.4, linewidth=0.9, zorder=1)

    p, eta = res['p'], res['partial_eta_sq']
    p_str = f'P < 0.001, η²p = {eta:.3f}' if p < 0.001 else f'P = {p:.3f}, η²p = {eta:.3f}'
    ax.text(0.02, 0.97, p_str,
            transform=ax.transAxes, fontsize=fontsize_anova, va='top', fontstyle='italic',
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                      edgecolor='black', linewidth=1.2, alpha=0.95))
    ax.text(0.98, 0.97, f'n = {len(wide)}', transform=ax.transAxes,
            fontsize=fontsize_anova + 1, va='top', ha='right', fontstyle='italic')

    ymax = wide.values.max()
    yr = ymax - wide.values.min()
    if sig_pairs:
        by = ymax + yr * 0.06
        bs = yr * 0.09
        for ib, (i, j, stars) in enumerate(sorted(sig_pairs, key=lambda x: x[1] - x[0])):
            add_bracket(ax, i, j, by + ib * bs, bs * 0.25, stars, fs=fontsize_anova + 2)
        ax.set_ylim(top=by + len(sig_pairs) * bs + bs * 0.7)

    # headroom so the stats box does not overlap the data
    ymin_ax, ymax_ax = ax.get_ylim()
    ax.set_ylim(top=ymax_ax + (ymax_ax - ymin_ax) * 0.12)

    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(TP_LABELS, fontsize=11)
    ax.set_xlabel('Timepoint', fontsize=12, fontweight='bold', labelpad=8)
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold', labelpad=8)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.1)
    ax.grid(False)
    ax.tick_params(labelsize=10)


# Load and clean
print(f"Loading: {DATA_FILE}")
raw = pd.read_excel(DATA_FILE, sheet_name=SHEET)
raw = raw.replace(MISSING_CODES, np.nan)

raw['timepoint_clean'] = raw['timepoint'].map(TP_MAP)
raw['participant'] = 'P' + raw['sequence_no'].astype(str).str.zfill(2)
raw['hr_carotid'] = 60000 / raw['carotid.waveform.rr.interval.ms']

os.makedirs(OUTPUT_DIR, exist_ok=True)
results = {col: analyse(col) for col, _, _ in MEASUREMENTS}
ylabels = {col: ylabel for col, ylabel, _ in MEASUREMENTS}


# One figure per variable
for col, ylabel, filename in MEASUREMENTS:
    fig, ax = plt.subplots(figsize=(6.5, 7.5))
    plot_panel(ax, *results[col], ylabel)
    fname = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(fname, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {fname}")


# 3-panel PWV figure
fig, axes = plt.subplots(1, 3, figsize=(20, 7.5))
for ax, col in zip(axes, PWV_COLUMNS):
    plot_panel(ax, *results[col], ylabels[col], fontsize_anova=9)

plt.tight_layout()
plt.subplots_adjust(wspace=0.22)
add_panel_labels(fig, list(zip(axes, ['A', 'B', 'C'])))

fname = os.path.join(OUTPUT_DIR, 'PWV_combined.png')
fig.savefig(fname, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"  Saved: {fname}")

print(f"\nAll figures saved to: {OUTPUT_DIR}")
