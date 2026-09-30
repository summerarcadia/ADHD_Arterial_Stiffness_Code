"""
Publication figure for carotid distensibility coefficient (DC), left
carotid. Uses the data and statistics functions in RM_ANOVA_CD_PP.py.

Box plot with individual trajectories, the RM-ANOVA P and partial
eta-squared, n, and brackets for pairs significant after Holm-Bonferroni
(drawn only when the omnibus test is significant; stars use the
Holm-adjusted P).

Input:  same as RM_ANOVA_CD_PP.py
Output: Graphs_pub/RM_ANOVA_CD_PP_boxplot_left.png (300 dpi) and .pdf
"""

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import RM_ANOVA_CD_PP as rm

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']

OUTPUT_DIR = os.path.join(rm.HERE, 'Graphs_pub')
TIMEPOINTS = rm.TIMEPOINTS
TICK_LABELS = ['Baseline', '30 min', '1 h', '5 h']


def p_to_stars(p):
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return 'n.s.'


def get_stats(df):
    wide, long_df = rm.complete_cases(df)
    _, s = rm.rm_anova(long_df)

    sig_pairs = []
    if s['P'] < 0.05:
        ph = rm.posthoc_paired_holm(wide)
        for comparison, p_adj in zip(ph['Comparison'], ph['p_holm_adjusted']):
            if p_adj < 0.05:
                a, b = comparison.split(' vs ')
                sig_pairs.append((a, b, p_adj))
    return wide, s, sig_pairs


def make_figure(wide, p_report, pes, sig_pairs):
    n = len(wide)
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    for _, r in wide.iterrows():
        ax.plot(range(4), r[TIMEPOINTS].values, color='grey', alpha=0.35,
                linewidth=0.8, zorder=1)

    box_data = [wide[tp].values for tp in TIMEPOINTS]
    bp = ax.boxplot(
        box_data, positions=range(4), widths=0.5, patch_artist=True,
        showmeans=True, meanline=True, showfliers=False, whis=(0, 100),
        meanprops=dict(color='black', linewidth=1.5, linestyle='--'),
        medianprops=dict(color='black', linewidth=1.5),
        whiskerprops=dict(color='black', linewidth=1),
        capprops=dict(color='black', linewidth=0), zorder=2,
    )
    for patch in bp['boxes']:
        patch.set_facecolor('none')
        patch.set_edgecolor('black')
        patch.set_linewidth(1.5)
    for cap in bp['caps']:
        cap.set_visible(False)

    anova_text = (f'P < 0.001, η²p = {pes:.3f}' if p_report < 0.001
                  else f'P = {p_report:.3f}, η²p = {pes:.3f}')
    ax.text(0.02, 0.97, anova_text, transform=ax.transAxes, fontsize=9,
            va='top', fontstyle='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='grey', alpha=0.8))
    ax.text(0.98, 0.97, f'n = {n}', transform=ax.transAxes, fontsize=10,
            ha='right', va='top', fontstyle='italic')

    y_max = max(d.max() for d in box_data)
    y_min = min(d.min() for d in box_data)
    rng = y_max - y_min
    gap = rng * 0.06
    top = y_max + rng * 0.10
    for i, (a, b, p_adj) in enumerate(sig_pairs):
        x1, x2 = TIMEPOINTS.index(a), TIMEPOINTS.index(b)
        y = top + i * gap
        ax.plot([x1, x1, x2, x2], [y - gap * 0.25, y, y, y - gap * 0.25],
                color='black', linewidth=1, zorder=3)
        ax.text((x1 + x2) / 2, y + gap * 0.05, p_to_stars(p_adj),
                ha='center', va='bottom', fontsize=11, fontweight='bold', zorder=3)
    # headroom for brackets and the stats box
    ax.set_ylim(bottom=y_min - rng * 0.05,
                top=top + max(len(sig_pairs), 1) * gap + rng * 0.12)

    ax.set_xticks(range(4))
    ax.set_xticklabels(TICK_LABELS, fontsize=11)
    ax.set_xlabel('Timepoint', fontsize=12, labelpad=8)
    ax.set_ylabel('Distensibility Coefficient\n(10$^{-3}$ $\\cdot$ mmHg$^{-1}$)',
                  fontsize=12, labelpad=8)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.1)
    ax.tick_params(labelsize=10)
    plt.tight_layout()
    return fig


def main():
    print('Calculating DC from carotid PP')
    lookup = rm.build_carotid_pp_lookup(rm.NIHEM_FILE)
    df = rm.extract_dc_values(rm.EXCEL_FILE, lookup)

    wide, s, sig_pairs = get_stats(df)
    print(f"\nRM-ANOVA: P = {s['P']:.4f} "
          f"({'GG-corrected' if s['GG_corrected'] else 'uncorrected, sphericity met'}), "
          f"partial eta-sq = {s['partial_eta_sq']:.3f}, n = {len(wide)}")
    print('Significant pairs (Holm-adjusted P):')
    for a, b, p_adj in sig_pairs:
        print(f'  {a} vs {b}: P = {p_adj:.4f} {p_to_stars(p_adj)}')
    if not sig_pairs:
        print('  none')

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fig = make_figure(wide, s['P'], s['partial_eta_sq'], sig_pairs)
    for ext in ('png', 'pdf'):
        path = os.path.join(OUTPUT_DIR, f'RM_ANOVA_CD_PP_boxplot_left.{ext}')
        fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f'Saved: {path}')
    plt.close(fig)


if __name__ == '__main__':
    main()
