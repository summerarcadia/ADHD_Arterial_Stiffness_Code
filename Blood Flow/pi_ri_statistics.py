#!/usr/bin/env python3
"""
Repeated-measures ANOVA across the four timepoints for carotid pulsatility
index (PI) and resistive index (RI), with Holm-Bonferroni post hoc tests
(run only when the omnibus test is significant).

Participants with a STILLS (still image only) or NONE (no usable clip)
note are excluded from both variables.

Input:  ADHD_velocities_05222026.xlsx (sheet Sheet1), same folder
Output: pi_ri_anova_results.csv, pi_ri_descriptive_stats.csv,
        pi_ri_posthoc_results.csv (only if a post hoc test was run)
"""

import os
from itertools import combinations

import numpy as np
import pandas as pd
import pingouin as pg
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(SCRIPT_DIR, 'ADHD_velocities_05222026.xlsx')
SHEET = 'Sheet1'

TP_MAP = {
    'BL': 'Baseline',
    '30min': '30min',
    '60min': '1hr',
    '5hr': '5hr'
}
TP_ORDER = ['Baseline', '30min', '1hr', '5hr']

MEASUREMENTS = [
    ('PI', 'Carotid Pulsatility Index (PI)'),
    ('RI', 'Carotid Resistive Index (RI)'),
]

EXCLUDE_FLAGS = {'STILLS', 'NONE'}


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
        'F': f_stat,
        'df_time': df_time,
        'df_error': df_error,
        'p_uncorrected': p_uncorrected,
        'mauchly_p': float(aov.loc[0, 'p_spher']),
        'gg_epsilon': float(aov.loc[0, 'eps']),
        'gg_corrected': not sphericity_met,
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


def posthoc_paired_holm(wide_df):
    rows = []
    for tp_i, tp_j in combinations(TP_ORDER, 2):
        diff = wide_df[tp_i].values - wide_df[tp_j].values
        t_stat, p_unc = stats.ttest_rel(wide_df[tp_i].values, wide_df[tp_j].values)
        sd_diff = diff.std(ddof=1)
        rows.append({
            'Comparison': f'{tp_i} vs {tp_j}',
            'Mean_Diff': diff.mean(),
            'SD_Diff': sd_diff,
            't_stat': t_stat,
            'p_uncorrected': p_unc,
            'Cohens_d': diff.mean() / sd_diff if sd_diff > 0 else 0,
        })

    posthoc = pd.DataFrame(rows)
    posthoc['p_holm_bonferroni'] = holm_bonferroni(posthoc['p_uncorrected'].tolist())
    posthoc['Significance'] = posthoc['p_holm_bonferroni'].map(sig_stars)
    return posthoc


# Load and clean
print(f"Loading data from: {DATA_FILE}")
raw = pd.read_excel(DATA_FILE, sheet_name=SHEET)
raw = raw[['ID', 'time', 'PI', 'RI', 'NOTES']].copy()
raw['timepoint_clean'] = raw['time'].map(TP_MAP)

print(f"Total observations: {len(raw)}")
print(f"Unique participants: {raw['ID'].nunique()}")


# Sample size per variable
print("\n" + "=" * 85)
print("SAMPLE SIZE PER VARIABLE")
print("=" * 85)
print(f"{'Variable':<30} {'N included':<12} {'N excluded':<12} {'Excluded (note)'}")
print("-" * 85)

flagged = raw.loc[raw['NOTES'].isin(EXCLUDE_FLAGS), 'ID'].unique()
reasons = []
for pid in sorted(flagged):
    flag = raw.loc[raw['ID'] == pid, 'NOTES'].dropna().unique()[0]
    reasons.append(f"{pid} ({flag})")
exclusion_str = '; '.join(reasons) or 'None'

raw_clean = raw[~raw['ID'].isin(flagged)].copy()

for col, title in MEASUREMENTS:
    print(f"{title:<30} {raw_clean['ID'].nunique():<12} {len(flagged):<12} {exclusion_str}")

print("=" * 85)

# other notes (e.g. PVCs) do not exclude a participant
other_notes = raw.loc[raw['ID'].isin(raw_clean['ID']) & raw['NOTES'].notna(),
                      ['ID', 'NOTES']].drop_duplicates()
if len(other_notes) > 0:
    print("\nIncluded participants with other notes:")
    for _, row in other_notes.iterrows():
        print(f"  {row['ID']}: {row['NOTES']}")


# RM-ANOVA, descriptives, and post hoc (only when the omnibus test is significant)
all_anova = []
all_posthoc = []
all_descriptive = []

print("\n\n" + "=" * 85)
print("REPEATED MEASURES ANOVA + HOLM-BONFERRONI POST HOC (PI / RI)")
print("=" * 85)

for col, title in MEASUREMENTS:
    wide = raw_clean.pivot_table(index='ID', columns='timepoint_clean', values=col)
    wide = wide[TP_ORDER].dropna()
    n = len(wide)

    res = rm_anova(wide)
    sig_label = sig_stars(res['p'])

    all_anova.append({
        'Variable': title, 'Column': col, 'N_complete': n,
        'F_stat': round(res['F'], 3), 'df_time': res['df_time'], 'df_error': res['df_error'],
        'p_value': round(res['p'], 4), 'p_uncorrected': round(res['p_uncorrected'], 4),
        'Mauchly_p': round(res['mauchly_p'], 4), 'GG_epsilon': round(res['gg_epsilon'], 3),
        'GG_corrected': res['gg_corrected'],
        'partial_eta_sq': round(res['partial_eta_sq'], 3),
        'Significance': sig_label
    })

    gg_note = ' (GG-corrected)' if res['gg_corrected'] else ''
    print(f"\n--- {title} (n={n}) ---")
    print(f"  F({res['df_time']},{res['df_error']}) = {res['F']:.3f}, "
          f"p = {res['p']:.4f}{gg_note}, eta_p^2 = {res['partial_eta_sq']:.3f} {sig_label}")

    for tp in TP_ORDER:
        vals = wide[tp]
        all_descriptive.append({
            'Variable': title, 'Timepoint': tp, 'N': n,
            'Mean': round(vals.mean(), 3), 'SD': round(vals.std(ddof=1), 3),
            'SEM': round(vals.sem(), 3), 'Median': round(vals.median(), 3),
            'Min': round(vals.min(), 3), 'Max': round(vals.max(), 3)
        })
        print(f"    {tp}: {vals.mean():.2f} +/- {vals.std(ddof=1):.2f}")

    if res['p'] >= 0.05:
        print("  Omnibus not significant, post hoc not run")
        continue

    ph = posthoc_paired_holm(wide)
    ph.insert(0, 'Variable', title)
    all_posthoc.append(ph)

    sig_ph = ph[ph['p_holm_bonferroni'] < 0.05]
    for _, row in sig_ph.iterrows():
        print(f"  {row['Comparison']}: diff={row['Mean_Diff']:.3f}+/-{row['SD_Diff']:.3f}, "
              f"d={row['Cohens_d']:.2f}, p_holm={row['p_holm_bonferroni']:.4f} {row['Significance']}")
    if sig_ph.empty:
        print("  No pairwise comparison significant after Holm-Bonferroni")


# Export
pd.DataFrame(all_anova).to_csv(os.path.join(SCRIPT_DIR, 'pi_ri_anova_results.csv'), index=False)
pd.DataFrame(all_descriptive).to_csv(os.path.join(SCRIPT_DIR, 'pi_ri_descriptive_stats.csv'), index=False)
print("\nSaved: pi_ri_anova_results.csv")
print("Saved: pi_ri_descriptive_stats.csv")

if all_posthoc:
    pd.concat(all_posthoc, ignore_index=True).to_csv(
        os.path.join(SCRIPT_DIR, 'pi_ri_posthoc_results.csv'), index=False)
    print("Saved: pi_ri_posthoc_results.csv")

print("\nDone.")
