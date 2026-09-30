#!/usr/bin/env python3
"""
Repeated-measures ANOVA across the four timepoints for the NiHEM measures
(blood pressure, PWV, AIx@75, heart rate), with Holm-Bonferroni post hoc
tests (run only when the omnibus test is significant).

Augmentation index is normalised to a heart rate of 75 bpm using the linear
correction from Wilkinson et al. (2000), "The influence of heart rate on
augmentation index and central arterial pressure in humans",
J Physiol 525(1):263-270:

    AIx@75 = AIx - 0.39 x (75 - HR),   HR = 60000 / carotid R-R interval (ms)

Participants are excluded from a variable only if they are missing a
timepoint for that variable.

Input:  NIHem_Data_Cleaned_3.xlsx (sheet tonometry_cleaned), same folder
Output: anova_results.csv, posthoc_results.csv,
        descriptive_stats.csv
"""

import os
from itertools import combinations

import numpy as np
import pandas as pd
import pingouin as pg
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(SCRIPT_DIR, 'NIHem_Data_Cleaned_3.xlsx')
SHEET = 'tonometry_cleaned'

MISSING_CODES = [-9999, -9998, -9997]

TP_MAP = {
    'Baseline': 'Baseline',
    'Condition 1': '30min',
    'Condition 2': '1hr',
    'Condition 3': '5hr'
}
TP_ORDER = ['Baseline', '30min', '1hr', '5hr']

MEASUREMENTS = [
    ('systolic.blood.pressure.mmhg',   'Systolic Blood Pressure (SBP)'),
    ('diastolic.blood.pressure.mmhg',  'Diastolic Blood Pressure (DBP)'),
    ('mean.arterial.pressure.mmhg',    'Mean Arterial Pressure (MAP)'),
    ('carotid.femoral.pwv.cm.s',       'Carotid-Femoral PWV'),
    ('carotid.brachial.pwv.cm.s',      'Carotid-Brachial PWV'),
    ('carotid.radial.pwv.cm.s',        'Carotid-Radial PWV'),
    ('AIx_75',                         'HR-Corrected Augmentation Index (AIx@75)'),
    ('hr_carotid',                     'Heart Rate'),
]


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


def complete_wide(col):
    wide = raw.pivot_table(index='participant', columns='timepoint_clean', values=col)
    return wide[TP_ORDER].dropna()


# Load and clean
print(f"Loading data from: {DATA_FILE}")
raw = pd.read_excel(DATA_FILE, sheet_name=SHEET)
raw = raw.replace(MISSING_CODES, np.nan)

raw['timepoint_clean'] = raw['timepoint'].map(TP_MAP)
raw['participant'] = 'P' + raw['sequence_no'].astype(str).str.zfill(2)
raw['hr_carotid'] = 60000 / raw['carotid.waveform.rr.interval.ms']
raw['AIx_75'] = raw['carotid.ai.percent'] - 0.39 * (75 - raw['hr_carotid'])

print(f"Total observations: {len(raw)}")
print(f"Unique participants: {raw['participant'].nunique()}")


# Sample size per variable
print("\n" + "=" * 85)
print("SAMPLE SIZE PER VARIABLE")
print("=" * 85)
print(f"{'Variable':<40} {'N included':<12} {'N excluded':<12} {'Excluded (missing timepoints)'}")
print("-" * 85)

for col, title in MEASUREMENTS:
    wide_all = raw.pivot_table(index='participant', columns='timepoint_clean', values=col)[TP_ORDER]
    complete = wide_all.dropna()
    incomplete = wide_all[wide_all.isna().any(axis=1)]

    reasons = []
    for pid in sorted(incomplete.index):
        missing = [tp for tp in TP_ORDER if pd.isna(wide_all.loc[pid, tp])]
        reasons.append(f"{pid} missing {', '.join(missing)}")
    print(f"{title:<40} {len(complete):<12} {len(incomplete):<12} {'; '.join(reasons) or 'None'}")

print("=" * 85)


# RM-ANOVA, descriptives, and post hoc (only when the omnibus test is significant)
all_anova = []
all_posthoc = []
all_descriptive = []

print("\n\n" + "=" * 85)
print("REPEATED MEASURES ANOVA + HOLM-BONFERRONI POST HOC")
print("=" * 85)

for col, title in MEASUREMENTS:
    wide = complete_wide(col)
    n = len(wide)

    if n < 3:
        print(f"\n  Skipped {title}: only {n} complete cases")
        continue

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
            'Mean': round(vals.mean(), 1), 'SD': round(vals.std(ddof=1), 1),
            'SEM': round(vals.sem(), 1), 'Median': round(vals.median(), 1),
            'Min': round(vals.min(), 1), 'Max': round(vals.max(), 1)
        })

    if res['p'] >= 0.05:
        print("  Omnibus not significant, post hoc not run")
        continue

    ph = posthoc_paired_holm(wide)
    ph.insert(0, 'Variable', title)
    all_posthoc.append(ph)

    sig_ph = ph[ph['p_holm_bonferroni'] < 0.05]
    for _, row in sig_ph.iterrows():
        print(f"  {row['Comparison']}: diff={row['Mean_Diff']:.1f}+/-{row['SD_Diff']:.1f}, "
              f"p_holm={row['p_holm_bonferroni']:.4f} {row['Significance']}")
    if sig_ph.empty:
        print("  No pairwise comparison significant after Holm-Bonferroni")


# Change from baseline at 5 hr
print("\n\n" + "=" * 85)
print("CHANGE FROM BASELINE AT 5HR")
print("=" * 85)

for col, title in MEASUREMENTS:
    wide = complete_wide(col)
    if len(wide) < 3:
        continue
    change = wide['5hr'] - wide['Baseline']
    print(f"{title}: {change.mean():.0f} +/- {change.std(ddof=1):.0f}")
    print(f"  5hr range: {wide['5hr'].min():.0f} - {wide['5hr'].max():.0f}")


# Export
pd.DataFrame(all_anova).to_csv(os.path.join(SCRIPT_DIR, 'anova_results.csv'), index=False)
pd.DataFrame(all_descriptive).to_csv(os.path.join(SCRIPT_DIR, 'descriptive_stats.csv'), index=False)
print("\nSaved: anova_results.csv")
print("Saved: descriptive_stats.csv")

if all_posthoc:
    pd.concat(all_posthoc, ignore_index=True).to_csv(
        os.path.join(SCRIPT_DIR, 'posthoc_results.csv'), index=False)
    print("Saved: posthoc_results.csv")

print("\nDone.")
