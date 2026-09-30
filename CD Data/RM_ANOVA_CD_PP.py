"""
Carotid distensibility coefficient (DC), left carotid, across the four
timepoints: one-way repeated-measures ANOVA with Holm-Bonferroni post hoc
tests (run only when the omnibus test is significant).

DC is calculated from the ultrasound diameters and carotid pulse pressure
(PP) from tonometry:

    PP = carotid.p.at.c0tmax.mmhg - carotid.pfoot.mmhg
    DC = (2 * dD * Dd + dD^2) / (PP * Dd^2) * 1000,   dD = Ds - Dd
    units: 10^-3 mmHg^-1

Only participants with all four timepoints are analysed. P is
Greenhouse-Geisser corrected when Mauchly's test is significant.

Input:  adhd_carotid_study_pp_corrected.xlsx (sheet CD Left), same folder
        ../NiHEM Data/NIHem_Data_Cleaned_3.xlsx (sheet tonometry_cleaned)
Output: RM_ANOVA_CD_PP_output.txt, RM_ANOVA_CD_PP_data_long.csv,
        RM_ANOVA_CD_PP_descriptives.csv, RM_ANOVA_CD_PP_anova_table.csv,
        RM_ANOVA_CD_PP_summary.csv, RM_ANOVA_CD_PP_posthoc_holm.csv
"""

import itertools
import os
import sys

import numpy as np
import openpyxl
import pandas as pd
import pingouin as pg
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
EXCEL_FILE = os.path.join(HERE, 'adhd_carotid_study_pp_corrected.xlsx')
NIHEM_FILE = os.path.join(os.path.dirname(HERE), 'NiHEM Data', 'NIHem_Data_Cleaned_3.xlsx')
NIHEM_SHEET = 'tonometry_cleaned'
PREFIX = 'RM_ANOVA_CD_PP'

TIMEPOINTS = ['Baseline', '30min', '1hr', '5hr']

NIHEM_TO_CD_TIMEPOINT = {
    'Baseline': 'Baseline',
    'Condition 1': '30min',
    'Condition 2': '1hr',
    'Condition 3': '5hr',
}

MISSING_CODES = [-9999, -9998, -9997]


def out_path(suffix):
    return os.path.join(HERE, f'{PREFIX}_{suffix}')


class Tee:
    """Write printed output to the screen and to a log file."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for s in self.streams:
            s.write(text)

    def flush(self):
        for s in self.streams:
            s.flush()


def participant_number(participant_id):
    digits = ''.join(ch for ch in str(participant_id) if ch.isdigit())
    return int(digits) if digits else None


def calc_dc(Ds, Dd, PP):
    delta_d = Ds - Dd
    return (2 * delta_d * Dd + delta_d ** 2) / (PP * Dd ** 2) * 1000


def build_carotid_pp_lookup(nihem_file_path):
    """Carotid PP keyed by (participant number, timepoint)."""
    df = pd.read_excel(nihem_file_path, sheet_name=NIHEM_SHEET)
    df = df.replace(MISSING_CODES, np.nan)
    df['carotid_PP'] = df['carotid.p.at.c0tmax.mmhg'] - df['carotid.pfoot.mmhg']

    lookup = {}
    for _, row in df.iterrows():
        cd_tp = NIHEM_TO_CD_TIMEPOINT.get(row['timepoint'])
        if cd_tp is not None:
            lookup[(int(row['sequence_no']), cd_tp)] = row['carotid_PP']

    print(f"  Carotid PP: {len(lookup)} (participant, timepoint) entries "
          f"from {df['sequence_no'].nunique()} participants")
    return lookup


def extract_dc_values(excel_file_path, carotid_pp_lookup):
    """Read Ds and Dd from the CD Left sheet and calculate DC."""
    wb = openpyxl.load_workbook(excel_file_path, data_only=True, read_only=True)
    ws = wb['CD Left']

    rows = []
    n_excluded = 0
    # headers are in row 7; rows without a number in column J are skipped
    for r in ws.iter_rows(min_row=8, max_col=10, values_only=True):
        participant, timepoint, ds, dd = r[0], r[1], r[2], r[3]
        if not participant or not timepoint or not isinstance(r[9], (int, float)):
            continue

        pnum = participant_number(participant)
        pp = carotid_pp_lookup.get((pnum, timepoint))
        if pp is None or pd.isna(pp) or pp <= 0:
            print(f"    No valid carotid PP, row excluded: {participant} / {timepoint}")
            n_excluded += 1
            continue

        rows.append({
            'Participant_ID': participant,
            'Timepoint': timepoint,
            'Ds': ds,
            'Dd': dd,
            'Carotid_PP': pp,
            'DC': calc_dc(ds, dd, pp),
        })

    wb.close()
    df = pd.DataFrame(rows)
    print(f"  Rows used: {len(df)}   excluded: {n_excluded}")
    print(f"  Participants: {df['Participant_ID'].nunique()}")
    return df


def complete_cases(df):
    """Wide (participant x timepoint) and long tables for complete cases."""
    wide = df.pivot_table(index='Participant_ID', columns='Timepoint', values='DC')[TIMEPOINTS]
    wide = wide.dropna()
    long_df = df[df['Participant_ID'].isin(wide.index) & df['Timepoint'].isin(TIMEPOINTS)].copy()
    long_df['Timepoint'] = pd.Categorical(long_df['Timepoint'], categories=TIMEPOINTS, ordered=True)
    return wide, long_df.sort_values(['Participant_ID', 'Timepoint'])


def rm_anova(long_df):
    """pingouin RM-ANOVA table plus a summary dict.

    P is Greenhouse-Geisser corrected when Mauchly's test is significant.
    """
    aov = pg.rm_anova(data=long_df, dv='DC', within='Timepoint', subject='Participant_ID',
                      detailed=True, effsize='np2', correction=True)
    # older pingouin versions use hyphens in column names (p-GG-corr)
    aov.columns = [c.replace('-', '_') for c in aov.columns]

    tp = aov[aov['Source'] == 'Timepoint'].iloc[0]
    err = aov[aov['Source'] == 'Error'].iloc[0]
    sphericity_met = float(tp['p_spher']) >= 0.05
    summary = {
        'F': float(tp['F']),
        'df1': int(tp['DF']),
        'df2': int(err['DF']),
        'P_uncorrected': float(tp['p_unc']),
        'partial_eta_sq': float(tp['np2']),
        'Mauchly_W': float(tp['W_spher']),
        'Mauchly_P': float(tp['p_spher']),
        'GG_epsilon': float(tp['eps']),
        'P_GG_corrected': float(tp['p_GG_corr']),
        'GG_corrected': not sphericity_met,
    }
    summary['P'] = summary['P_uncorrected'] if sphericity_met else summary['P_GG_corrected']
    return aov, summary


def holm_bonferroni(p_values):
    """Holm-adjusted P values (step-down, monotone)."""
    p = np.asarray(p_values, dtype=float)
    m = len(p)
    adjusted = np.empty(m)
    running_max = 0.0
    for rank, idx in enumerate(np.argsort(p)):
        running_max = max(running_max, (m - rank) * p[idx])
        adjusted[idx] = min(1.0, running_max)
    return adjusted


def posthoc_paired_holm(wide):
    rows = []
    for tp1, tp2 in itertools.combinations(TIMEPOINTS, 2):
        a = wide[tp1].values
        b = wide[tp2].values
        t, p = stats.ttest_rel(a, b)
        d = a - b
        rows.append({
            'Comparison': f'{tp1} vs {tp2}',
            'Mean_1': a.mean(), 'Mean_2': b.mean(),
            'Difference': d.mean(),
            't': t, 'df': len(a) - 1,
            'p_uncorrected': p,
            'cohens_dz': d.mean() / d.std(ddof=1),
        })
    ph = pd.DataFrame(rows)
    ph['p_holm_adjusted'] = holm_bonferroni(ph['p_uncorrected'].values)
    ph['Holm_significant'] = np.where(ph['p_holm_adjusted'] < 0.05, 'Yes', 'No')
    return ph


def run():
    line = "=" * 70
    print(line)
    print("Carotid distensibility (left): one-way repeated-measures ANOVA")
    print(line)

    print("\nCalculating DC from carotid PP")
    lookup = build_carotid_pp_lookup(NIHEM_FILE)
    df = extract_dc_values(EXCEL_FILE, lookup)

    wide, df_rm = complete_cases(df)
    n_all = df['Participant_ID'].nunique()
    print(f"\nComplete cases (all 4 timepoints): n = {len(wide)} (dropped {n_all - len(wide)})")
    df_rm.to_csv(out_path('data_long.csv'), index=False)

    print("\nDescriptive statistics (DC, 10^-3 mmHg^-1)")
    desc = df_rm.groupby('Timepoint', observed=True)['DC'].agg(
        n='count', mean='mean', SD='std',
        SE=lambda x: x.std() / np.sqrt(len(x)),
    ).reindex(TIMEPOINTS)
    desc['CI95_lower'] = desc['mean'] - 1.96 * desc['SE']
    desc['CI95_upper'] = desc['mean'] + 1.96 * desc['SE']
    print(desc.round(4).to_string())
    desc.to_csv(out_path('descriptives.csv'))

    print("\nRepeated-measures ANOVA (within = Timepoint, subject = Participant_ID)")
    aov, s = rm_anova(df_rm)
    print(aov.to_string())
    aov.to_csv(out_path('anova_table.csv'), index=False)

    print(f"\n  F({s['df1']}, {s['df2']}) = {s['F']:.3f}")
    print(f"  Uncorrected P = {s['P_uncorrected']:.4f}")
    print(f"  Partial eta-squared = {s['partial_eta_sq']:.3f}")
    print(f"  Mauchly's W = {s['Mauchly_W']:.3f}, P = {s['Mauchly_P']:.3f}")
    print(f"  Greenhouse-Geisser epsilon = {s['GG_epsilon']:.3f}, "
          f"GG-corrected P = {s['P_GG_corrected']:.4f}")
    print(f"  Reported P = {s['P']:.4f} "
          f"({'GG-corrected' if s['GG_corrected'] else 'sphericity met, uncorrected'})")

    summary = pd.DataFrame([{'Outcome': 'Carotid DC (carotid PP), left', 'n': len(wide), **s}])
    summary.to_csv(out_path('summary.csv'), index=False)

    print("\nPost hoc")
    if s['P'] < 0.05:
        print("  Paired t-tests with Holm-Bonferroni correction")
        ph = posthoc_paired_holm(wide)
        print(ph.round(4).to_string(index=False))
        ph.to_csv(out_path('posthoc_holm.csv'), index=False)
    else:
        print("  Omnibus not significant, post hoc not run")

    print("\n" + line)
    print(f"Results saved to {HERE} (files start with {PREFIX}_)")
    print(line)


if __name__ == '__main__':
    with open(out_path('output.txt'), 'w', encoding='utf-8') as log:
        original_stdout = sys.stdout
        sys.stdout = Tee(original_stdout, log)
        try:
            run()
        finally:
            sys.stdout = original_stdout
