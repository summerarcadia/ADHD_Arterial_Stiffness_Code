# ADHD_Arterial_Stiffness_Code

Analysis code for a study of acute arterial stiffness responses to stimulant medication in adults with ADHD.

The manuscript associated with this code is in preparation. A citation and link will be added once it is available.

## Study overview

- Within-subjects design (N = 19)
- Arterial stiffness and hemodynamic measures collected at four timepoints: baseline, 30 min, 1 hr, and 5 hr after a dose of stimulant medication
- Measures analysed include:
  - Carotid distensibility
  - Augmentation index normalised to a heart rate of 75 bpm (AIx@75)
  - Pulse wave velocity (carotid-femoral, carotid-brachial, carotid-radial)
  - Heart rate and blood pressure
  - Carotid blood flow pulsatility index (PI) and resistive index (RI)

## Statistical pipeline

1. Repeated-measures ANOVA across the four timepoints
2. Mauchly's test of sphericity. When it is significant (P < 0.05), the Greenhouse-Geisser corrected P is reported. F and degrees of freedom are reported uncorrected, and the output tables include the uncorrected P, Mauchly's P, and the Greenhouse-Geisser epsilon.
3. Holm-Bonferroni corrected post hoc comparisons (paired t-tests), run only when the omnibus test is significant
4. Effect sizes: partial eta squared (η²p) for the omnibus test and Cohen's d for pairwise comparisons

## Repository structure

The code is split into three analysis modules. Each has a statistics script and one or more graphing scripts.

```
ADHD_Arterial_Stiffness_Code/
├── NiHEM Data/
│   ├── adhd_statistics.py
│   ├── adhd_graphing_pub_version.py
│   └── adhd_graphing_aix75.py
├── CD Data/
│   ├── RM_ANOVA_CD_PP.py
│   └── plot_RM_ANOVA_CD_PP.py
├── Blood Flow/
│   ├── pi_ri_statistics.py
│   └── pi_ri_graphing_pub_version.py
├── requirements.txt
├── LICENSE
└── README.md
```

Keep the folder names exactly as shown. `CD Data/RM_ANOVA_CD_PP.py` looks for the NiHEM file in the sibling `NiHEM Data` folder.

### Module guide

| Module | Script | Purpose |
| --- | --- | --- |
| NiHEM Data | `adhd_statistics.py` | Statistics for NiHEM outputs: blood pressure, pulse wave velocity (PWV), augmentation index (AIx@75), and heart rate |
| NiHEM Data | `adhd_graphing_pub_version.py` | Publication figures for blood pressure, PWV, and heart rate |
| NiHEM Data | `adhd_graphing_aix75.py` | Publication figures for AIx@75 |
| CD Data | `RM_ANOVA_CD_PP.py` | Calculates carotid distensibility (CD) from ultrasound diameters and carotid pulse pressure (PP), then runs the repeated-measures ANOVA and Holm post hoc tests |
| CD Data | `plot_RM_ANOVA_CD_PP.py` | Publication figure for carotid distensibility (left carotid) |
| Blood Flow | `pi_ri_statistics.py` | Statistics for pulsatility index (PI) and resistive index (RI) |
| Blood Flow | `pi_ri_graphing_pub_version.py` | Publication figures for PI and RI |

## Requirements

- Python 3.13
- pandas
- numpy
- scipy
- pingouin
- matplotlib
- openpyxl

Install the tested versions with:

```
pip install -r requirements.txt
```

## Data

Participant data are **not** included in this repository.

- The data were collected under ethics approval that does not cover public release.
- To run the code, place your own Excel files in the module folders described below.
- Each script finds its input file in its own folder, so the file names, sheet names, and column names must match exactly.

### Timepoint labels

Each file uses its own labels for the same four timepoints.

| Timepoint | NiHEM file | Carotid distensibility file | Blood flow file |
| --- | --- | --- | --- |
| Baseline | `Baseline` | `Baseline` | `BL` |
| 30 min | `Condition 1` | `30min` | `30min` |
| 1 hr | `Condition 2` | `1hr` | `60min` |
| 5 hr | `Condition 3` | `5hr` | `5hr` |

### 1. NiHEM data (`NiHEM Data/NIHem_Data_Cleaned_3.xlsx`)

- One sheet named `tonometry_cleaned`
- One row per participant per timepoint (19 participants x 4 timepoints = 76 rows)
- Used by all three NiHEM scripts and by both carotid distensibility scripts

Columns the code reads:

| Column | Description |
| --- | --- |
| `sequence_no` | Participant number (for example `00001`). Becomes participant ID `P01` |
| `timepoint` | Timepoint label (see table above) |
| `systolic.blood.pressure.mmhg` | Systolic blood pressure |
| `diastolic.blood.pressure.mmhg` | Diastolic blood pressure |
| `mean.arterial.pressure.mmhg` | Mean arterial pressure |
| `carotid.femoral.pwv.cm.s` | Carotid-femoral PWV |
| `carotid.brachial.pwv.cm.s` | Carotid-brachial PWV |
| `carotid.radial.pwv.cm.s` | Carotid-radial PWV |
| `carotid.ai.percent` | Augmentation index (AIx), used to calculate AIx@75 |
| `carotid.waveform.rr.interval.ms` | Carotid R-R interval. Heart rate is calculated as 60000 divided by this value |
| `carotid.pfoot.mmhg` | Carotid foot pressure (used for pulse pressure) |
| `carotid.p.at.c0tmax.mmhg` | Carotid maximum pressure (used for pulse pressure) |

Other notes:

- The values -9999, -9998, and -9997 are treated as missing.
- Other columns exported by the device are ignored by the code.
- Columns that identify people (for example initials and study dates) are not used. Remove them from any file you share.
- A participant is left out of a variable only if they are missing a timepoint for that variable.

### 2. Carotid distensibility data (`CD Data/adhd_carotid_study_pp_corrected.xlsx`)

This workbook holds the left carotid ultrasound diameter measurements. It has these sheets:

| Sheet | Content |
| --- | --- |
| `Instructions` | Data entry protocol |
| `Baseline_Left`, `30min_Left`, `1hr_Left`, `5hr_Left` | Raw diameter measurements, one sheet per timepoint |
| `Analysis` | Summary per participant and timepoint |
| `CD Left` | Diameters per participant and timepoint, used for the distensibility calculation |
| `Long_Format_Python` | Long format copy of the measurements |

Raw measurement sheets (for example `Baseline_Left`):

- Columns: `Participant`, `Image`, `Phase`, then `Cycle 1` to `Cycle 10` with measurements `M1`, `M2`, `M3` under each cycle, then `Average`
- Per participant and per timepoint: 3 images, 10 cardiac cycles per image, 3 measurements per cycle in systole and 3 in diastole
- Each participant has 6 rows (systole and diastole for each of the 3 images)
- Participant IDs look like `P001`

`CD Left` sheet (the only sheet the code reads):

- Notes are in rows 1 to 6, column headers are in row 7, and data start at row 8
- Columns A to D are read by the code:

| Column | Name | Description |
| --- | --- | --- |
| A | `Participant_ID` | For example `P001`. The number is matched to `sequence_no` in the NiHEM file |
| B | `Timepoint` | `Baseline`, `30min`, `1hr`, or `5hr` |
| C | `Ds (cm)` | Mean systolic diameter |
| D | `Dd (cm)` | Mean diastolic diameter |

- Column J must hold a number for a row to be kept (the value itself is not used)

Carotid distensibility is calculated in the script:

```
delta_D = Ds - Dd
DC = (2 * delta_D * Dd + delta_D^2) / (PP * Dd^2) * 1000
```

- Units are 10^-3 mmHg^-1
- PP is the carotid pulse pressure from the NiHEM file: `carotid.p.at.c0tmax.mmhg` minus `carotid.pfoot.mmhg`
- Pulse pressure is taken from carotid tonometry, not from the brachial cuff

### 3. Blood flow data (`Blood Flow/ADHD_velocities_05222026.xlsx`)

- One sheet named `Sheet1`, one row per participant per timepoint
- Columns:

| Column | Description |
| --- | --- |
| `ID` | Participant ID, for example `ADHD01` |
| `time` | Timepoint label (see table above) |
| `PI` | Pulsatility index |
| `RI` | Resistive index |
| `NOTES` | Optional flag. `STILLS` (still image only) or `NONE` (no usable clip) excludes that participant. Other notes are kept and printed for review |

Other columns (for example `Vmax` and `Vmin`) are ignored by the code.

## Usage

Run each command from the repository root. The scripts can be run in any order: each graphing script calculates the statistics it displays.

NiHEM measures:

```
python "NiHEM Data/adhd_statistics.py"
python "NiHEM Data/adhd_graphing_pub_version.py"
python "NiHEM Data/adhd_graphing_aix75.py"
```

Carotid distensibility (needs the NiHEM file as well):

```
python "CD Data/RM_ANOVA_CD_PP.py"
python "CD Data/plot_RM_ANOVA_CD_PP.py"
```

Blood flow (PI and RI):

```
python "Blood Flow/pi_ri_statistics.py"
python "Blood Flow/pi_ri_graphing_pub_version.py"
```

### Outputs

All outputs are written next to the script that makes them.

| Script | Output |
| --- | --- |
| `adhd_statistics.py` | `anova_results.csv`, `descriptive_stats.csv`, `posthoc_results.csv` |
| `adhd_graphing_pub_version.py` | `Graphs_pub/` (one PNG per blood pressure, PWV, and heart rate measure, and `PWV_combined.png`) |
| `adhd_graphing_aix75.py` | `Graphs_aix75/AIx_75.png` |
| `RM_ANOVA_CD_PP.py` | `RM_ANOVA_CD_PP_output.txt` (full printed log) and `RM_ANOVA_CD_PP_*.csv` (long data, descriptives, ANOVA table, summary, post hoc) |
| `plot_RM_ANOVA_CD_PP.py` | `Graphs_pub/RM_ANOVA_CD_PP_boxplot_left.png` and `.pdf` |
| `pi_ri_statistics.py` | `pi_ri_anova_results.csv`, `pi_ri_descriptive_stats.csv`, `pi_ri_posthoc_results.csv` |
| `pi_ri_graphing_pub_version.py` | `Graphs_pub/` (`PI.png`, `RI.png`, `PI_RI_combined.png`) |

Post hoc files are written only when at least one omnibus test is significant. Rerunning a script overwrites its previous outputs.

## Notes on methods

- Carotid distensibility is calculated from carotid tonometry pressure.
- AIx is normalised to a heart rate of 75 bpm using a linear correction: AIx@75 = AIx - 0.39 x (75 - HR). The coefficient comes from Wilkinson et al. (2000, J Physiol). See the header of `adhd_statistics.py` for details.
- Effect size values in the output files may be signed. Direction is reported in the manuscript text, so magnitudes are reported unsigned.

## Author

Hayley So  
University of Waterloo  
Supervisor: Jason Au

## License

Released under the [MIT License](LICENSE).

## Contact

h2so@uwaterloo.ca
