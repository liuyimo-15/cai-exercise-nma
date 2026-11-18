# Data dictionary and joins

All CSV files use UTF-8 with a header row. Blank/NA means unavailable, not zero. Identifiers and categorical strings are case sensitive. Input rows remain in their original order to preserve numerical benchmark ordering.

## Keys and provenance

| Field | Meaning |
|---|---|
| `outcome` / `key` | CAIT, FAAM_ADL, Sports, Anterior, Posteromedial, Posterolateral in executable inputs; extraction tables may use the corresponding Chinese names |
| `id` / `analysis_id` | Model unit. Usually RPTnnn; RPT092_306 combines reports and must not be interpreted as two independent trials |
| `report_id` / `report_ids` | Original paper identifiers. Multiple report IDs are separated by semicolons |
| `t` / `treatment` / `node` | Regimen code, potentially several source component codes joined by `+`; NAT is the reference node |
| `original_arm` | Arm label in the source extraction, mapped to a regimen node |
| `r050` | Whether the row/report is used in the reconstructed R050 network |
| `source_pages` / `page` | Retained source-page locator; PDF-page numbers and printed page numbers may differ, as noted in the extraction |
| `source_file`, `title`, `doi` | Locate the report through `data/study_bibliography.csv`; journal full texts are not bundled |

Primary arm rows are keyed by `(outcome,id,t)`. Demographic/extraction joins may use `(outcome,analysis_id,node,report_id,original_arm)`. Covariates are keyed by `(outcome,id)`, not ID alone. A single trial can appear in multiple outcome networks. Do not add participants across outcomes to obtain a unique-person total.

Source component codes are BAL, FAM, STR, PAG, CST and GFM, preserving the fixed workbook classification. `data/extraction/component_matrix_from_workbook.csv` records their exact 0/1 incidence for each source regimen/outcome. A `:` term in the CNMA specification is a product interaction of the named component indicators, not a new regimen. This repository does not reinterpret or invent missing clinical definitions of these codes.

## Arm-level analysis inputs

| Field | Type and unit | Use |
|---|---|---|
| `mean` | Continuous; original change-score scale | Actual fixed analysis input. CAIT points, FAAM-ADL percentage-point scale as retained; sports/balance preserve source scale before SMD standardization |
| `sd` | Positive continuous; same scale as mean | Actual input SD, including source-reported, reconstructed or model-derived SDs already present in the workbook. Row-specific source notes are retained in the source cell map |
| `n` | Integer >1 | Arm-specific analysis denominator; not automatically the randomized n or baseline demographic n |
| `sm` | MD or SMD | Meta-analysis effect metric; executable SMD uses exact Hedges correction |

The raw six input CSVs contain 356 arm rows before connectivity/subset selection. `meta_regression/arms.csv` and `subgroups/arms.csv` contain the 325 selected R050 arm rows. These are arm/outcome records, not counts of independent trials. The SDs are used as supplied; baseline/post SDs and a within-person correlation are not additional inputs to this runner.

## Continuous covariates

| Executable field | Source extraction field | Unit / model increment |
|---|---|---|
| `weeks` | `weeks` | Fixed planned duration in weeks; slope per +2 weeks |
| `cait` | `cait_baseline_adult` | Selected-arm demographic-n-weighted adult CAIT baseline mean; slope per -5 points |
| `age` | `age_mean_years` | Selected-arm demographic-n-weighted mean age in years; slope per +5 years |
| `female` | `female_prop` | Female count / corresponding demographic denominator, range 0–1; slope per +0.10 |
| `frequency` | `sessions_per_week` | Fixed planned sessions/week; slope per +1 session/week |
| `increment` | Locked scaling in `variables.csv` | Includes the negative CAIT increment; preserve its sign |
| `core` | 1=weeks/cait, 0=age/female/frequency | Defines the separate BH families used in each model kind/scenario |

Ranges, uncertain values, proxy timing estimates and pediatric CAIT scores are not substituted for verified main regression values. The extraction tables retain separate reasons/conditional fields. The baseline severity covariate can be mathematically related to the change outcome; study-level regression remains observational/ecological.

## Subgroup fields

| Field | Executable levels |
|---|---|
| `duration` | Short_lt6, Long_ge6, Unknown |
| `tool` | SEBT, YBT, Other_NA / unknown as recorded |
| `population_group` | Recreational_active, Athlete, Unknown |
| `supervision_group` | Onsite, Home_remote, Hybrid, Unknown |
| `age_scope` | Adult_verified, Includes_under18, Age_unconfirmed |
| `frequency_group` | Less3, Atleast3, Unknown |
| `sex_composition` | Male_only, Female_only, Mixed, Unknown |

Category definitions and per-report reasons are retained in the extraction tables and `subgroup_definitions.json`. RPT005's cross-boundary duration is not assigned using the title alone. `Adult_verified` requires evidence supporting adult-only recruitment; a study mean age ≥18 alone is insufficient. Joint subgroups require at least five connected report IDs in each stratum, shared active nodes and sufficient residual df; an unsupported candidate has a status, not an invented estimate.

## Risk, specifications and outputs

| Field | Meaning |
|---|---|
| `D1`–`D5`, `Overall` | Supplied outcome-specific user-scenario RoB2: Low risk, Some concerns, High risk |
| `Original_Overall` | Previously retained original Overall judgment, supplied for comparison; not silently substituted in the high-risk sensitivity |
| `terms` | Semicolon-separated locked CNMA main/interaction terms |
| `estimate`, `se`, `lower`, `upper` | Estimated contrast/coefficient, standard error, nominal 95% CI |
| `pi_lower`, `pi_upper` | Primary Gaussian prediction interval, incorporating the same tau² but not tau²-estimation uncertainty |
| `tau2` | Between-study variance on the applicable MD/SMD scale |
| `Pscore`, `rank` | Regimen ranking statistic and rank, larger P-score better; not a best-treatment probability |
| `estimable`, `status` | Explicit availability/identifiability; missing estimates are not zero |
| `scenario`, `type`, `removed` | Primary / exclusion / LOO model and omitted IDs |
| `point_reversal` | Estimated sign changed relative to its primary comparison |
| `significance_changed` | CI zero-crossing classification changed relative to primary |
| `significant_reversal` | Significant benefit changed to significant harm or vice versa |
| `p_global`, `p_interaction`, `q_BH` | Global unadjusted interaction P and its recorded BH-adjusted q; not a within-subgroup P-value comparison |
| `kind` | Treatment_specific interaction model or Common_active simplifying slope model |
| `scenario=R050` | Main retained report scope; Cohort_cautious is a separate restriction |
| `x0` | Median complete-case covariate centering value, model-specific |
| `residual_df` | Contrast-level residual df used for main t/F inference |
| `cluster_residual_df` | Report-cluster k minus fitted rank, for CR1 sensitivity eligibility |
| `inside_node_range` | Whether prediction is supported by the observed covariate range of the node |
| `grade` | Conditional comparison-specific modified-CINeMA certainty, not a rating of a component or P-score |

`docs/input_file_inventory.csv` enumerates every actual input file and its columns. `results/plot_data/` provides exact plotting numbers; original paper figures were not redigitized during packaging.
