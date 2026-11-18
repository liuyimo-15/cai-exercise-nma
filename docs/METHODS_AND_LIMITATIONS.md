# Recorded analysis specification and unresolved limitations

## Inputs, outcomes and effect direction

The six outcomes are CAIT, FAAM-ADL, sports function, anterior balance, posteromedial balance and posterolateral balance. CAIT and FAAM-ADL use MD; the remaining outcomes use Hedges g. Positive regimen effects favour the named treatment over NAT in the supplied improvement-oriented input. Input means/SDs are preserved from the fixed-value change-score dataset; this bundle does not re-extract papers or re-impute their SDs. In particular, R050 in a filename does not establish that every row's SD was calculated with a correlation of 0.50; original row notes remain authoritative.

For SMD, the pooled within-trial SD uses all selected trial arms with df=sum(n)-number_of_arms; exact Hedges correction is used. Sampling variance uses total N in its SMD term, preserving the reconstructed R050 specification. Shared multi-arm sampling covariance is retained. Between-study heterogeneity is fitted by Gaussian REML with multi-arm random-effect correlation 0.5. The latter correlation is not a baseline-follow-up correlation.

The raw CSV data are exactly the first five populated numeric input columns of the corresponding six workbook sheets. Rows containing comments/headers rather than numeric mean, SD and n are excluded. A source-cell map and executable workbook comparison are included. Selection of the NAT-connected R050 network is a separate step. Balance subset matching is implemented in `match_balance.R`, not hidden in the source data. The selected and excluded IDs are recorded in `results/tables/primary/diagnostics/*/analysed_arms.csv`, `balance_subset_matching.csv`, and the archived conditional assessment scope.

## NMA and CNMA

Primary treatment effects and 95% normal CIs are reconstructed and checked against the supplied two-decimal R050 PDF. `netmeta` independently checks effects/SE at the same tau. SIDE uses back-calculation. Study contributions use `netcontrib(method='shortestpath', study=TRUE)`. Global design-by-treatment Wald tests use the fitted primary heterogeneity; they are not automatically interchangeable with every package's inconsistency test.

CNMA terms/interaction terms are read from `data/analysis_inputs/sensitivity/spec.csv`, as displayed in the primary PDF. Each sensitivity model refits REML without reselecting terms. Non-identifiable coefficients are missing with an explicit status. This reproduces fixed specifications; the original process of selecting interactions has not been recovered. Negative interactions represent departures from additivity and cannot alone be interpreted as regimen-level harm.

## Sensitivity analyses and rankings

All arms of each omitted report/analysis ID are removed for leave-one-out fits; heterogeneity is re-estimated. High risk means `Overall == 'High risk'` for the supplied outcome-specific scenario. Unrated RPT092_306 is not silently classified as low risk; a separate high-risk-plus-unrated exclusion is tabulated. Supplemental joint deletions of linked ID pairs are separate scenarios.

P-scores are recomputed; higher ranks better and P-score is not the probability of being best. Ranking changes use the same retained node set. Missing/disconnected nodes do not receive a neutral rank. A point-sign reversal is distinct from a change in nominal CI significance and from a statistically significant reversal.

## Subgroups and network meta-regression

Subgroups use treatment-specific subgroup effects with common REML tau² and a joint global interaction test; only shared NAT-connected nodes are displayed. Primary candidate variables are duration (<6/≥6 weeks), balance tool, and activity/athlete population. Supervision is evaluated only when estimable; age restrictions are sensitivity analyses, not adult-youth interaction tests. Unknown or cross-boundary values are not arbitrarily assigned.

Univariable regressions: weeks (increment +2), baseline adult CAIT (-5), mean age (+5 years), female proportion (+0.10), planned frequency (+1 session/week). Study-level demographic covariates use selected-arm baseline demographic denominators, not necessarily final outcome n. Complete cases are refitted for each covariate. Treatment-specific interactions are the principal models; common-active-vs-NAT slopes are simplifying sensitivity models. No automatic cutpoints, splines or multivariable variable selection are introduced.

Models require at least 10 complete NAT-connected report IDs and at least 5 residual contrast df. Rank-reduced parameterizations identify estimable functions. CIs use t and global tests use F with residual contrast df; these are approximate and are not claimed as Knapp-Hartung/CR2. Core and exploratory global tests use separate BH families within model kind/scenario. Independent metafor REML fits, core-variable leave-one-out fits and report-cluster CR1 sensitivity checks are included. Predictions are confined to node-specific covariate support and estimable functions. Bubble sizes reflect inverse sampling variance; shaded bands are fitted-mean CIs, not prediction intervals.

## Optional conditional CINeMA/GRADE

This step reproduces the actual earlier, modified assessment and its caveats, rather than relabelling it as a standard website export. The unit is outcome × regimen comparison, not a single grade per outcome, CNMA coefficient or ranking. RCT evidence starts high. Within-study bias is weighted by comparison-specific study contributions (Low=1, Some=2, High=3); weighted score <1.5/1.5–2.5/>2.5 produces 0/1/2 levels. Reporting bias and indirectness are fixed at zero by user instruction. Total downgrade = within-study bias + max(imprecision, heterogeneity, incoherence), capped to four certainty categories. These are local operational rules.

The recorded precision rule is permissive: 0 if CI is wholly on one side of zero or entirely within ±delta, 2 if it crosses the clinically important threshold opposite the point estimate, otherwise 1. Heterogeneity considers a changed interpretation using PI=estimate±1.96*sqrt(SE²+tau²), without adding the same concern twice. Local SIDE P≤0.10 prompts direct/indirect clinical-direction review; severe concern requires opposite directions beyond delta. Where local SIDE is not estimable, the recorded global-P fallback is used. This is not proof of coherence in an untestable network. The code, complete reasons and threshold sensitivities are retained.

Thresholds: CAIT 3 points; FAAM-ADL 8 percentage points; sports original scale 9 points re-expressed as delta=9/12.272635666673093=0.7333388071 using the recorded reference SD; balance directions operational SMD 0.5, with 0.2/0.8 sensitivity. The latter is not a validated CAI direction-specific MCID. Individual change MCIDs used as between-group thresholds are an explicit review assumption. Sources and original assessment references are in `reference_outputs/conditional_grade_assessment.json`. No independent dual assessment or ROB-MEN execution is asserted.

## Remaining scientific issues

- RPT023/RPT184 may share a cohort/active arm. RPT092_306 and RPT141 may share a larger project and differ in tested side. Main R050 units are retained for numerical reproduction; separate cautious scenarios do not resolve trial identity.
- RPT005/RPT048 may require cluster-design variance review. Direction/tool mapping and affected versus untrained-side eligibility need adjudication. Some studies include adolescents.
- The supplied latest RoB2 ratings were previously adjusted at the user's request. They are scenario inputs, not evidence-verified revised judgments. `Original_Overall` supports comparison but is not automatically a final adjudication either.
- Original R050 inclusion/exclusion configuration, component selection process, and exact figure-producing program remain unavailable. Matching rounded outputs does not prove unique model identification or validate inclusion decisions.
- Original Excel report formatting, article extraction from journal PDFs, and the original 60-page PDF layout are archived rather than regenerated. New tables/figures cover the explicitly listed analyses.
- The local numerical rerun is distinct from a fresh internet dependency restore, Windows/Linux execution and GitHub-hosted execution; report those only after they have actually passed.
