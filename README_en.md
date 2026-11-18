# Exercise rehabilitation for chronic ankle instability: reproducible R analyses

This repository bundles the actual six-outcome analysis inputs and previously executed R scripts, with a portable runner and frozen numerical benchmarks. It covers regimen NMA, fixed-specification CNMA, leave-one-report-out/high-risk exclusion sensitivity analyses, network subgroups and univariable network meta-regression.

**Provenance:** the original script/model that created the supplied R050 PDF was not recovered. The primary code here reconstructs its displayed treatment effects, heterogeneity and SIDE results. Matching a PDF does not recover original exclusion rationales or reproduce the complete original 60-page layout. The original fixed-value workbook and analysis PDF are included unchanged.

## Run

Use R 4.4.0, open the `.Rproj` or enter the repository root, then:

```bash
Rscript setup.R
Rscript run_all.R
```

`setup.R` restores the recorded CRAN package versions into `.Rlibrary` using `renv.lock`. Internet and, for source builds, system compilation tools may be needed. Dependency records were assembled from the installed analysis environment, not a historical renv snapshot. Version comparison normalizes R's equivalent hyphen/dot version notation.

`run_all.R` rebuilds model objects from CSV inputs, creates tables/PDFs, and checks them against previously frozen CSV benchmarks. No archived model objects are used to bypass fitting. Options: `--no-figures`, `--allow-version-mismatch`, `--with-grade`. An existing R library can be supplied via `NMA_R_LIBRARY`; exact-version checks still apply.

Optional Python steps:

```bash
python3 -m pip install -r dependencies/requirements-grade.txt
python3 scripts/check_source_workbook.py
Rscript run_all.R --with-grade
```

Set `PYTHON` to a working interpreter if required. The grade assessment is a **conditional, user-specified modified CINeMA scenario**, not a completed standard CINeMA/ROB-MEN assessment.

## Contents

- `data/analysis_inputs`: locked arm-level means, SDs, sample sizes, covariates, subgroup assignments and CNMA terms.
- `data/source_snapshots`: unchanged fixed-value workbook, original R050 reference PDF and user-scenario RoB2 PDF.
- `data/extraction`, `data/study_bibliography.csv`: numeric extraction mappings, provenance pages, missingness notes, titles and DOIs; no full-text journal PDFs.
- `work`: executable analysis scripts and intermediates. Old relative folder names are intentionally retained.
- `R`, `scripts`: portable runner helpers, plotting/data checks, optional grade assessment.
- `results`: regenerated tables, vector figures, plot data and numerical validation reports.
- `expected`: frozen original numerical outputs used only after refitting, never as fitted estimates.
- `reference_outputs`: archived review workbooks and original figures; Excel report formatting is not regenerated.
- `archive/original_R_scripts`: unmodified pre-packaging scripts, retained for comparison, not direct execution.
- `dependencies`, `renv.lock`: recorded dependency closure, R version and optional Python versions.
- `docs`: methods, dictionaries, script changes and unresolved source limitations.
- `validation/REPRODUCIBILITY_REPORT.md`: actual local validation coverage and untested environments.
- `SHA256SUMS.csv`: distribution integrity manifest.

`mean`/`sd` preserve the fixed-value analysis input (change scores and their input SDs), including any pre-existing estimated/model-derived values annotated in the workbook. They are not newly extracted baseline/post-treatment observations. Missing values remain missing. Primary raw-input CSV values can be checked against exact source workbook cells.

Primary plot pages are newly arranged networks, NMA forests, P-scores and CNMA forests. Sensitivity, subgroup and regression plots use portable PDF devices/Helvetica with platform-specific native PNG devices. Sensitivity PNGs are preview resolution; original vector PDFs are archived. Compare numeric plot data rather than PDF binary hashes.

## Interpretation and publication

Report IDs are not adjudicated independent-trial counts. Potential shared cohorts/arms, affected versus untrained side, cluster designs and direction mappings remain documented limitations. The RoB2 data contain a user-adjusted rating scenario; `Original_Overall` is preserved. The regression/subgroup analyses are post hoc exploratory analyses. Numerical reproducibility does not validate these scientific decisions.

Upload the **contents** of this extracted directory, preserving subfolders and dotfiles, to a GitHub repository. A manual Actions workflow is supplied; cloud execution is not claimed until actually run. No repository was created or published as part of preparing this bundle. No licence or author identity is assigned on the study owner's behalf; see `LICENSE-NOTICE.md`.
