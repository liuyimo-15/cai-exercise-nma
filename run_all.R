#!/usr/bin/env Rscript
# CSV inputs -> models -> tables -> figures -> comparison with frozen results.
cmd <- commandArgs(trailingOnly = FALSE)
filearg <- grep("^--file=", cmd, value = TRUE)
if (length(filearg)) setwd(dirname(normalizePath(sub("^--file=", "", filearg[1]))))
args <- commandArgs(trailingOnly = TRUE)
source("R/bootstrap.R")
pk <- read.csv("dependencies/package_versions.csv")
checks <- lapply(seq_len(nrow(pk)), function(i) {
  p <- pk$Package[i]
  have <- if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "MISSING"
  same <- have != "MISSING" && package_version(have) == package_version(pk$Version[i])
  data.frame(package = p, expected = pk$Version[i], installed = have, match = same)
})
checks <- do.call(rbind, checks)
if (any(checks$installed == "MISSING")) {
  print(checks[checks$installed == "MISSING", ]); stop("Run Rscript setup.R first.")
}
if (any(!checks$match) && !"--allow-version-mismatch" %in% args) {
  print(checks[!checks$match, ]); stop("Version mismatch. Restore renv.lock or explicitly use --allow-version-mismatch.")
}
for (d in c("results/logs", "results/validation", "results/tables", "results/figures/main",
            "results/figures/sensitivity/Figures", "results/figures/subgroups", "results/figures/meta_regression"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
write.csv(checks, "results/validation/package_versions_checked.csv", row.names = FALSE)
map <- c(primary = "cinema_20261112", sensitivity = "sensitivity_20261112", subgroups = "subgroups_20261113", meta_regression = "nmr_20261113")
for (g in names(map)) {
  target <- file.path("work", map[g]); dir.create(target, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(file.path("data/analysis_inputs", g), full.names = TRUE)
  stopifnot(all(file.copy(files, target, overwrite = TRUE)))
}
Rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run <- function(script, extra = character()) {
  log <- file.path("results/logs", paste0(gsub("[^A-Za-z0-9]", "_", script), ".log"))
  cat("Running", script, paste(extra, collapse = " "), "\n"); flush.console()
  t0 <- Sys.time()
  status <- system2(Rscript, c("--vanilla", shQuote(script), extra), stdout = log, stderr = log)
  if (status != 0) {
    cat(tail(readLines(log, warn = FALSE), 40), sep = "\n")
    stop("Stage failed; see ", log)
  }
  cat("  Passed in", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "seconds\n")
}
run("R/check_inputs.R")
run("work/cinema_20261112/reproduce.R")
run("work/cinema_20261112/match_balance.R")
run("work/cinema_20261112/diagnostics.R")
run("work/sensitivity_20261112/analyse.R")
run("work/subgroups_20261113/analyse.R")
run("work/nmr_20261113/analyse.R")
run("work/nmr_20261113/robust.R")
run("work/nmr_20261113/check_models.R")
if (!"--no-figures" %in% args) {
  run("R/main_figures.R")
  run("work/sensitivity_20261112/figures.R")
  run("work/subgroups_20261113/figures.R", "--clean")
  run("work/nmr_20261113/figures.R")
}
run("R/export_results.R")
run("R/verify_results.R", if ("--no-figures" %in% args) "--no-figures" else character())
if ("--with-grade" %in% args) {
  python <- Sys.getenv("PYTHON", unset = Sys.which("python3"))
  if (!nzchar(python)) stop("Python 3 required for optional grade step. See README.")
  status <- system2(python, shQuote("scripts/conditional_grade.py"), stdout = "results/logs/conditional_grade.log", stderr = "results/logs/conditional_grade.log")
  if (status != 0) stop("Conditional grade failed; see results/logs/conditional_grade.log")
  status <- system2(python, shQuote("scripts/verify_grade.py"))
  if (status != 0) stop("Conditional grade verification failed.")
}
writeLines(c(capture.output(sessionInfo()), "", paste("Completed:", format(Sys.time(), tz = "UTC")),
             "This run verifies numerical reproducibility, not trial eligibility or risk-of-bias validity."), "results/validation/sessionInfo.txt")
cat("Complete. Tables: results/tables; figures: results/figures; checks: results/validation\n")
