source("R/bootstrap.R")
manifest <- read.csv("validation/expected_files.csv")
if ("--no-figures" %in% commandArgs(trailingOnly = TRUE))
  manifest <- manifest[!basename(manifest$path) %in% c("bubble_points.csv", "curve_values.csv"), ]
rows <- list()
for (i in seq_len(nrow(manifest))) {
  path <- manifest$path[i]
  # PNG/PDF binary checksums are not meaningful numerical reproducibility tests.
  ref <- read.csv(file.path("expected", path), check.names = FALSE, na.strings = c("", "NA"))
  got <- read.csv(path, check.names = FALSE, na.strings = c("", "NA"))
  if (!identical(dim(ref), dim(got)) || !identical(names(ref), names(got))) stop("Schema mismatch: ", path)
  err <- 0; numeric_cells <- 0; ok <- TRUE
  for (j in seq_along(ref)) {
    a <- ref[[j]]; b <- got[[j]]
    if (!identical(is.na(a), is.na(b))) { ok <- FALSE; break }
    k <- !is.na(a)
    if (is.numeric(a) && is.numeric(b)) {
      finite <- k & is.finite(a) & is.finite(b)
      if (any(k & !finite) && !identical(a[k & !finite], b[k & !finite])) ok <- FALSE
      if (any(finite)) {
        delta <- abs(a[finite] - b[finite]); err <- max(err, delta)
        ok <- ok && all(delta <= manifest$absolute_tolerance[i] + manifest$relative_tolerance[i] * abs(a[finite]))
        numeric_cells <- numeric_cells + sum(finite)
      }
    } else ok <- ok && identical(as.character(a[k]), as.character(b[k]))
  }
  rows[[path]] <- data.frame(path = path, rows = nrow(got), numeric_cells = numeric_cells, max_abs_error = err, passed = ok)
}
report <- do.call(rbind, rows); rownames(report) <- NULL
write.csv(report, "results/validation/numeric_reproduction.csv", row.names = FALSE)
if (!all(report$passed)) { print(report[!report$passed, ]); stop("Reference result comparison failed.") }
# Check selected primary network estimates against the displayed R050 values.
pdfref <- read.csv("work/cinema_20261112/pdf_reference.csv")
est <- read.csv("work/cinema_20261112/all_estimates.csv")
j <- match(paste(pdfref$outcome, pdfref$treatment, "NAT"), paste(est$outcome, est$treatment, est$reference))
stopifnot(!anyNA(j))
delta <- abs(as.matrix(pdfref[, c("estimate", "lower", "upper")]) - as.matrix(est[j, c("estimate", "lower", "upper")]))
stopifnot(max(delta) < .00501)
cnma <- read.csv("work/sensitivity_20261112/CNMA_pdf_validation.csv")
write.csv(data.frame(kind = c("NMA", "CNMA"), coefficients = c(nrow(pdfref), cnma$reference_rows),
                     max_display_error = c(max(delta), cnma$max_difference), tolerance = .00501, passed = TRUE),
          "results/validation/R050_display_agreement.csv", row.names = FALSE)
cat(nrow(report), "reference CSV tables passed;", sum(report$numeric_cells), "numeric cells checked.\n")
