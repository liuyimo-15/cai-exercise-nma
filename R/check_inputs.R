source("R/bootstrap.R")
if (file.exists("validation/input_md5.csv")) {
  integrity <- read.csv("validation/input_md5.csv")
  got <- unname(tools::md5sum(integrity$path))
  stopifnot(!anyNA(got), identical(got, integrity$md5))
}
keys <- c("CAIT", "FAAM_ADL", "Sports", "Anterior", "Posteromedial", "Posterolateral")
a <- read.csv("work/nmr_20261113/arms.csv")
s <- read.csv("work/subgroups_20261113/arms.csv")
makekey <- function(d) paste(d$outcome, d$id, d$t, sep = "|")
stopifnot(!anyDuplicated(makekey(a)), !anyDuplicated(makekey(s)))
stopifnot(setequal(makekey(a), makekey(s)))
i <- match(makekey(a), makekey(s))
stopifnot(max(abs(as.matrix(a[, c("n", "mean", "sd")]) - as.matrix(s[i, c("n", "mean", "sd")]))) < 1e-12)
stopifnot(all(is.finite(a$mean)), all(a$sd > 0), all(a$n > 1), all(a$n == as.integer(a$n)))
cv <- read.csv("work/nmr_20261113/covariates.csv")
stopifnot(!anyDuplicated(paste(cv$outcome, cv$id)), all(is.na(cv$female) | (cv$female >= 0 & cv$female <= 1)))
mapping <- read.csv("data/extraction/arm_mapping_and_demographics.csv")
cn <- c(CAIT = "CAIT", `FAAM-ADL` = "FAAM_ADL", `运动功能` = "Sports", `前向平衡` = "Anterior", `后内侧平衡` = "Posteromedial", `后外侧平衡` = "Posterolateral")
mapping$outcome <- unname(cn[mapping$outcome])
mapping <- mapping[mapping$r050 %in% c(TRUE, "True", "TRUE"), ]
mapping$id <- mapping$analysis_id; mapping$t <- mapping$node
j <- match(makekey(a), makekey(mapping))
stopifnot(!anyNA(j), nrow(a) == nrow(mapping))
stopifnot(max(abs(a$mean - mapping$input_change_mean[j])) < 1e-12,
          max(abs(a$sd - mapping$input_change_sd[j])) < 1e-12,
          identical(as.numeric(a$n), as.numeric(mapping$n_analysis[j])))
rows <- lapply(keys, function(k) {
  d <- a[a$outcome == k, ]
  r <- read.csv(file.path("data/analysis_inputs/primary", paste0(k, "_arms.csv")))
  m <- match(paste(d$id, d$t), paste(r$id, r$t))
  stopifnot(!anyNA(m), max(abs(as.matrix(d[, c("n", "mean", "sd")]) - as.matrix(r[m, c("n", "mean", "sd")]))) < 1e-12)
  data.frame(outcome = k, analysis_ids = length(unique(d$id)), arm_rows = nrow(d), analysis_n_sum = sum(d$n), input_agreement = TRUE)
})
write.csv(do.call(rbind, rows), "results/validation/input_checks.csv", row.names = FALSE)
cat("Inputs agree across primary, sensitivity, subgroup, regression and extracted arm mapping.\n")
