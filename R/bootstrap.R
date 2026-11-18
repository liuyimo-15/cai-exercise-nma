# All paths resolve from the repository root. No workstation-specific paths.
options(stringsAsFactors = FALSE, repos = c(CRAN = "https://cloud.r-project.org"))
if (!grepl("UTF-8|utf8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
  for (loc in c("C.UTF-8", "en_US.UTF-8", "English_United States.utf8")) {
    ok <- suppressWarnings(Sys.setlocale("LC_CTYPE", loc))
    if (nzchar(ok)) break
  }
}
libs <- c(Sys.getenv("NMA_R_LIBRARY", unset = ""), ".Rlibrary")
libs <- libs[nzchar(libs) & dir.exists(libs)]
.libPaths(c(normalizePath(libs, winslash = "/"), .libPaths()))
set.seed(20261118)

# Choose a usable native PNG device; capabilities('cairo') alone can be TRUE
# even when a Mac lacks the external XQuartz shared libraries.
nma_png_type <- function() {
  if (.Platform$OS.type == "windows") return("windows")
  if (Sys.info()[["sysname"]] == "Darwin") return("quartz")
  "cairo"
}
