# ============================================================
# 00_theory_config.R
# Stand-alone theoretical model of taxonomic uncertainty
# ============================================================
#
# This module is intentionally independent from the empirical
# Collembola workflow. It does not source 00_config.R from the
# empirical project and does not write into its output directory.
#
# The model implements:
#   true community N
#   -> identification process K
#   -> reporting / workflow process C
#   -> observed community Z = N K C
#
# It is designed for cross-taxon simulations (e.g. Collembola,
# earthworms, ants), not for calibration on a single taxon.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(vegan)
})

.theory_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)

THEORY_OUT_DIR <- file.path(THEORY_DIR, "outputs")
dir.create(THEORY_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

THEORY_SEED <- getOption("taxo.theory_seed", 123)

assert_scalar <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  ok <- length(x) == 1 && is.finite(x) && x >= lower && x <= upper
  if (integer) ok <- ok && abs(x - round(x)) < .Machine$double.eps^0.5
  if (!ok) {
    stop(
      sprintf("`%s` must be one finite scalar in [%s, %s].", name, lower, upper),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

assert_named_list <- function(x, name) {
  if (!is.list(x) || is.null(names(x))) {
    stop(sprintf("`%s` must be a named list.", name), call. = FALSE)
  }
  invisible(TRUE)
}

clamp <- function(x, lower = 0, upper = 1) pmin(pmax(x, lower), upper)

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)
}

safe_quantile <- function(x, p) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else {
    as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
  }
}

relative_change <- function(value, baseline) {
  ifelse(is.finite(baseline) & baseline != 0, 100 * (value - baseline) / baseline, NA_real_)
}

write_theory_csv <- function(x, filename) {
  readr::write_csv(x, file.path(THEORY_OUT_DIR, filename))
}
