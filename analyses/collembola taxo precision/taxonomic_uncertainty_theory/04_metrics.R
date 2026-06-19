# ============================================================
# 04_metrics.R
# Alpha, gamma, occupancy, beta-diversity and Blowes contrasts
# ============================================================

hill_numbers_one_site <- function(x) {
  x <- as.numeric(x)
  x <- x[x > 0]

  if (length(x) == 0) {
    return(c(q0 = NA_real_, q1 = NA_real_, q2 = NA_real_))
  }

  p <- x / sum(x)
  c(
    q0 = length(x),
    q1 = exp(-sum(p * log(p))),
    q2 = 1 / sum(p^2)
  )
}

mean_pairwise_distance <- function(comm, method, binary = FALSE) {
  if (nrow(comm) < 2 || ncol(comm) < 1 || all(rowSums(comm) == 0)) return(NA_real_)

  d <- tryCatch(
    vegan::vegdist(comm, method = method, binary = binary),
    error = function(e) NULL
  )

  if (is.null(d) || length(d) == 0) NA_real_ else mean(d, na.rm = TRUE)
}

community_metrics <- function(comm) {
  if (!is.matrix(comm)) comm <- as.matrix(comm)
  storage.mode(comm) <- "numeric"

  hills <- t(apply(comm, 1, hill_numbers_one_site))
  active <- colSums(comm) > 0

  gamma <- sum(active)
  mean_alpha_q0 <- safe_mean(hills[, "q0"])
  mean_occupancy <- if (gamma > 0) {
    mean(colMeans(comm[, active, drop = FALSE] > 0))
  } else {
    NA_real_
  }

  tibble(
    n_sites = nrow(comm),
    total_abundance = sum(comm),
    alpha_q0 = mean_alpha_q0,
    alpha_q1 = safe_mean(hills[, "q1"]),
    alpha_q2 = safe_mean(hills[, "q2"]),
    gamma = gamma,
    mean_occupancy = mean_occupancy,
    alpha_over_gamma = ifelse(gamma > 0, mean_alpha_q0 / gamma, NA_real_),
    mean_bray_curtis = mean_pairwise_distance(comm, method = "bray", binary = FALSE),
    mean_jaccard = mean_pairwise_distance(comm, method = "jaccard", binary = TRUE),
    # Bray-Curtis on binary data is Sørensen dissimilarity.
    mean_sorensen = mean_pairwise_distance(comm, method = "bray", binary = TRUE)
  )
}

metric_stability <- function(comm_baseline, comm_target) {
  if (nrow(comm_baseline) != nrow(comm_target)) {
    stop("Baseline and target matrices must contain the same sites in the same order.", call. = FALSE)
  }

  alpha_base <- rowSums(comm_baseline > 0)
  alpha_target <- rowSums(comm_target > 0)

  beta_cor <- function(mat1, mat2, method, binary = FALSE) {
    d1 <- tryCatch(vegan::vegdist(mat1, method = method, binary = binary), error = function(e) NULL)
    d2 <- tryCatch(vegan::vegdist(mat2, method = method, binary = binary), error = function(e) NULL)
    if (is.null(d1) || is.null(d2) || length(d1) < 3 || length(d2) < 3) return(NA_real_)
    suppressWarnings(stats::cor(as.numeric(d1), as.numeric(d2), method = "spearman"))
  }

  tibble(
    alpha_q0_spearman = suppressWarnings(stats::cor(alpha_base, alpha_target, method = "spearman")),
    bray_spearman = beta_cor(comm_baseline, comm_target, method = "bray", binary = FALSE),
    jaccard_spearman = beta_cor(comm_baseline, comm_target, method = "jaccard", binary = TRUE),
    sorensen_spearman = beta_cor(comm_baseline, comm_target, method = "bray", binary = TRUE)
  )
}

compare_to_baseline <- function(baseline_metrics, target_metrics, tolerance_pct = 0.5) {
  metric_names <- c(
    "alpha_q0", "alpha_q1", "alpha_q2", "gamma",
    "mean_occupancy", "mean_bray_curtis", "mean_jaccard", "mean_sorensen"
  )

  out <- purrr::map_dfr(metric_names, function(metric) {
    tibble(
      metric = metric,
      baseline = baseline_metrics[[metric]],
      target = target_metrics[[metric]],
      delta = target_metrics[[metric]] - baseline_metrics[[metric]],
      relative_delta_pct = relative_change(target_metrics[[metric]], baseline_metrics[[metric]])
    )
  })

  delta_alpha <- out %>% filter(metric == "alpha_q0") %>% pull(relative_delta_pct)
  delta_gamma <- out %>% filter(metric == "gamma") %>% pull(relative_delta_pct)
  delta_occupancy <- out %>% filter(metric == "mean_occupancy") %>% pull(relative_delta_pct)

  signature <- dplyr::case_when(
    !is.finite(delta_alpha) || !is.finite(delta_gamma) ~ NA_character_,
    abs(delta_alpha) < tolerance_pct && abs(delta_gamma) < tolerance_pct ~ "no_meaningful_change",
    delta_occupancy > tolerance_pct ~ "apparent_homogenisation",
    delta_occupancy < -tolerance_pct ~ "apparent_differentiation",
    TRUE ~ "no_meaningful_change"
  )

  list(
    long = out,
    blowes = tibble(
      delta_alpha_pct = delta_alpha,
      delta_gamma_pct = delta_gamma,
      delta_occupancy_pct = delta_occupancy,
      beta_signature = signature
    )
  )
}
