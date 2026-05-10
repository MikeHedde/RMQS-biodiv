# ============================================================
# 14_reviewer_analyses.R
# Additional reviewer-requested analyses and tables
#
# To be sourced after:
#   03_read_parse.R, 04_detection_effort_reviewed.R, 05_occupancy_fit.R,
#   10_diversity_indices.R
#
# Main outputs:
#   reviewer_response/coverage_095_summary.csv
#   reviewer_response/coverage_095_by_habitat.csv
#   reviewer_response/occupancy_species_representativeness.csv
#   reviewer_response/gpd_il_sensitivity.csv
#   reviewer_response/pitfall_subset_spatial_sensitivity_*.csv
# ============================================================

stopifnot(exists("dat0"), exists("out_dir"))
reviewer_dir <- file.path(out_dir, "reviewer_response")
dir.create(reviewer_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Helpers ----------------------------------------------------------
row_coverage_rr <- function(v) {
  n <- sum(v, na.rm = TRUE)
  if (n <= 0) return(0)
  f1 <- sum(v == 1, na.rm = TRUE)
  1 - f1 / n
}

read_if_exists <- function(path) {
  if (file.exists(path)) readr::read_csv(path, show_col_types = FALSE) else NULL
}

# ---------- A. Coverage >= 0.95 without extrapolation ------------------------
# Uses observed coverage from 10_diversity_indices.R. If div_cov2 exists in memory,
# use it; otherwise read it from the output folder.
if (exists("div_cov2")) {
  div_cov_reviewer <- div_cov2
} else {
  f1 <- file.path(out_dir, "diversity", "Hill_div_results.csv")
  f2 <- "data/derived-data/Hill_div_results.csv"
  if (file.exists(f1)) div_cov_reviewer <- read.csv(f1)
  else if (file.exists(f2)) div_cov_reviewer <- read.csv(f2)
  else div_cov_reviewer <- NULL
}

if (!is.null(div_cov_reviewer)) {
  div_cov_reviewer <- div_cov_reviewer %>%
    dplyr::mutate(
      method = as.character(method),
      reached_095 = coverage >= 0.95
    )

  coverage_095_summary <- div_cov_reviewer %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      n_cases = dplyr::n(),
      n_reached_095 = sum(reached_095, na.rm = TRUE),
      pct_reached_095 = 100 * mean(reached_095, na.rm = TRUE),
      coverage_median = stats::median(coverage, na.rm = TRUE),
      coverage_q1 = stats::quantile(coverage, 0.25, na.rm = TRUE),
      coverage_q3 = stats::quantile(coverage, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  readr::write_csv(coverage_095_summary, file.path(reviewer_dir, "coverage_095_summary.csv"))

  habitat_col <- intersect(c("HABITAT_CODE", "HABITAT_OC", "habitat", "HABITAT"), names(div_cov_reviewer))[1]
  if (!is.na(habitat_col)) {
    coverage_095_by_habitat <- div_cov_reviewer %>%
      dplyr::group_by(method, habitat = .data[[habitat_col]]) %>%
      dplyr::summarise(
        n_cases = dplyr::n(),
        n_reached_095 = sum(reached_095, na.rm = TRUE),
        pct_reached_095 = 100 * mean(reached_095, na.rm = TRUE),
        coverage_median = stats::median(coverage, na.rm = TRUE),
        .groups = "drop"
      )
    readr::write_csv(coverage_095_by_habitat, file.path(reviewer_dir, "coverage_095_by_habitat.csv"))
  }
}

# ---------- B. Representativeness of occupancy-modeled species ---------------
# Compares modeled species with the full species pool using frequency, abundance,
# and optionally traits when trait file is available.
occ_p_path <- file.path(out_dir, "occ_model_output", "p_hat_by_method_unmarked_full.csv")
tab_p_rr <- read_if_exists(occ_p_path)
modeled_species <- if (!is.null(tab_p_rr) && "species" %in% names(tab_p_rr)) unique(tab_p_rr$species) else character(0)

species_base <- dat0 %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_sites_detected = dplyr::n_distinct(site_id[det == 1]),
    total_abundance = sum(as.numeric(abund), na.rm = TRUE),
    n_records = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    occupancy_modelled = species %in% modeled_species
  )

trait_path <- "data/raw-data/1.faune/trait_carabidae.csv"
if (file.exists(trait_path)) {
  trait0_rr <- read.csv(trait_path, header = TRUE, sep = ";")
  # Try flexible joins; adjust manually if your trait file uses different columns.
  sp_trait_col <- intersect(c("species", "LB_NOM", "taxon", "Species", "Taxon"), names(trait0_rr))[1]
  if (!is.na(sp_trait_col)) {
    trait0_rr <- trait0_rr %>% dplyr::rename(species = dplyr::all_of(sp_trait_col))
    species_base <- species_base %>% dplyr::left_join(trait0_rr, by = "species")
  }
}

readr::write_csv(species_base, file.path(reviewer_dir, "occupancy_species_representativeness_by_species.csv"))

repr_summary <- species_base %>%
  dplyr::group_by(occupancy_modelled) %>%
  dplyr::summarise(
    n_species = dplyr::n(),
    median_sites_detected = stats::median(n_sites_detected, na.rm = TRUE),
    q1_sites_detected = stats::quantile(n_sites_detected, 0.25, na.rm = TRUE),
    q3_sites_detected = stats::quantile(n_sites_detected, 0.75, na.rm = TRUE),
    median_total_abundance = stats::median(total_abundance, na.rm = TRUE),
    q1_total_abundance = stats::quantile(total_abundance, 0.25, na.rm = TRUE),
    q3_total_abundance = stats::quantile(total_abundance, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(repr_summary, file.path(reviewer_dir, "occupancy_species_representativeness_summary.csv"))

# ---------- C. GPD interception-length sensitivity ---------------------------
# This does not refit occupancy models; it documents how alternative IL values
# affect the effort-equivalent positioning of GPD relative to Pitfall k.
il_values <- if (exists("gpd_il_values")) gpd_il_values else c(0.80, 1.00, 1.16, 1.30, 1.50)

if (exists("effort_wide") && "days_gpd" %in% names(effort_wide)) {
  pit_eff_cols <- paste0("eff_p", pit_reps)
  pit_eff_cols <- pit_eff_cols[pit_eff_cols %in% names(effort_wide)]

  gpd_il_sensitivity <- purrr::map_dfr(il_values, function(il) {
    tmp <- effort_wide %>%
      dplyr::mutate(eff_gpd_alt = il * days_gpd * dplyr::coalesce(n_gpd_units, 1L))

    med_gpd <- stats::median(tmp$eff_gpd_alt, na.rm = TRUE)
    med_pits <- sapply(pit_eff_cols, function(cc) stats::median(tmp[[cc]], na.rm = TRUE))
    closest <- names(which.min(abs(med_pits - med_gpd)))

    tibble::tibble(
      IL_GPD_m = il,
      median_eff_gpd = med_gpd,
      closest_pitfall_effort = sub("eff_p", "Pitfall", closest),
      median_eff_closest_pitfall = unname(med_pits[closest]),
      ratio_gpd_to_closest = med_gpd / unname(med_pits[closest])
    )
  })

  readr::write_csv(gpd_il_sensitivity, file.path(reviewer_dir, "gpd_il_sensitivity.csv"))
}

# ---------- D. Pitfall subset spatial sensitivity ----------------------------
# Requires trap coordinates. Define in 02_config.R for your machine, e.g.:

if (exists("trap_coord_path") && file.exists(trap_coord_path)) {
  trap_xy0 <- read.csv(trap_coord_path, header = TRUE, sep = ifelse(grepl("\\.csv$", trap_coord_path), ";", ","))

  # Flexible column standardisation
  nm <- names(trap_xy0)
  site_col <- intersect(c("site_id", "NOM_STATION", "station", "site"), nm)[1]
  trap_col <- intersect(c("trap", "TRAP", "trap_id", "PB", "piege", "piège"), nm)[1]
  x_col <- intersect(c("x", "X", "coord_x", "x_m"), nm)[1]
  y_col <- intersect(c("y", "Y", "coord_y", "y_m"), nm)[1]

  if (is.na(trap_col) || is.na(x_col) || is.na(y_col)) {
    warning("Trap coordinate file found, but trap/x/y columns were not recognised. Spatial sensitivity skipped.")
  } else {
    trap_xy <- trap_xy0 %>%
      dplyr::transmute(
        site_id = if (!is.na(site_col)) as.character(.data[[site_col]]) else NA_character_,
        trap = as.integer(.data[[trap_col]]),
        x = as.numeric(.data[[x_col]]),
        y = as.numeric(.data[[y_col]])
      ) %>%
      dplyr::filter(!is.na(trap), !is.na(x), !is.na(y))

    # If coordinates are generic, expand to all sites.
    if (all(is.na(trap_xy$site_id))) {
      trap_xy <- tidyr::crossing(site_id = sort(unique(dat0$site_id)), trap_xy %>% dplyr::select(-site_id))
    }

    pairwise_stats <- function(df) {
      if (nrow(df) <= 1) {
        return(tibble::tibble(mean_pairwise_distance = 0, max_pairwise_distance = 0, convex_hull_area = 0,
                              centroid_distance = sqrt(mean(df$x)^2 + mean(df$y)^2)))
      }
      d <- as.matrix(stats::dist(df[, c("x", "y")]))
      upper <- d[upper.tri(d)]
      hull_area <- 0
      if (nrow(df) >= 3) {
        h <- grDevices::chull(df$x, df$y)
        xx <- df$x[h]; yy <- df$y[h]
        hull_area <- 0.5 * abs(sum(xx * c(yy[-1], yy[1]) - yy * c(xx[-1], xx[1])))
      }
      tibble::tibble(
        mean_pairwise_distance = mean(upper, na.rm = TRUE),
        max_pairwise_distance = max(upper, na.rm = TRUE),
        convex_hull_area = hull_area,
        centroid_distance = sqrt(mean(df$x, na.rm = TRUE)^2 + mean(df$y, na.rm = TRUE)^2)
      )
    }

    # Species x site abundance by individual trap
    pit_trap_sp <- dat0 %>%
      dplyr::filter(method == "Pitfall", !is.na(trap), trap %in% 1:10) %>%
      dplyr::group_by(site_id, trap, species) %>%
      dplyr::summarise(abund = sum(as.numeric(abund), na.rm = TRUE), .groups = "drop")

    ks <- c(2, 4, 6, 8, 10)
    
    draw_one_site <- function(site) {
      traps_site <- trap_xy %>% dplyr::filter(site_id == site, trap %in% 1:10) %>% dplyr::distinct(trap, .keep_all = TRUE)
      if (nrow(traps_site) == 0) return(tibble::tibble())
      available <- sort(traps_site$trap)

      purrr::map_dfr(ks, function(k) {
        if (length(available) < k) return(tibble::tibble())
        
        # Enumerate all possible combinations of k traps among available traps.
        # For 10 traps: C(10,2)=45, C(10,4)=210, C(10,6)=210,
        # C(10,8)=45, C(10,10)=1.
        draws <- combn(available, k, simplify = FALSE)

        
      purrr::imap_dfr(draws, function(tr_set, draw_id) {
          xy_sub <- traps_site %>% dplyr::filter(trap %in% tr_set)
          stats_sub <- pairwise_stats(xy_sub)

          sp_sub <- pit_trap_sp %>%
            dplyr::filter(site_id == site, trap %in% tr_set) %>%
            dplyr::group_by(species) %>%
            dplyr::summarise(abund = sum(abund, na.rm = TRUE), detected = as.integer(abund > 0), .groups = "drop")

          tibble::tibble(
            site_id = site,
            k = k,
            draw_id = draw_id,
            trap_set = paste(tr_set, collapse = ","),
            n_species = sum(sp_sub$detected > 0, na.rm = TRUE),
            total_abundance = sum(sp_sub$abund, na.rm = TRUE),
            coverage = row_coverage_rr(sp_sub$abund)
          ) %>%
            dplyr::bind_cols(stats_sub)
        })
      })
    }

    pitfall_subset_sensitivity <- purrr::map_dfr(sort(unique(dat0$site_id)), draw_one_site)
    readr::write_csv(pitfall_subset_sensitivity, file.path(reviewer_dir, "pitfall_subset_spatial_sensitivity_draws.csv"))

    pitfall_subset_sensitivity_summary <- pitfall_subset_sensitivity %>%
      dplyr::group_by(k) %>%
      dplyr::summarise(
        n_site_draws = dplyr::n(),
        richness_median = stats::median(n_species, na.rm = TRUE),
        richness_q1 = stats::quantile(n_species, 0.25, na.rm = TRUE),
        richness_q3 = stats::quantile(n_species, 0.75, na.rm = TRUE),
        coverage_median = stats::median(coverage, na.rm = TRUE),
        mean_pairwise_distance_median = stats::median(mean_pairwise_distance, na.rm = TRUE),
        convex_hull_area_median = stats::median(convex_hull_area, na.rm = TRUE),
        .groups = "drop"
      )
    readr::write_csv(pitfall_subset_sensitivity_summary, file.path(reviewer_dir, "pitfall_subset_spatial_sensitivity_summary.csv"))

    # How much do spatial metrics explain richness/coverage among subsets?
    pitfall_subset_models <- list()
    if (nrow(pitfall_subset_sensitivity) > 0) {
      pitfall_subset_models$richness <- stats::lm(n_species ~ k + mean_pairwise_distance + convex_hull_area,
                                                  data = pitfall_subset_sensitivity)
      pitfall_subset_models$coverage <- stats::lm(coverage ~ k + mean_pairwise_distance + convex_hull_area,
                                                  data = pitfall_subset_sensitivity)
      capture.output(summary(pitfall_subset_models$richness),
                     file = file.path(reviewer_dir, "pitfall_subset_model_richness.txt"))
      capture.output(summary(pitfall_subset_models$coverage),
                     file = file.path(reviewer_dir, "pitfall_subset_model_coverage.txt"))
    }
  }
} else {
  message("No trap_coord_path found. Define trap_coord_path in 02_config.R to run pitfall spatial sensitivity.")
}

# --- SPATIAL EXTENT ANALYSIS (from saved draws) ---

# 1. Read subsets 
subset_draws <- read.csv(
  file.path(out_dir, "reviewer_response/pitfall_subset_spatial_sensitivity_draws.csv")
)

spatial_vs_detection <- subset_draws %>%
  group_by(k) %>%
  summarise(
    mean_coverage = mean(coverage, na.rm = TRUE),
    sd_coverage = sd(coverage, na.rm = TRUE),
    mean_richness = mean(n_species, na.rm = TRUE),
    sd_richness = sd(n_species, na.rm = TRUE),
    mean_max_dist = mean(max_pairwise_distance, na.rm = TRUE),
    sd_max_dist = sd(max_pairwise_distance, na.rm = TRUE),
    mean_pairwise_dist = mean(mean_pairwise_distance, na.rm = TRUE),
    mean_hull_area = mean(convex_hull_area, na.rm = TRUE),
    .groups = "drop"
  )

spatial_cor <- subset_draws %>%
  summarise(
    cor_coverage_maxdist = cor(coverage, max_pairwise_distance, use = "complete.obs"),
    cor_richness_maxdist = cor(n_species, max_pairwise_distance, use = "complete.obs"),
    cor_abundance_maxdist = cor(total_abundance, max_pairwise_distance, use = "complete.obs"),
    cor_coverage_hull = cor(coverage, convex_hull_area, use = "complete.obs"),
    cor_richness_hull = cor(n_species, convex_hull_area, use = "complete.obs")
  )

write.csv(
  spatial_vs_detection,
  file.path(out_dir, "reviewer_response/pitfall_spatial_extent_by_k.csv"),
  row.names = FALSE
)

write.csv(
  spatial_cor,
  file.path(out_dir, "reviewer_response/pitfall_spatial_extent_correlations.csv"),
  row.names = FALSE
)

spatial_vs_detection
spatial_cor


# Figures
library(ggplot2)
library(dplyr)

df_site <- pitfall_subset_sensitivity %>%
  dplyr::group_by(site_id, k) %>%
  dplyr::summarise(
    richness_mean = mean(n_species, na.rm = TRUE),
    richness_sd   = sd(n_species, na.rm = TRUE),
    coverage_mean = mean(coverage, na.rm = TRUE),
    .groups = "drop"
  )

df_k <- df_site %>%
  dplyr::group_by(k) %>%
  dplyr::summarise(
    mean_richness = mean(richness_mean, na.rm = TRUE),
    q25 = quantile(richness_mean, 0.25, na.rm = TRUE),
    q75 = quantile(richness_mean, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

nested_labels <- tibble::tibble(
  k = c(2, 4, 6, 8, 10),
  trap_set = c(
    "1,2",
    "1,2,3,4",
    "1,2,3,4,5,6",
    "1,2,3,4,5,6,7,8",
    "1,2,3,4,5,6,7,8,9,10"
  )
)

nested_site <- pitfall_subset_sensitivity %>%
  dplyr::inner_join(nested_labels, by = c("k", "trap_set")) %>%
  dplyr::transmute(
    site_id,
    k,
    richness_nested = n_species,
    coverage_nested = coverage
  )


ggplot() +
  geom_line(
    data = df_site,
    aes(x = k, y = richness_mean, group = site_id),
    color = "grey40",
    alpha = 0.30,
    linewidth = 0.7) +
  geom_point(
    data = df_site,
    aes(x = k, y = richness_mean, group = site_id),
    color = "grey40",
    alpha = 0.45,
    size = 1.2) +
  geom_line(
    data = nested_site,
    aes(x = k, y = richness_nested, group = site_id),
    color = "blue",
    alpha = 0.2)+
  geom_ribbon(
    data = df_k,
    aes(x = k, ymin = q25, ymax = q75),
    fill = "grey20",
    alpha = 0.4) +
  geom_line(
    data = df_k,
    aes(x = k, y = mean_richness),
    color = "red",
    linewidth = 1.4) +
  geom_point(
    data = df_k,
    aes(x = k, y = mean_richness),
    color = "red",
    size = 3) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  labs(
    x = "Number of pitfall traps",
    y = "Mean species richness per site",
    title = "Robustness of the pitfall effort gradient across sites"
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)

# df = pitfall_subset_spatial_sensitivity_draws.csv
# colonnes attendues : site_id, k, trap_set, n_species, coverage


# --- Résumés par station ---
df_site_k <- pitfall_subset_sensitivity %>%
  mutate(
    k = as.numeric(k),
    site_lab = str_extract(site_id, "[^_]+$")
  ) %>%
  group_by(site_id, site_lab, k) %>%
  summarise(
    mean_richness = mean(n_species, na.rm = TRUE),
    q25 = quantile(n_species, 0.25, na.rm = TRUE),
    q75 = quantile(n_species, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

nested_site <- pitfall_subset_sensitivity %>%
  mutate(
    k = as.numeric(k),
    site_lab = str_extract(site_id, "[^_]+$")
  ) %>%
  inner_join(nested_labels, by = c("k", "trap_set")) %>%
  transmute(
    site_id,
    site_lab,
    k,
    richness_nested = n_species
  )

# --- Figure par station ---
p_sites <- ggplot() +
  geom_ribbon(
    data = df_site_k,
    aes(x = k, ymin = q25, ymax = q75),
    fill = "grey70",
    alpha = 0.35
  ) +
  geom_line(
    data = df_site_k,
    aes(x = k, y = mean_richness),
    color = "black",
    linewidth = 0.45
  ) +
  geom_point(
    data = df_site_k,
    aes(x = k, y = mean_richness),
    color = "black",
    size = 0.8
  ) +
  geom_line(
    data = nested_site,
    aes(x = k, y = richness_nested),
    color = "blue",
    linewidth = 0.55,
    alpha = 0.8
  ) +
  geom_point(
    data = nested_site,
    aes(x = k, y = richness_nested),
    color = "blue",
    size = 0.9,
    alpha = 0.9
  ) +
  facet_wrap(~ site_lab) +
  scale_x_continuous(breaks = c(2, 6, 10)) +
  labs(
    x = "Number of pitfall traps",
    y = "Species richness",
    title = "Site-level robustness of the pitfall effort gradient"
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_blank(),   
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.position = "none"
  )

# --- Figure générale ---
df_global <- df_site_k %>%
  group_by(k) %>%
  summarise(
    grand_mean = mean(mean_richness, na.rm = TRUE),
    q25 = quantile(mean_richness, 0.25, na.rm = TRUE),
    q75 = quantile(mean_richness, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

nested_global <- nested_site %>%
  group_by(k) %>%
  summarise(
    nested_mean = mean(richness_nested, na.rm = TRUE),
    nested_q25 = quantile(richness_nested, 0.25, na.rm = TRUE),
    nested_q75 = quantile(richness_nested, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p_global <- ggplot() +
  geom_ribbon(
    data = df_global,
    aes(x = k, ymin = q25, ymax = q75),
    fill = "grey60",
    alpha = 0.35
  ) +
  geom_line(
    data = df_global,
    aes(x = k, y = grand_mean),
    color = "black",
    linewidth = 1.2
  ) +
  geom_point(
    data = df_global,
    aes(x = k, y = grand_mean),
    color = "black",
    size = 2.5
  ) +
  geom_line(
    data = nested_global,
    aes(x = k, y = nested_mean),
    color = "blue",
    linewidth = 1.2
  ) +
  geom_point(
    data = nested_global,
    aes(x = k, y = nested_mean),
    color = "blue",
    size = 2.5
  ) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  labs(
    x = "Number of pitfall traps",
    y = "Mean species richness",
    title = "Overall pattern across sites"
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# --- Figure combinée ---
p_final <- p_global / p_sites +
  plot_layout(heights = c(1, 3))

p_final
ggsave(
  file.path(out_dir, "reviewer_response/FigS_random_species_richness.png"),
  p_final,
  width = 9,
  height = 10,
  dpi = 300
)


p1 <- ggplot(pitfall_subset_sensitivity, aes(x = factor(k), y = n_species)) +
  geom_violin(fill = "grey85", color = NA) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
  labs(y = "Richness", x = "k") +
  theme_minimal()

p2 <- ggplot(pitfall_subset_sensitivity, aes(x = factor(k), y = coverage)) +
  geom_violin(fill = "grey85", color = NA) +
  stat_summary(fun = mean, geom = "point", color = "blue", size = 3) +
  labs(y = "Coverage", x = "k") +
  theme_minimal()

library(patchwork)
p3 <- p1 + p2

ggsave(
  file.path(out_dir, "reviewer_response/FigS_random_coverage.png"),
  p3,
  width = 9,
  height = 10,
  dpi = 300
)

#Coverage 95%
library(scales)
n_habitat <- coverage_095_by_habitat %>%
  group_by(habitat) %>%
  summarise(
    n_sites = max(n_cases)
  )

coverage_by_habitat <- coverage_095_by_habitat %>%
  filter(!is.na(method), method != "NA") %>%
  mutate(
    method = factor(
      method,
      levels = c("Pitfall2", "Pitfall4", "Pitfall6", "Pitfall8", "Pitfall10", "GPD")
    ),
    habitat = factor(
      habitat,
      levels = c("Crop", "Grassland / Meadow", "Forest", "Schrubland")
    ),
    prop_095 = pct_reached_095 / 100
  )

coverage_by_habitat1 <- coverage_by_habitat %>%
  left_join(n_habitat) %>%
  mutate(
    habitat_lab = paste0(habitat, " (n = ", n_sites, ")")
  ) %>%
  filter(!is.na(habitat))

cov_by_habitat <- ggplot(coverage_by_habitat1, aes(x = method, y = pct_reached_095, fill = method)) +
  geom_col(width = 0.75, alpha = 0.9) +
  facet_wrap(~ habitat_lab, ncol = 2) +
  scale_y_continuous(,
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_fill_manual(
    values = c(
      "Pitfall2"  = "#F4A261",
      "Pitfall4"  = "#E76F51",
      "Pitfall6"  = "#2A9D8F",
      "Pitfall8"  = "#457B9D",
      "Pitfall10" = "#1D3557",
      "GPD"       = "#9B5DE5"
    )
  ) +
  labs(
    x = "Protocol",
    y = "Proportion reaching coverage ≥ 0.95",
    title = "Coverage attainment by habitat and protocol"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "reviewer_response/FigS_random_coverage_by_habitat.png"),
  cov_by_habitat,
  width = 9,
  height = 10,
  dpi = 300
)

# ============================================================
# APPENDIX — REPRESENTATIVENESS OF MODELED SPECIES
# ============================================================

dir.create(file.path(out_dir, "trait"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Full observed community
# -----------------------------
# Adapter seulement si ton fichier source a un autre nom
comm_file <- "data/raw-data/1.faune/carabidae.csv"

carab_raw <- read.csv(
  comm_file,
  h = TRUE,
  sep = ";",
  fileEncoding = "latin1",
  check.names = FALSE
)

# Détection automatique des colonnes
sp_col <- intersect(
  c("species", "Species", "NOM_VALIDE", "Nom_valide", "nom_valide",
    "taxon", "Taxon", "nom_taxon", "NOM_TAXON"),
  names(carab_raw)
)[1]

ab_col <- intersect(
  c("ABONDANCE", "Abondance", "abondance",
    "ABONDANCE_TOTALE", "abondance_totale", "count", "n"),
  names(carab_raw)
)[1]

site_col <- intersect(
  c("STATION", "NOM_STATION", "Station", "station", "SITE", "Site", "site",
    "ID_SITE", "id_site"),
  names(carab_raw)
)[1]

stopifnot(!is.na(sp_col))
stopifnot(!is.na(site_col))

if (is.na(ab_col)) {
  carab_raw$.__abundance__ <- 1
  ab_col <- ".__abundance__"
}

# -----------------------------
# 2) Canonise species names
# -----------------------------
canonise_species <- function(x) {
  tmp <- tibble(species_raw = unique(x)) %>%
    mutate(parsed = rgnparser::gn_parse(species_raw)) %>%
    tidyr::unnest_wider(parsed)
  
  tmp %>%
    transmute(
      species_raw,
      species = purrr::map_chr(canonical, \(z) {
        if (is.null(z) || length(z) == 0) return(NA_character_)
        if (!is.null(z$simple)  && length(z$simple)  > 0 && nzchar(z$simple[1]))  return(z$simple[1])
        if (!is.null(z$full)    && length(z$full)    > 0 && nzchar(z$full[1]))    return(z$full[1])
        if (!is.null(z$stemmed) && length(z$stemmed) > 0 && nzchar(z$stemmed[1])) return(z$stemmed[1])
        NA_character_
      })
    )
}

sp_map_comm <- canonise_species(carab_raw[[sp_col]])

comm_sp <- carab_raw %>%
  rename(
    species_raw = all_of(sp_col),
    site        = all_of(site_col),
    abundance   = all_of(ab_col)
  ) %>%
  left_join(sp_map_comm, by = "species_raw") %>%
  filter(!is.na(species), species != "") %>%
  mutate(abundance = as.numeric(abundance)) %>%
  group_by(species) %>%
  summarise(
    n_sites = n_distinct(site[abundance > 0]),
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

modeled_species <- unique(tab_p$species)

repr_dat <- comm_sp %>%
  left_join(BL, by = "species") %>%
  left_join(trait_sp, by = "species") %>%
  mutate(
    modeled = ifelse(species %in% modeled_species,
                     "Modelled species",
                     "Non-modelled species"),
    log_total_abundance = log10(total_abundance + 1),
    log_n_sites = log10(n_sites + 1),
    log_BL = log10(mean_BL),
    full_wing = case_when(
      pct_winged == 1 ~ "fully winged",
      pct_winged <  1 ~ "not fully winged",
      TRUE ~ NA_character_
    )
  )

readr::write_csv(
  repr_dat,
  file.path(out_dir, "trait/species_representativeness_dataset.csv")
)

# -----------------------------
# 3) Summary table
# -----------------------------
repr_summary <- repr_dat %>%
  group_by(modeled) %>%
  summarise(
    n_species = n(),
    total_abundance = sum(total_abundance, na.rm = TRUE),
    median_n_sites = median(n_sites, na.rm = TRUE),
    median_body_length = median(mean_BL, na.rm = TRUE),
    min_body_length = min(mean_BL, na.rm = TRUE),
    max_body_length = max(mean_BL, na.rm = TRUE),
    prop_fully_winged = mean(full_wing == "fully winged", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prop_total_abundance =
      total_abundance / sum(total_abundance)
  )

readr::write_csv(
  repr_summary,
  file.path(out_dir, "trait/species_representativeness_summary.csv")
)

# -----------------------------
# 4) Statistical comparisons
# -----------------------------
sink(file.path(out_dir, "trait/species_representativeness_tests.txt"))

cat("Representativeness of modelled species\n")
cat("=====================================\n\n")

cat("Number of species:\n")
print(table(repr_dat$modeled))

cat("\nContribution to total abundance:\n")
print(repr_summary[, c("modeled", "n_species", "prop_total_abundance")])

cat("\nBody length comparison:\n")
print(wilcox.test(mean_BL ~ modeled, data = repr_dat))

cat("\nNumber of occupied sites comparison:\n")
print(wilcox.test(n_sites ~ modeled, data = repr_dat))

cat("\nTotal abundance comparison:\n")
print(wilcox.test(total_abundance ~ modeled, data = repr_dat))

cat("\nWing category comparison:\n")
print(table(repr_dat$modeled, repr_dat$full_wing))
print(chisq.test(table(repr_dat$modeled, repr_dat$full_wing)))

sink()

# -----------------------------
# 5) Figure S — Representativeness
# -----------------------------
repr_long <- repr_dat %>%
  select(species, modeled, log_BL, log_n_sites, log_total_abundance) %>%
  pivot_longer(
    cols = c(log_BL, log_n_sites, log_total_abundance),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = dplyr::recode(
      variable,
      log_BL = "Trait: Body length\nlog10(mm)",
      log_n_sites = "Occurrence frequency\nlog10(number of sites + 1)",
      log_total_abundance = "Total abundance\nlog10(N + 1)"
    )
  ) %>%
  filter(!is.na(value), is.finite(value))

p_repr <- ggplot(repr_long, aes(x = value, fill = modeled)) +
  geom_density(alpha = 0.35) +
  facet_wrap(~ variable, scales = "free", ncol = 1) +
  labs(
    x = NULL,
    y = "Density",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.position = "top"
  )

ggsave(
  file.path(out_dir, "trait/FigS_representativeness_modeled_species.png"),
  p_repr,
  width = 9,
  height = 10,
  dpi = 300
)

# -----------------------------
# 6) Figure S — Modeled species on occurrence-abundance space
# -----------------------------
p_occ_abund <- ggplot(
  repr_dat,
  aes(x = log_n_sites, y = log_total_abundance)
) +
  geom_point(aes(colour = modeled), alpha = 0.35, size = 3) +
  labs(
    x = "Occurrence frequency log10(number of sites + 1)",
    y = "Total abundance log10(N + 1)",
    shape = NULL
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.position = "top"
  )

ggsave(
  file.path(out_dir, "trait/FigS_occurrence_abundance_modeled_species.png"),
  p_occ_abund,
  width = 9,
  height = 7,
  dpi = 300
)

# -----------------------------
# 7) Figure S — Wing development representativeness
# -----------------------------

wing_dat <- repr_dat %>%
  filter(!is.na(full_wing)) %>%
  mutate(
    modeled = factor(modeled, levels = c("Non-modelled species", "Modelled species")),
    full_wing = factor(full_wing, levels = c("not fully winged", "fully winged"))
  ) %>%
  count(modeled, full_wing) %>%
  group_by(modeled) %>%
  mutate(
    prop = n / sum(n),
    label = paste0(n, " sp.\n", scales::percent(prop, accuracy = 1))
  ) %>%
  ungroup()

p_wing_repr <- ggplot(
  wing_dat,
  aes(x = modeled, y = prop, fill = full_wing)
) +
  geom_col(width = 0.7, alpha = 0.35) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 5
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Proportion of species",
    fill = "Wing development"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.position = "top"
  )

ggsave(
  file.path(out_dir, "trait/FigS_wing_development_representativeness.png"),
  p_wing_repr,
  width = 8,
  height = 7,
  dpi = 300
)

message("=== FIN 14_reviewer_analyses ===")
