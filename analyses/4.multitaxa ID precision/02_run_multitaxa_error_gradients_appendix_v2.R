# =============================================================================
# 02_run_multitaxa_error_gradients_appendix_v2.R
# Appendix: response curves from 1% to 20% taxonomic-identification error.
#
# This script is deliberately separate from the main analysis:
#   - the main paper uses a single transparent 10% reference scenario;
#   - this script quantifies how responses scale from 1% to 20%;
#   - it uses coherent source -> target maps, fixed within each iteration across
#     all error rates, so curves isolate error intensity rather than redraw
#     taxonomic targets at each rate.
#
# Prerequisites:
#   00_prepare_multitaxa_inputs_2024.R
#   05_build_multitaxa_taxref_pools.R
#
# It does NOT fit GDMs; the appendix focuses on alpha/beta response curves.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(vegan)
})

# ---- 0. Settings -------------------------------------------------------------
INPUT_DIR <- "outputs_multitaxa_2024"
OUT_DIR <- file.path(INPUT_DIR, "uncertainty_results", "appendix_error_gradients")
REGIONAL_POOL_DIR <- "regional_pools"
EXPERT_MAP_DIR <- "expert_confusion_maps"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

ASSEMBLAGES <- c(
  "collembola_soil_core",
  "araneae_pitfall",
  "carabidae_pitfall",
  "formicidae_pitfall",
  "isopoda_pitfall",
  "isopoda_hand_sorting",
  "diplopoda_pitfall",
  "diplopoda_hand_sorting"
)

ERROR_RATES <- seq(0.01, 0.20, by = 0.01)
N_SIM <- 30L
RARE_WEIGHTED_MAX_ERROR <- 0.25
MIN_BETA_SITES <- 8L

set.seed(20260619)

# ---- 1. Helpers --------------------------------------------------------------
safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || dplyr::n_distinct(x[keep]) < 2L || dplyr::n_distinct(y[keep]) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

make_matrix <- function(comm_long, station_frame) {
  stations <- as.character(station_frame$station)
  taxa <- sort(unique(comm_long$taxon_unit))
  
  if (!length(taxa)) {
    return(matrix(0, nrow = length(stations), ncol = 0,
                  dimnames = list(stations, character(0))))
  }
  
  comm_long %>%
    filter(!is.na(station), !is.na(taxon_unit), abundance > 0) %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    tidyr::complete(station = stations, taxon_unit = taxa, fill = list(abundance = 0)) %>%
    pivot_wider(names_from = taxon_unit, values_from = abundance, values_fill = 0) %>%
    arrange(match(station, stations)) %>%
    tibble::column_to_rownames("station") %>%
    as.matrix()
}

alpha_metrics <- function(mat) {
  total <- rowSums(mat)
  q0 <- rowSums(mat > 0)
  q1 <- exp(vegan::diversity(mat, index = "shannon"))
  q2 <- vegan::diversity(mat, index = "invsimpson")
  q1[!is.finite(q1)] <- 0
  q2[!is.finite(q2)] <- 0
  tibble(station = rownames(mat), q0 = q0, q1 = q1, q2 = q2, total = total)
}

prepare_beta_matrix <- function(mat) {
  keep <- rowSums(mat) > 0
  mat <- mat[keep, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (nrow(mat) < MIN_BETA_SITES || ncol(mat) < 2L) return(NULL)
  mat
}

mean_dissimilarity <- function(mat, metric) {
  mat <- prepare_beta_matrix(mat)
  if (is.null(mat)) return(NA_real_)
  if (metric == "bray") {
    d <- vegan::vegdist(mat, method = "bray")
  } else {
    d <- vegan::vegdist((mat > 0) * 1L, method = "bray", binary = TRUE)
  }
  mean(as.numeric(d), na.rm = TRUE)
}

distance_stability <- function(base, scenario, metric) {
  common <- intersect(rownames(base), rownames(scenario))
  if (length(common) < MIN_BETA_SITES) return(NA_real_)
  
  b <- base[common, , drop = FALSE]
  s <- scenario[common, , drop = FALSE]
  
  if (metric == "bray") {
    db <- vegan::vegdist(b, method = "bray")
    ds <- vegan::vegdist(s, method = "bray")
  } else {
    db <- vegan::vegdist((b > 0) * 1L, method = "bray", binary = TRUE)
    ds <- vegan::vegdist((s > 0) * 1L, method = "bray", binary = TRUE)
  }
  
  safe_cor(as.numeric(db), as.numeric(ds))
}

read_assemblage <- function(assemblage_id) {
  base <- file.path(INPUT_DIR, "assemblages", assemblage_id)
  list(
    station_frame = readr::read_csv(paste0(base, "__station_frame.csv"), show_col_types = FALSE) %>%
      mutate(station = as.character(station)),
    species = readr::read_csv(paste0(base, "__species_long.csv"), show_col_types = FALSE) %>%
      mutate(station = as.character(station)),
    lookup = readr::read_csv(paste0(base, "__taxon_lookup.csv"), show_col_types = FALSE)
  )
}

species_meta_from_lookup <- function(lookup) {
  lookup %>%
    transmute(
      taxon_unit = taxon_unit_species,
      genus = genus
    ) %>%
    filter(!is.na(taxon_unit), !is.na(genus)) %>%
    distinct()
}

read_regional_pool <- function(assemblage_id) {
  path <- file.path(REGIONAL_POOL_DIR, paste0(assemblage_id, "__regional_pool.csv"))
  if (!file.exists(path)) return(NULL)

  pool <- readr::read_csv(path, show_col_types = FALSE) %>%
    janitor::clean_names()
  if (!all(c("taxon_unit", "genus") %in% names(pool))) return(NULL)

  pool %>%
    transmute(
      taxon_unit = as.character(taxon_unit),
      genus = as.character(genus),
      cd_ref = if ("cd_ref" %in% names(pool)) as.character(cd_ref) else NA_character_
    ) %>%
    filter(!is.na(taxon_unit), !is.na(genus)) %>%
    distinct()
}

read_regional_source_audit <- function(assemblage_id) {
  path <- file.path(REGIONAL_POOL_DIR, paste0(assemblage_id, "__observed_taxref_match_audit.csv"))
  if (!file.exists(path)) return(tibble(taxon_unit = character(), cd_ref = character()))

  audit <- readr::read_csv(path, show_col_types = FALSE) %>% janitor::clean_names()
  if (!all(c("observed_taxon_unit", "cd_ref") %in% names(audit))) {
    return(tibble(taxon_unit = character(), cd_ref = character()))
  }

  audit %>%
    transmute(
      taxon_unit = as.character(observed_taxon_unit),
      cd_ref = as.character(cd_ref)
    ) %>%
    filter(!is.na(taxon_unit)) %>%
    distinct(taxon_unit, .keep_all = TRUE)
}

read_expert_map <- function(assemblage_id) {
  path <- file.path(EXPERT_MAP_DIR, paste0(assemblage_id, "__expert_confusions.csv"))
  if (!file.exists(path)) return(NULL)
  
  dat <- readr::read_csv(path, show_col_types = FALSE) %>%
    janitor::clean_names()
  
  if (!all(c("source_taxon_unit", "target_taxon_unit") %in% names(dat))) return(NULL)
  if (!"enabled" %in% names(dat)) dat$enabled <- TRUE
  if (!"weight" %in% names(dat)) dat$weight <- 1
  
  dat %>%
    mutate(
      enabled = tolower(as.character(enabled)) %in% c("true", "t", "1", "yes", "y"),
      weight = suppressWarnings(as.numeric(weight)),
      weight = coalesce(weight, 1)
    ) %>%
    filter(enabled, source_taxon_unit != target_taxon_unit) %>%
    select(source_taxon_unit, target_taxon_unit, weight) %>%
    distinct()
}

make_rare_probs <- function(species_long, rate) {
  tmp <- species_long %>%
    group_by(taxon_unit) %>%
    summarise(total_abundance = sum(abundance), .groups = "drop") %>%
    mutate(raw_weight = 1 / sqrt(total_abundance))
  
  scale_factor <- rate / weighted.mean(tmp$raw_weight, w = tmp$total_abundance)
  tmp %>%
    transmute(taxon_unit, p_error = pmin(raw_weight * scale_factor, RARE_WEIGHTED_MAX_ERROR))
}

make_coherent_map <- function(source_meta, candidate_pool, seed) {
  set.seed(seed)

  if (!"cd_ref" %in% names(source_meta)) source_meta$cd_ref <- NA_character_
  if (!"cd_ref" %in% names(candidate_pool)) candidate_pool$cd_ref <- NA_character_

  src <- source_meta %>%
    transmute(
      source_taxon_unit = taxon_unit,
      source_genus = genus,
      source_cd_ref = as.character(cd_ref)
    ) %>%
    distinct()

  cand <- candidate_pool %>%
    transmute(
      target_taxon_unit = taxon_unit,
      target_genus = genus,
      target_cd_ref = as.character(cd_ref)
    ) %>%
    distinct()

  purrr::map_dfr(seq_len(nrow(src)), function(i) {
    x <- src[i, , drop = FALSE]
    available <- cand %>%
      filter(
        target_genus == x$source_genus,
        target_taxon_unit != x$source_taxon_unit,
        is.na(x$source_cd_ref) | is.na(target_cd_ref) | target_cd_ref != x$source_cd_ref
      )

    if (!nrow(available)) {
      return(tibble(
        source_taxon_unit = x$source_taxon_unit,
        target_taxon_unit = NA_character_
      ))
    }

    tibble(
      source_taxon_unit = x$source_taxon_unit,
      target_taxon_unit = sample(available$target_taxon_unit, 1L)
    )
  })
}

make_expert_map <- function(source_meta, expert_map, seed) {
  if (is.null(expert_map) || !nrow(expert_map)) return(NULL)
  set.seed(seed)
  
  allowed_sources <- source_meta %>% pull(taxon_unit)
  expert_map %>%
    filter(source_taxon_unit %in% allowed_sources) %>%
    group_by(source_taxon_unit) %>%
    group_modify(~ {
      selected <- slice_sample(.x, n = 1L, weight_by = weight)
      tibble(target_taxon_unit = selected$target_taxon_unit)
    }) %>%
    ungroup()
}

simulate_map <- function(species_long, map, probs, seed) {
  set.seed(seed)
  
  species_long %>%
    left_join(map, by = c("taxon_unit" = "source_taxon_unit")) %>%
    left_join(probs, by = "taxon_unit") %>%
    mutate(
      p_error = coalesce(p_error, 0),
      abundance = as.integer(round(abundance))
    ) %>%
    purrr::pmap_dfr(function(station, taxon_unit, abundance, target_taxon_unit, p_error, ...) {
      can_swap <- abundance > 0 &&
        !is.na(target_taxon_unit) &&
        is.finite(p_error) && p_error > 0
      
      if (!can_swap) {
        return(tibble(station = station, taxon_unit = taxon_unit, abundance = abundance))
      }
      
      n_swap <- rbinom(1L, abundance, pmin(p_error, 1))
      bind_rows(
        if (abundance - n_swap > 0) tibble(station = station, taxon_unit = taxon_unit, abundance = abundance - n_swap),
        if (n_swap > 0) tibble(station = station, taxon_unit = target_taxon_unit, abundance = n_swap)
      )
    }) %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop")
}


map_diagnostics <- function(species_long, map, probs) {
  source_abundance <- species_long %>%
    group_by(taxon_unit) %>%
    summarise(source_total_abundance = sum(abundance), .groups = "drop")

  total_abundance <- sum(source_abundance$source_total_abundance)

  diag <- source_abundance %>%
    left_join(
      map %>% select(source_taxon_unit, target_taxon_unit),
      by = c("taxon_unit" = "source_taxon_unit")
    ) %>%
    left_join(probs, by = "taxon_unit") %>%
    mutate(
      p_error = coalesce(p_error, 0),
      eligible = !is.na(target_taxon_unit) & p_error > 0
    )

  tibble(
    eligible_individual_share = if (total_abundance > 0) {
      sum(diag$source_total_abundance[diag$eligible]) / total_abundance
    } else NA_real_,
    expected_reassigned_pct = if (total_abundance > 0) {
      100 * sum(diag$source_total_abundance * diag$p_error * diag$eligible) / total_abundance
    } else NA_real_,
    n_eligible_source_species = sum(diag$eligible),
    n_source_species = nrow(diag)
  )
}

summarise_gradient <- function(data) {
  data %>%
    group_by(assemblage_id, mechanism, error_rate) %>%
    summarise(
      across(
        c(
          eligible_individual_share, expected_reassigned_pct,
          gamma_change_pct, mean_q0_change_pct, mean_q1_change_pct,
          mean_q2_change_pct, bray_mean_change_pct, sorensen_mean_change_pct,
          q0_stability, bray_stability, sorensen_stability
        ),
        list(
          median = ~ median(.x, na.rm = TRUE),
          p10 = ~ quantile(.x, .10, na.rm = TRUE, names = FALSE),
          p90 = ~ quantile(.x, .90, na.rm = TRUE, names = FALSE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n_iter = n(),
      .groups = "drop"
    ) %>%
    mutate(across(ends_with("_median") | ends_with("_p10") | ends_with("_p90"), ~ ifelse(is.nan(.x), NA_real_, .x)))
}

# ---- 2. Gradient simulations -------------------------------------------------
if (!dir.exists(file.path(INPUT_DIR, "assemblages"))) {
  stop("Prepared assemblies missing. Run 00_prepare_multitaxa_inputs_2024.R first.")
}

manifest_path <- file.path(INPUT_DIR, "assemblage_manifest.csv")
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE) %>%
  filter(assemblage_id %in% ASSEMBLAGES)

out_rows <- list()
index <- 1L

for (assemblage_id in manifest$assemblage_id) {
  message("Gradient analysis: ", assemblage_id)
  
  dat <- read_assemblage(assemblage_id)
  meta <- species_meta_from_lookup(dat$lookup)
  
  if (!nrow(dat$species) || !nrow(meta)) {
    message("  Skipped: no species-level community.")
    next
  }
  
  base_mat <- make_matrix(dat$species, dat$station_frame)
  base_alpha <- alpha_metrics(base_mat)
  base_gamma <- sum(colSums(base_mat) > 0)
  base_bray <- mean_dissimilarity(base_mat, "bray")
  base_sor <- mean_dissimilarity(base_mat, "sorensen")
  
  regional_pool <- read_regional_pool(assemblage_id)
  regional_meta <- meta %>%
    left_join(read_regional_source_audit(assemblage_id), by = "taxon_unit")
  expert_map <- read_expert_map(assemblage_id)
  
  for (iter_i in seq_len(N_SIM)) {
    # Crucial design: maps are drawn once per iteration and reused at every rate.
    observed_map <- make_coherent_map(meta, meta, seed = 100000L + iter_i)
    regional_map <- if (!is.null(regional_pool) && nrow(regional_pool)) {
      make_coherent_map(regional_meta, regional_pool, seed = 200000L + iter_i)
    } else NULL
    coherent_expert_map <- make_expert_map(meta, expert_map, seed = 300000L + iter_i)
    
    mechanisms <- list(
      observed_pool = observed_map,
      rare_weighted_observed_pool = observed_map
    )
    if (!is.null(regional_map)) mechanisms$regional_pool <- regional_map
    if (!is.null(coherent_expert_map) && nrow(coherent_expert_map)) mechanisms$expert_map <- coherent_expert_map
    
    for (rate in ERROR_RATES) {
      for (mechanism in names(mechanisms)) {
        probs <- if (mechanism == "rare_weighted_observed_pool") {
          make_rare_probs(dat$species, rate)
        } else {
          meta %>% transmute(taxon_unit, p_error = rate)
        }
        
        map_diag <- map_diagnostics(
          dat$species,
          mechanisms[[mechanism]],
          probs
        )

        scen_long <- simulate_map(
          dat$species, mechanisms[[mechanism]], probs,
          seed = 400000L + iter_i * 1000L + round(rate * 100) * 10L + match(mechanism, names(mechanisms))
        )
        scen_mat <- make_matrix(scen_long, dat$station_frame)
        scen_alpha <- alpha_metrics(scen_mat)
        scen_gamma <- sum(colSums(scen_mat) > 0)
        
        common <- intersect(base_alpha$station, scen_alpha$station)
        a0 <- base_alpha %>% filter(station %in% common) %>% arrange(match(station, common))
        a1 <- scen_alpha %>% filter(station %in% common) %>% arrange(match(station, common))
        
        row <- tibble(
          assemblage_id = assemblage_id,
          mechanism = mechanism,
          error_rate = rate,
          iter = iter_i,
          eligible_individual_share = map_diag$eligible_individual_share,
          expected_reassigned_pct = map_diag$expected_reassigned_pct,
          n_eligible_source_species = map_diag$n_eligible_source_species,
          n_source_species = map_diag$n_source_species,
          gamma_change_pct = 100 * (scen_gamma / base_gamma - 1),
          mean_q0_change_pct = 100 * (mean(a1$q0) / mean(a0$q0) - 1),
          mean_q1_change_pct = 100 * (mean(a1$q1) / mean(a0$q1) - 1),
          mean_q2_change_pct = 100 * (mean(a1$q2) / mean(a0$q2) - 1),
          q0_stability = safe_cor(a0$q0, a1$q0),
          bray_stability = distance_stability(base_mat, scen_mat, "bray"),
          sorensen_stability = distance_stability(base_mat, scen_mat, "sorensen"),
          bray_mean_change_pct = 100 * (mean_dissimilarity(scen_mat, "bray") / base_bray - 1),
          sorensen_mean_change_pct = 100 * (mean_dissimilarity(scen_mat, "sorensen") / base_sor - 1)
        )
        
        out_rows[[index]] <- row
        index <- index + 1L
      }
    }
  }
}

gradient_by_iter <- bind_rows(out_rows)
gradient_summary <- summarise_gradient(gradient_by_iter)

readr::write_csv(gradient_by_iter, file.path(OUT_DIR, "error_gradient_by_iter.csv"))
readr::write_csv(gradient_summary, file.path(OUT_DIR, "error_gradient_summary.csv"))

eligibility_summary <- gradient_summary %>%
  select(
    assemblage_id, mechanism, error_rate,
    eligible_individual_share_median, eligible_individual_share_p10, eligible_individual_share_p90,
    expected_reassigned_pct_median, expected_reassigned_pct_p10, expected_reassigned_pct_p90,
    n_eligible_source_species_median, n_source_species_median
  )
readr::write_csv(
  eligibility_summary,
  file.path(OUT_DIR, "error_gradient_eligibility_audit.csv")
)

# ---- 3. Appendix figures -----------------------------------------------------
assemblage_labels <- c(
  collembola_soil_core = "Collembola\nsoil cores",
  araneae_pitfall = "Araneae\npitfall traps",
  carabidae_pitfall = "Carabidae\npitfall traps",
  formicidae_pitfall = "Formicidae\npitfall traps",
  isopoda_pitfall = "Isopoda\npitfall traps",
  isopoda_hand_sorting = "Isopoda\nhand sorting",
  diplopoda_pitfall = "Diplopoda\npitfall traps",
  diplopoda_hand_sorting = "Diplopoda\nhand sorting"
)

mechanism_labels <- c(
  observed_pool = "Observed-pool congeneric",
  rare_weighted_observed_pool = "Rare-weighted observed-pool",
  regional_pool = "Regional-pool congeneric",
  expert_map = "Expert confusion map"
)

mechanism_colours <- c(
  observed_pool = "#0072B2",
  rare_weighted_observed_pool = "#CC79A7",
  regional_pool = "#D55E00",
  expert_map = "#009E73"
)

theme_paper <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      legend.position = "bottom",
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.major.y = element_blank()
    )
}

gradient_long_state <- gradient_summary %>%
  select(
    assemblage_id, mechanism, error_rate,
    gamma_change_pct_median, gamma_change_pct_p10, gamma_change_pct_p90,
    bray_mean_change_pct_median, bray_mean_change_pct_p10, bray_mean_change_pct_p90,
    sorensen_mean_change_pct_median, sorensen_mean_change_pct_p10, sorensen_mean_change_pct_p90
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, mechanism, error_rate),
    names_to = c("metric", ".value"),
    names_pattern = "(.*)_(median|p10|p90)"
  ) %>%
  mutate(
    metric = recode(
      metric,
      gamma_change_pct = "Regional richness (γ)",
      bray_mean_change_pct = "Mean Bray–Curtis",
      sorensen_mean_change_pct = "Mean Sørensen"
    ),
    assemblage = recode(assemblage_id, !!!assemblage_labels),
    mechanism = recode(mechanism, !!!mechanism_labels)
  )

p_state <- ggplot(
  gradient_long_state,
  aes(x = 100 * error_rate, y = median, colour = mechanism, fill = mechanism)
) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.3) +
  geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.6) +
  facet_grid(metric ~ assemblage, scales = "free_y") +
  scale_colour_manual(values = mechanism_colours, name = "Error mechanism") +
  scale_fill_manual(values = mechanism_colours, name = "Error mechanism") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  labs(
    title = "Response of inventory and beta diversity to error rates from 1% to 20%",
    subtitle = "Shaded bands show 10th–90th percentiles; each source→target map is fixed within an iteration across the full gradient.",
    x = "Identification-error rate (%)",
    y = "Change relative to the species baseline (%)"
  ) +
  theme_paper(base_size = 7.5) +
  theme(axis.text.x = element_text(angle = 0))

ggsave(file.path(OUT_DIR, "figures", "FigS_error_gradient_directional_effects.pdf"),
       p_state, width = 250, height = 165, units = "mm", device = grDevices::pdf)
ggsave(file.path(OUT_DIR, "figures", "FigS_error_gradient_directional_effects.png"),
       p_state, width = 250, height = 165, units = "mm", dpi = 450)

gradient_long_stability <- gradient_summary %>%
  select(
    assemblage_id, mechanism, error_rate,
    q0_stability_median, q0_stability_p10, q0_stability_p90,
    bray_stability_median, bray_stability_p10, bray_stability_p90,
    sorensen_stability_median, sorensen_stability_p10, sorensen_stability_p90
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, mechanism, error_rate),
    names_to = c("metric", ".value"),
    names_pattern = "(.*)_(median|p10|p90)"
  ) %>%
  mutate(
    metric = recode(
      metric,
      q0_stability = "q0 richness stability",
      bray_stability = "Bray–Curtis stability",
      sorensen_stability = "Sørensen stability"
    ),
    assemblage = recode(assemblage_id, !!!assemblage_labels),
    mechanism = recode(mechanism, !!!mechanism_labels)
  )

p_stability <- ggplot(
  gradient_long_stability,
  aes(x = 100 * error_rate, y = median, colour = mechanism, fill = mechanism)
) +
  geom_hline(yintercept = 1, colour = "grey35", linewidth = 0.3) +
  geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.6) +
  facet_grid(metric ~ assemblage, scales = "free_y") +
  scale_colour_manual(values = mechanism_colours, name = "Error mechanism") +
  scale_fill_manual(values = mechanism_colours, name = "Error mechanism") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  labs(
    title = "Stability of alpha and beta metrics along a 1%–20% identification-error gradient",
    subtitle = "Stability is the Spearman correlation with the species-level baseline.",
    x = "Identification-error rate (%)",
    y = "Stability"
  ) +
  theme_paper(base_size = 7.5)

ggsave(file.path(OUT_DIR, "figures", "FigS_error_gradient_stability.pdf"),
       p_stability, width = 250, height = 165, units = "mm", device = grDevices::pdf)
ggsave(file.path(OUT_DIR, "figures", "FigS_error_gradient_stability.png"),
       p_stability, width = 250, height = 165, units = "mm", dpi = 450)

message("\nDone.")
message("Tables: ", OUT_DIR)
message("Figures: ", file.path(OUT_DIR, "figures"))
message("Eligibility audit: ", file.path(OUT_DIR, "error_gradient_eligibility_audit.csv"))
