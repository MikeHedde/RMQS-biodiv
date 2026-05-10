# ============================================================
# 09_traits_detectability_refond.R
# Traits -> detectability from occupancy-derived p-hat
#
# Objectives
#   1) Build a single coherent trait-detectability dataset from:
#      - occupancy predictions: p_hat_by_method_unmarked_full.csv
#      - raw BETSI-like Carabidae trait table: trait_carabidae.csv
#   2) Fit statistical models consistent with the figures:
#      - Body length: method + log10(body length) vs method * log10(body length)
#      - Wing development: method + log10(body length) + wing vs method * wing + log10(body length)
#      - Optional cross-trait model: log_BL * full_wing, controlling for method
#   3) Produce figures based on the fitted models, not just raw summaries.
#
# Required upstream object
#   out_dir, normally defined by previous scripts.
#   If absent, default below is used.
# ============================================================

# ------------------------------------------------------------
# 0) Packages and paths
# ------------------------------------------------------------
occ_p_file   <- file.path(out_dir, "occ_model_output", "p_hat_by_method_unmarked_full.csv")
trait_file   <- "data/raw-data/1.faune/trait_carabidae.csv"
trait_out_dir <- file.path(out_dir, "trait")
dir.create(trait_out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(occ_p_file)) stop("Cannot find: ", occ_p_file)
if (!file.exists(trait_file)) stop("Cannot find: ", trait_file)

method_levels <- c("Pitfall2", "Pitfall4", "Pitfall6", "Pitfall8", "Pitfall10", "GPD")

# Optional switches
use_weight_cap <- TRUE
weight_cap_prob <- 0.99

# ------------------------------------------------------------
# 1) Read occupancy p-hat and raw trait table
# ------------------------------------------------------------

tab_p <- readr::read_csv(occ_p_file, show_col_types = FALSE)
trait0 <- read.csv(trait_file, header = TRUE, sep = ";", stringsAsFactors = FALSE)

required_p_cols <- c("species", "method", "p_hat", "lcl", "ucl")
missing_p_cols <- setdiff(required_p_cols, names(tab_p))
if (length(missing_p_cols) > 0) {
  stop("Missing columns in p-hat table: ", paste(missing_p_cols, collapse = ", "))
}

required_trait_cols <- c("species", "raw_trait_value", "trait_value", "attribute_trait")
missing_trait_cols <- setdiff(required_trait_cols, names(trait0))
if (length(missing_trait_cols) > 0) {
  stop("Missing columns in trait table: ", paste(missing_trait_cols, collapse = ", "))
}

# ------------------------------------------------------------
# 2) Canonical species-name map from trait table
# ------------------------------------------------------------

dir.create(
  file.path(out_dir, "trait", "reframed_trait_datasets"),
  recursive = TRUE,
  showWarnings = FALSE
)

message("Parsing trait-table species names...")

tmp <- trait0 %>%
  dplyr::distinct(species) %>%
  dplyr::mutate(parsed = rgnparser::gn_parse(species))

tmp2 <- tmp %>%
  tidyr::unnest_wider(parsed)

if (!"canonical" %in% names(tmp2)) {
  stop("rgnparser output does not contain a 'canonical' column. Inspect tmp2.")
}

map_species <- tmp2 %>%
  dplyr::transmute(
    species,
    sp = purrr::map_chr(canonical, \(x) {
      if (is.null(x) || length(x) == 0) return(NA_character_)
      if (!is.null(x$simple)  && length(x$simple)  > 0 && nzchar(x$simple[1]))  return(x$simple[1])
      if (!is.null(x$full)    && length(x$full)    > 0 && nzchar(x$full[1]))    return(x$full[1])
      if (!is.null(x$stemmed) && length(x$stemmed) > 0 && nzchar(x$stemmed[1])) return(x$stemmed[1])
      NA_character_
    })
  )

bad_parsed <- tmp2 %>%
  dplyr::filter(purrr::map_lgl(canonical, ~ is.null(.x) || length(.x) == 0))

readr::write_csv(map_species, file.path(trait_out_dir, "reframed_trait_datasets/trait_species_name_map.csv"))

# ------------------------------------------------------------
# 3) Body length and wing-development trait tables
# ------------------------------------------------------------

BL <- trait0 %>%
  dplyr::left_join(map_species, by = "species") %>%
  dplyr::filter(raw_trait_value == "Body_length") %>%
  dplyr::mutate(trait_value = as.numeric(trait_value)) %>%
  dplyr::group_by(sp) %>%
  dplyr::summarise(mean_BL = mean(trait_value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::rename(species = sp) %>%
  dplyr::filter(!is.na(species), species != "", is.finite(mean_BL), mean_BL > 0)

trait_wing_sp <- trait0 %>%
  dplyr::left_join(map_species, by = "species") %>%
  dplyr::filter(raw_trait_value == "Wing_development") %>%
  dplyr::mutate(
    wing_bin = dplyr::case_when(
      stringr::str_detect(attribute_trait, stringr::regex("macropterous", ignore_case = TRUE)) ~ "fully_winged",
      stringr::str_detect(attribute_trait, stringr::regex("brachypterous|apterous", ignore_case = TRUE)) ~ "not_winged",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::group_by(sp) %>%
  dplyr::summarise(
    n_records = sum(!is.na(wing_bin)),
    n_winged = sum(wing_bin == "fully_winged", na.rm = TRUE),
    n_notwinged = sum(wing_bin == "not_winged", na.rm = TRUE),
    pct_winged = dplyr::if_else(n_records > 0, n_winged / n_records, NA_real_),
    pct_notwinged = dplyr::if_else(n_records > 0, n_notwinged / n_records, NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::rename(species = sp) %>%
  dplyr::select(species, pct_winged, n_records, n_winged, n_notwinged)

trait_wing_gn <- trait_wing_sp %>%
  dplyr::mutate(genus = stringr::word(species, 1)) %>%
  dplyr::group_by(genus) %>%
  dplyr::summarise(
    pct_winged_gn = if (all(is.na(pct_winged))) NA_real_ else stats::median(pct_winged, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(BL, file.path(trait_out_dir, "reframed_trait_datasets/trait_body_length_species.csv"))
readr::write_csv(trait_wing_sp, file.path(trait_out_dir, "reframed_trait_datasets/trait_wing_species.csv"))
readr::write_csv(trait_wing_gn, file.path(trait_out_dir, "reframed_trait_datasets/trait_wing_genus_imputation.csv"))

# ------------------------------------------------------------
# 4) Build common analysis dataset
# ------------------------------------------------------------

# We use the logit-scale CI to derive SEs, which is more coherent with logit_p models.
# If CI bounds are exactly 0/1, they are clipped before transformation.

dat_traits <- tab_p %>%
  dplyr::mutate(
    species = as.character(species),
    method = factor(method, levels = method_levels)
  ) %>%
  dplyr::left_join(BL, by = "species") %>%
  dplyr::left_join(trait_wing_sp, by = "species") %>%
  dplyr::mutate(genus = stringr::word(species, 1)) %>%
  dplyr::left_join(trait_wing_gn, by = "genus") %>%
  dplyr::mutate(
    wing_imputed_genus = is.na(pct_winged) & !is.na(pct_winged_gn),
    pct_winged_final = dplyr::coalesce(pct_winged, pct_winged_gn),
    full_wing = dplyr::case_when(
      pct_winged_final == 1 ~ "fully_winged",
      pct_winged_final <  1 ~ "not_fully_winged",
      TRUE ~ NA_character_
    ),
    full_wing = factor(full_wing, levels = c("not_fully_winged", "fully_winged")),
    log_BL = log10(mean_BL),
    p_clip = pmin(pmax(p_hat, 1e-5), 1 - 1e-5),
    lcl_clip = pmin(pmax(lcl, 1e-5), 1 - 1e-5),
    ucl_clip = pmin(pmax(ucl, 1e-5), 1 - 1e-5),
    logit_p = qlogis(p_clip),
    se_logit = (qlogis(ucl_clip) - qlogis(lcl_clip)) / (2 * 1.96),
    weight_raw = 1 / se_logit^2,
    species = factor(species)
  ) %>%
  dplyr::mutate(
    weight = dplyr::case_when(
      is.finite(weight_raw) & weight_raw > 0 ~ weight_raw,
      TRUE ~ 1
    )
  )

if (use_weight_cap) {
  cap <- stats::quantile(dat_traits$weight[is.finite(dat_traits$weight)], weight_cap_prob, na.rm = TRUE)
  dat_traits <- dat_traits %>% dplyr::mutate(weight = pmin(weight, cap))
}

# Separate datasets are useful because body length can be available when wing is missing.
dat_BL <- dat_traits %>%
  dplyr::filter(
    !is.na(method),
    !is.na(log_BL), is.finite(log_BL),
    is.finite(logit_p),
    is.finite(weight), weight > 0
  ) %>%
  dplyr::mutate(method = relevel(method, ref = "Pitfall4"))

dat_wing <- dat_traits %>%
  dplyr::filter(
    !is.na(method),
    !is.na(log_BL), is.finite(log_BL),
    !is.na(full_wing),
    is.finite(logit_p),
    is.finite(weight), weight > 0
  ) %>%
  dplyr::mutate(
    method = relevel(method, ref = "Pitfall4"),
    full_wing = relevel(full_wing, ref = "not_fully_winged")
  )

readr::write_csv(dat_traits, file.path(trait_out_dir, "reframed_trait_datasets/traits_detectability_dataset_all.csv"))
readr::write_csv(dat_BL, file.path(trait_out_dir, "reframed_trait_datasets/traits_detectability_dataset_body_length.csv"))
readr::write_csv(dat_wing, file.path(trait_out_dir, "reframed_trait_datasets/traits_detectability_dataset_wing.csv"))

summary_traits <- tibble::tibble(
  dataset = c("all_joined", "body_length", "wing"),
  n_rows = c(nrow(dat_traits), nrow(dat_BL), nrow(dat_wing)),
  n_species = c(dplyr::n_distinct(dat_traits$species), dplyr::n_distinct(dat_BL$species), dplyr::n_distinct(dat_wing$species)),
  n_species_with_genus_wing_imputation = c(
    dplyr::n_distinct(dat_traits$species[dat_traits$wing_imputed_genus], na.rm = TRUE),
    dplyr::n_distinct(dat_BL$species[dat_BL$wing_imputed_genus], na.rm = TRUE),
    dplyr::n_distinct(dat_wing$species[dat_wing$wing_imputed_genus], na.rm = TRUE)
  )
)
readr::write_csv(summary_traits, file.path(trait_out_dir, "reframed_trait_datasets/traits_detectability_dataset_summary.csv"))
print(summary_traits)

# ------------------------------------------------------------
# 5) Utility functions
# ------------------------------------------------------------

save_txt <- function(..., file) {
  capture.output(..., file = file)
}

save_varcorr <- function(model, filename) {
  vc <- as.data.frame(lme4::VarCorr(model))
  readr::write_csv(vc, file.path(trait_out_dir, filename))
  invisible(vc)
}


dir.create(
  file.path(out_dir, "trait", "models"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 6) Body-length models
# ------------------------------------------------------------

m_BL_add_ML <- lmer(
  logit_p ~ method + log_BL + (1 | species),
  data = dat_BL,
  weights = weight,
  REML = FALSE
)

m_BL_int_ML <- lmer(
  logit_p ~ method * log_BL + (1 | species),
  data = dat_BL,
  weights = weight,
  REML = FALSE
)

BL_model_LRT <- anova(m_BL_add_ML, m_BL_int_ML)

m_BL_add <- update(m_BL_add_ML, REML = TRUE)
m_BL_int <- update(m_BL_int_ML, REML = TRUE)

BL_anova_add <- car::Anova(m_BL_add, type = 3)
BL_anova_int <- car::Anova(m_BL_int, type = 3)

BL_trends <- emmeans::emtrends(m_BL_int, ~ method, var = "log_BL")
BL_trends_pairs <- pairs(BL_trends, adjust = "tukey")

save_txt(
  summary(m_BL_add),
  BL_anova_add,
  BL_model_LRT,
  file = file.path(trait_out_dir, "models/BL_additive_summary.txt")
)

save_txt(
  summary(m_BL_int),
  BL_anova_int,
  BL_model_LRT,
  BL_trends,
  BL_trends_pairs,
  file = file.path(trait_out_dir, "models/BL_method_interaction_summary.txt")
)

readr::write_csv(as.data.frame(BL_model_LRT), file.path(trait_out_dir, "models/BL_additive_vs_interaction_LRT.csv"))
readr::write_csv(as.data.frame(BL_anova_add), file.path(trait_out_dir, "models/BL_additive_Anova_typeIII.csv"))
readr::write_csv(as.data.frame(BL_anova_int), file.path(trait_out_dir, "models/BL_interaction_Anova_typeIII.csv"))
readr::write_csv(as.data.frame(BL_trends), file.path(trait_out_dir, "models/BL_emtrends_by_method.csv"))
readr::write_csv(as.data.frame(BL_trends_pairs), file.path(trait_out_dir, "models/BL_emtrends_pairwise.csv"))
save_varcorr(m_BL_int, "models/BL_interaction_random_effect_variance.csv")

# ------------------------------------------------------------
# 7) Wing-development models
# ------------------------------------------------------------

m_wing_add_ML <- lmer(
  logit_p ~ method + log_BL + full_wing + (1 | species),
  data = dat_wing,
  weights = weight,
  REML = FALSE
)

m_wing_int_ML <- lmer(
  logit_p ~ method * full_wing + log_BL + (1 | species),
  data = dat_wing,
  weights = weight,
  REML = FALSE
)

wing_model_LRT <- anova(m_wing_add_ML, m_wing_int_ML)

m_wing_add <- update(m_wing_add_ML, REML = TRUE)
m_wing_int <- update(m_wing_int_ML, REML = TRUE)

wing_anova_add <- car::Anova(m_wing_add, type = 3)
wing_anova_int <- car::Anova(m_wing_int, type = 3)

wing_emm <- emmeans::emmeans(
  m_wing_int,
  ~ full_wing | method,
  at = list(log_BL = stats::median(dat_wing$log_BL, na.rm = TRUE))
)
wing_contrasts <- pairs(wing_emm, adjust = "none")
wing_effect_by_method <- contrast(
  wing_emm,
  method = "revpairwise",
  by = "method",
  adjust = "none"
) %>%
  as.data.frame()

save_txt(
  summary(m_wing_add),
  wing_anova_add,
  wing_model_LRT,
  file = file.path(trait_out_dir, "models/wing_additive_summary.txt")
)

save_txt(
  summary(m_wing_int),
  wing_anova_int,
  wing_model_LRT,
  wing_emm,
  wing_contrasts,
  file = file.path(trait_out_dir, "models/wing_method_interaction_summary.txt")
)

readr::write_csv(as.data.frame(wing_model_LRT), file.path(trait_out_dir, "models/wing_additive_vs_interaction_LRT.csv"))
readr::write_csv(as.data.frame(wing_anova_add), file.path(trait_out_dir, "models/wing_additive_Anova_typeIII.csv"))
readr::write_csv(as.data.frame(wing_anova_int), file.path(trait_out_dir, "models/wing_interaction_Anova_typeIII.csv"))
readr::write_csv(as.data.frame(wing_emm), file.path(trait_out_dir, "models/wing_emmeans_by_method.csv"))
readr::write_csv(as.data.frame(wing_contrasts), file.path(trait_out_dir, "models/wing_contrasts_by_method.csv"))
readr::write_csv(wing_effect_by_method, file.path(trait_out_dir, "models/wing_effect_by_method.csv"))
save_varcorr(m_wing_int, "models/wing_interaction_random_effect_variance.csv")

# ------------------------------------------------------------
# 8) Optional cross-trait model: body length × wing development
#    This tests whether the body-size effect differs between wing categories,
#    while controlling for method. It is not the main protocol-filtering model.
# ------------------------------------------------------------

m_cross_add_ML <- lmer(
  logit_p ~ method + log_BL + full_wing + (1 | species),
  data = dat_wing,
  weights = weight,
  REML = FALSE
)

m_cross_int_ML <- lmer(
  logit_p ~ method + log_BL * full_wing + (1 | species),
  data = dat_wing,
  weights = weight,
  REML = FALSE
)

cross_model_LRT <- anova(m_cross_add_ML, m_cross_int_ML)
m_cross_int <- update(m_cross_int_ML, REML = TRUE)
cross_anova_int <- car::Anova(m_cross_int, type = 3)

cross_trends <- emmeans::emtrends(m_cross_int, ~ full_wing, var = "log_BL")
cross_trends_pairs <- pairs(cross_trends, adjust = "none")

save_txt(
  summary(m_cross_int),
  cross_anova_int,
  cross_model_LRT,
  cross_trends,
  cross_trends_pairs,
  file = file.path(trait_out_dir, "models/cross_BL_by_wing_summary.txt")
)
readr::write_csv(as.data.frame(cross_model_LRT), file.path(trait_out_dir, "models/cross_BL_by_wing_LRT.csv"))
readr::write_csv(as.data.frame(cross_trends), file.path(trait_out_dir, "models/cross_BL_trends_by_wing.csv"))
readr::write_csv(as.data.frame(cross_trends_pairs), file.path(trait_out_dir, "models/cross_BL_trends_by_wing_pairs.csv"))


dir.create(
  file.path(out_dir, "trait", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)
# ------------------------------------------------------------
# 9) Figure A: body length effect by method
#     Uses predictions from m_BL_int, i.e. same model as BL emtrends.
# ------------------------------------------------------------
method_levels <- c("Pitfall2", "Pitfall4", "Pitfall6", "Pitfall8", "Pitfall10", "GPD")

log_seq <- seq(
  stats::quantile(dat_BL$log_BL, 0.05, na.rm = TRUE),
  stats::quantile(dat_BL$log_BL, 0.95, na.rm = TRUE),
  length.out = 100
)

pred_BL <- emmeans::emmeans(
  m_BL_int,
  ~ method | log_BL,
  at = list(log_BL = log_seq)
) %>%
  as.data.frame() %>%
  dplyr::mutate(
    p_fit = plogis(emmean),
    p_lcl = plogis(lower.CL),
    p_ucl = plogis(upper.CL),
    method = factor(method, levels = method_levels)
  ) %>%
  dplyr::arrange(method, log_BL)

dat_BL <- dat_BL %>%
  mutate(method = factor(method, levels = method_levels))

pred_BL <- pred_BL %>%
  mutate(method = factor(method, levels = method_levels))

p_BL <- ggplot(dat_BL, aes(x = log_BL, y = p_clip)) +
  geom_point(alpha = 0.5, size = 1.7) +
  geom_ribbon(
    data = pred_BL,
    aes(x = log_BL, ymin = p_lcl, ymax = p_ucl),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_BL,
    aes(x = log_BL, y = p_fit),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ method, ncol = 6) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "log10(body length, mm)", y = "Detectability p̂") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 22)
  )

print(p_BL)

ggsave(file.path(trait_out_dir, "figures/Detectability_vs_BL_by_method.png"), p_BL, width = 9, height = 6, dpi = 200)
ggsave(file.path(trait_out_dir, "figures/Fig_trait_BL_detectability_by_method.png"), p_BL, width = 12, height = 5, dpi = 300)

# ------------------------------------------------------------
# 10) Figure B: wing development by method
#     Raw points/violins + model-based predictions from m_wing_int
#     at median body length.
# ------------------------------------------------------------

pred_wing <- emmeans::emmeans(
  m_wing_int,
  ~ method * full_wing,
  at = list(log_BL = stats::median(dat_wing$log_BL, na.rm = TRUE))
) %>%
  as.data.frame() %>%
  dplyr::mutate(
    p_fit = plogis(emmean),
    p_lcl = plogis(lower.CL),
    p_ucl = plogis(upper.CL),
    method = factor(method, levels = method_levels),
    full_wing_plot = factor(
      full_wing,
      levels = c("not_fully_winged", "fully_winged"),
      labels = c("Not fully", "Full")
    )
  )

dat_wing_plot <- dat_wing %>%
  dplyr::mutate(
    full_wing_plot = factor(
      full_wing,
      levels = c("not_fully_winged", "fully_winged"),
      labels = c("Not fully", "Full")
    )
  )

dat_wing_plot <- dat_wing_plot %>%
  mutate(method = factor(method, levels = method_levels))

p_wing <- ggplot(dat_wing_plot, aes(x = full_wing_plot, y = p_clip)) +
  geom_violin(alpha = 0.45) +
  geom_jitter(alpha = 0.5, width = 0.15, size = 1.7) +
  geom_pointrange(
    data = pred_wing,
    aes(x = full_wing_plot, y = p_fit, ymin = p_lcl, ymax = p_ucl),
    inherit.aes = FALSE,
    linewidth = 0.7,
    size = 0.8
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ method, ncol = 6) +
  labs(x = "Wing development", y = "Detectability p̂") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 22),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(p_wing)

ggsave(file.path(trait_out_dir, "figures/detectability_vs_wing_by_method.png"), p_wing, width = 12, height = 6, dpi = 200)
ggsave(file.path(trait_out_dir, "figures/Fig_trait_wing_binary_detectability_by_method.png"), p_wing, width = 12, height = 5, dpi = 300)

# ------------------------------------------------------------
# 11) Combined figure
# ------------------------------------------------------------

p_traits <- cowplot::plot_grid(
  p_BL,
  p_wing,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1),
  labels = c("A", "B")
)

print(p_traits)

ggsave(
  file.path(trait_out_dir, "figures/detectability_traits_by_method.png"),
  p_traits,
  width = 28,
  height = 12,
  dpi = 200
)

ggsave(
  file.path(trait_out_dir, "figures/Fig_traits_detectability_combined.png"),
  p_traits,
  width = 12,
  height = 10,
  dpi = 300
)

# ------------------------------------------------------------
# 12) Compact text summary for manuscript/revision notes
# ------------------------------------------------------------

sink(file.path(trait_out_dir, "trait_detectability_key_results.txt"))
cat("TRAIT-DETECTABILITY ANALYSIS\n")
cat("============================\n\n")
cat("Dataset summary\n")
print(summary_traits)
cat("\nBody length: additive vs method interaction LRT\n")
print(BL_model_LRT)
cat("\nBody length: Type III Anova for interaction model\n")
print(BL_anova_int)
cat("\nBody length: method-specific slopes\n")
print(BL_trends)
cat("\nBody length: pairwise slope contrasts\n")
print(BL_trends_pairs)
cat("\nWing development: additive vs method interaction LRT\n")
print(wing_model_LRT)
cat("\nWing development: Type III Anova for interaction model\n")
print(wing_anova_int)
cat("\nWing development: estimated means by method at median body length\n")
print(wing_emm)
cat("\nWing development: contrasts within method\n")
print(wing_contrasts)
cat("\nCross-trait model: log_BL * full_wing LRT\n")
print(cross_model_LRT)
sink()

message("Trait-detectability analysis completed. Outputs written to: ", trait_out_dir)
