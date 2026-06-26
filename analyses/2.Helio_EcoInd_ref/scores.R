# ============================================================
# Delta scores vs reference locale (BioIndicateurs2)
# Input: analyses/Helio/delta_score_comp.csv
# Output: plot + tables
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(ggplot2)

in_path <- "analyses/Helio/delta_score_comp.csv"

# 1) Read
dat <- read_delim(in_path, delim = ";", show_col_types = FALSE)

# 2) Identify score columns and pivot long (index + framework)
score_cols <- names(dat)[str_detect(names(dat), "^score_.*_gamme_")]

long_scores <- dat %>%
  select(TTT, Nom_site, Annee, Traitement, Usage, Projet, all_of(score_cols)) %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = c("index", "framework"),
    names_pattern = "^score_(.*)_gamme_(.*)$",
    values_to = "score"
  ) %>%
  mutate(
    score = as.numeric(score),
    index = as.character(index),
    framework = as.character(framework)
  )

# 3) Compute delta vs local reference within each (Nom_site, Annee, index, framework)
delta_plot <- long_scores %>%
  group_by(Nom_site, Annee, index, framework) %>%
  mutate(score_ref = score[TTT == "ref"][1]) %>%
  ungroup() %>%
  mutate(delta = score - score_ref) %>%
  filter(TTT != "ref") %>%
  filter(!is.na(score_ref), !is.na(delta)) %>%
  filter(
    Nom_site %in% c("QualiAgro", "Thil", "Yvetot", "MetalEurope"),
    index %in% c("densi", "RSr")
  ) %>%
  mutate(
    index = case_when(
      index == "densi" ~ "Density",
      index == "RSr"   ~ "Species richness",
      TRUE ~ index
    ),
    abs_delta = abs(delta)
  )


# ordre des labels : Site -> TTT -> Traitement
row_order <- delta_plot %>%
  distinct(Nom_site, TTT, Traitement) %>%
  arrange(Nom_site, TTT, Traitement) 


# ============================================================
# 6) Best / co-best (sert à l’alpha)
# ============================================================

eps <- 0.01
case_vars <- c("Nom_site", "Traitement", "index")

rank_table <- delta_plot %>%
  filter(!is.na(abs_delta)) %>%
  group_by(across(all_of(case_vars))) %>%
  mutate(
    max_abs = max(abs_delta, na.rm = TRUE),
    is_best = abs(abs_delta - max_abs) <= eps,
    n_best  = sum(is_best, na.rm = TRUE),
    best_type = ifelse(is_best & n_best == 1, "best_solo",
                       ifelse(is_best & n_best > 1, "best_co", NA_character_))
  ) %>%
  ungroup() %>%
  mutate(is_best_any = best_type %in% c("best_solo", "best_co"))

# réinjecter le flag best dans l'objet de plot
delta_plot2 <- delta_plot %>%
  left_join(
    rank_table) %>%
  mutate(is_best_any = ifelse(is.na(is_best_any), FALSE, is_best_any))

# ============================================================
# 7) Forest plot (alpha = best/co-best uniquement)
# ============================================================
delta_plot2 <- delta_plot2 %>%
  group_by(Nom_site, index, Traitement) %>%
  mutate(panel_max = max(abs_delta, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Nom_site, index) %>%
  mutate(Traitement = forcats::fct_reorder(Traitement, panel_max, .fun = max, .desc = TRUE)) %>%
  ungroup()

pd <- position_dodge2(width = 0.7, preserve = "single")

p <- ggplot(
  delta_plot2,
  aes(
    x = abs_delta,
    y = Traitement,
    group = framework
  )
) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  
  # --- head
  geom_point(
    aes(
      colour = framework,
      shape  = is_best_any,
      size   = is_best_any,
      alpha  = is_best_any
    ),
    position = pd,
    stroke = 0.9
  ) +
  
    # --- stick
  geom_errorbar(
    aes(xmin = 0, xmax = pmax(0, abs_delta),
        colour = framework,
        alpha  = is_best_any
    ),
    linewidth = 0.6,
    position = pd
  )+

  facet_grid(Nom_site ~ index, scales = "free_y", space = "free_y") +
  
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 16),
    guide = "none"
  ) +
  scale_size_manual(
    values = c(`FALSE` = 2.0, `TRUE` = 3.0),
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c(`FALSE` = 0.35, `TRUE` = 1.0),
    guide = "none"
  ) +
  
  coord_cartesian(xlim = c(0, 5)) +
  theme_bw() +
  labs(
    x = "absolute Δ score\n(traitement − contrôle local)",
    y = NULL,
    colour = "Framework"
  ) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.y = element_text(size = 5),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p

ggsave(
  "analyses/Helio/forestplot_delta_scores_site_treatment.png",
  p, width = 14, height = 10, dpi = 200
)

# ============================================================
# 8) Summary best/worst (si tu veux garder, inchangé sauf alignement n_cases)
# ============================================================

best_counts <- rank_table %>%
  filter(!is.na(best_type)) %>%
  count(framework, best_type, name = "n")

n_cases <- rank_table %>%
  distinct(across(all_of(case_vars))) %>%
  nrow()

summary_rank <- best_counts %>%
  group_by(framework, best_type) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(pct = 100 * n / n_cases) %>%
  arrange(best_type, desc(n))

summary_rank

report_rank <- summary_rank %>%
  mutate(best_type = factor(best_type, levels = c("best_solo", "best_co"))) %>%
  tidyr::pivot_wider(names_from = best_type, values_from = c(n, pct), values_fill = 0)

report_rank
