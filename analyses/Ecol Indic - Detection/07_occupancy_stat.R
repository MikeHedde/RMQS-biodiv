# ============================================================
# STATS p̂
# Requiert: tab_p
# ============================================================

tab_p     <- readr::read_csv(file.path(out_dir, "occ_model_output/p_hat_by_method_unmarked_full.csv"), show_col_types = FALSE)

# Data handling
## Sécurisation des p_hat
tab_p2 <- tab_p %>%
  mutate(
    p_hat_clip = pmin(pmax(p_hat, 1e-5), 1 - 1e-5),
    logit_p    = qlogis(p_hat_clip),
    # SE approx à partir des IC95%
    se_p       = (ucl - lcl) / (2 * 1.96),
    se_logit   = se_p / (p_hat_clip * (1 - p_hat_clip)),
    weight     = ifelse(is.finite(se_logit) & se_logit > 0,
                        1 / se_logit^2, 1)
  )

## Ordre des méthodes
tab_p2$method <- factor(
  tab_p2$method,
  levels = c("Pitfall2", "Pitfall4", "Pitfall6",
             "Pitfall8", "Pitfall10", "GPD")
)

## Variable numérique pour le gradient Pitfall
tab_p2 <- tab_p2 %>%
  mutate(
    ntrap = case_when(
      method == "Pitfall2"  ~ 2,
      method == "Pitfall4"  ~ 4,
      method == "Pitfall6"  ~ 6,
      method == "Pitfall8"  ~ 8,
      method == "Pitfall10" ~ 10,
      TRUE ~ NA_real_
    ),
    is_gpd = method == "GPD"
  )

# Tests stats
## Effet global de la méthode
m0 <- lmer(logit_p ~ 1 + (1 | species),
           data = tab_p2, weights = weight)

m1 <- lmer(logit_p ~ method + (1 | species),
           data = tab_p2, weights = weight)

anova(m0, m1) 

##Test de tendance Pitfall2 → Pitfall10
m_trend <- lmer(logit_p ~ ntrap + (1 | species),
                data = filter(tab_p2, !is_gpd),
                weights = weight)

summary(m_trend)

##GPD vs Pitfall (à effort contrôlé)
m_gpd <- lmer(logit_p ~ method + (1 | species),
              data = filter(tab_p2,
                            method %in% c("Pitfall6", "GPD")),
              weights = weight)

summary(m_gpd)

##Contrastes emmeans + table “différences de méthodes”
emm <- emmeans(m1, ~ method)
emm

contr <- contrast(
  emm,
  method = list(
    "Pitfall4 vs Pitfall2"   = c(-1, 1, 0, 0, 0, 0),
    "Pitfall6 vs Pitfall4"   = c(0, -1, 1, 0, 0, 0),
    "Pitfall8 vs Pitfall6"   = c(0, 0, -1, 1, 0, 0),
    "Pitfall10 vs Pitfall8"  = c(0, 0, 0, -1, 1, 0),
    "GPD vs Pitfall4"        = c(0, -1, 0, 0, 0, 1),
    "GPD vs Pitfall6"        = c(0, 0, -1, 0, 0, 1),
    "GPD vs Pitfall8"        = c(0, 0, 0, -1, 0, 1)
  )
)

contr_res <- summary(contr, infer = TRUE)
contr_res

tab_results <- contr_res %>%
  as.data.frame() %>%
  mutate(
    odds_ratio = exp(estimate),
    lower_OR   = exp(lower.CL),
    upper_OR   = exp(upper.CL)
  ) %>%
  select(contrast, estimate, lower.CL, upper.CL,
         odds_ratio, lower_OR, upper_OR, p.value)

tab_results

dir.create(
  file.path(out_dir, "occ_model_output", "detection_stats"),
  recursive = TRUE,
  showWarnings = FALSE
)
readr::write_csv(tab_results, file.path(out_dir, "occ_model_output/detection_stats/p_hat_by_method_unmarked_full.csv"))
