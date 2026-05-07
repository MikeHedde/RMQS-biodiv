# ============================================================
# TRAITS → DÉTECTABILITÉ (p̂) ~ méthode * log(length)
# Sorties: modèles, figures
# ============================================================

tab_p <- readr::read_csv(file.path(out_dir, "occ_model_output/p_hat_by_method_unmarked_full.csv"), show_col_types = FALSE)

trait0 <- read.csv("data/raw-data/1.faune/trait_carabidae.csv", h=T, sep=";")

# Max Body length
tmp <- trait0 %>%
  distinct(species) %>%
  mutate(parsed = rgnparser::gn_parse(species))

# Aplatit parsed en colonnes
tmp2 <- tmp %>%
  tidyr::unnest_wider(parsed)

canon_col <- intersect(
  c("canonicalName", "canonical", "canonical_name", "canonicalname"),
  names(tmp2)
)[1]

map <- tmp2 %>%
  transmute(
    species,
    sp = map_chr(canonical, \(x) {
      if (is.null(x) || length(x) == 0) return(NA_character_)
      if (!is.null(x$simple)  && length(x$simple)  > 0 && nzchar(x$simple[1])) return(x$simple[1])
      if (!is.null(x$full)    && length(x$full)    > 0 && nzchar(x$full[1]))   return(x$full[1])
      if (!is.null(x$stemmed) && length(x$stemmed) > 0 && nzchar(x$stemmed[1]))return(x$stemmed[1])
      NA_character_
    })
  )

bad <- tmp2 %>% filter(map_lgl(canonical, ~ is.null(.x) || length(.x) == 0))
bad$species %>% head(50)

BL <- trait0 %>%
  left_join(map, by = "species") %>%
  filter(raw_trait_value == "Body_length") %>%
  mutate(trait_value = as.numeric(trait_value)) %>%
  group_by(sp) %>%
  summarise(mean_BL = mean(trait_value, na.rm = TRUE), .groups = "drop") %>%
  rename(species = sp) %>%
  filter(!is.na(species), species != "")

dat <- tab_p %>% left_join(BL, by = "species")

# Mise à l’échelle logit + poids méta depuis IC95%
dat2 <- dat %>%
  mutate(
    p_clip  = pmin(pmax(p_hat, 1e-5), 1 - 1e-5),
    logit_p = logit(p_clip),
    se_p    = ifelse(!is.na(lcl) & !is.na(ucl), (ucl - lcl) / (2 * 1.96), NA_real_),
    se_logit = ifelse(!is.na(se_p), se_p / (p_clip * (1 - p_clip)), NA_real_),
    weight   = ifelse(is.finite(1/(se_logit^2)), 1/(se_logit^2), 1),
    log_BL = log10(mean_BL)
  ) %>%
  filter(!is.na(log_BL), is.finite(log_BL)) %>%
  mutate(species = factor(species))

dat2_tr <- dat2 %>%
  filter(method != "Pitfall10") %>%
  mutate(
    method = factor(
      method,
      levels = c("Pitfall2","Pitfall4","Pitfall6","Pitfall8","GPD")
    )
  )

dat2_tr$method <- relevel(dat2_tr$method, ref = "Pitfall4")

m_tr_int <- lmer(
  logit_p ~ method * log_BL + (1 | species),
  data = dat2_tr,
  weights = weight,
  REML = TRUE
)

s <- summary(m_tr_int)

trends <- emmeans::emtrends(m_tr_int, ~ method, var = "log_BL")
pairs(trends)

capture.output(s, file = file.path(out_dir, "trait/model_traits_detectability_summary.txt"))
capture.output(trends, file = file.path(out_dir, "trait/model_traits_emtrends.txt"))
capture.output(pairs(trends), file = file.path(out_dir, "trait/model_traits_emtrends_pairs.txt"))

# Prédictions lissées sur gamme de masses par méthode
# grid global (évite les plages différentes par méthode)
log_seq <- seq(
  quantile(dat2_tr$log_BL, 0.05, na.rm = TRUE),
  quantile(dat2_tr$log_BL, 0.95, na.rm = TRUE),
  length.out = 100
)

emm_grid <- emmeans::emmeans(
  m_tr_int,
  ~ method | log_BL,
  at = list(log_BL = log_seq)
)

pred_df <- as.data.frame(emm_grid) %>%
  dplyr::mutate(
    p_fit = plogis(emmean),
    p_lcl = plogis(lower.CL),
    p_ucl = plogis(upper.CL)
  )

dat2_tr$method <- factor(dat2_tr$method, levels = c("Pitfall2", "Pitfall4", "Pitfall6", "Pitfall8", "GPD"))
pred_df$method <- factor(pred_df$method, levels = c("Pitfall2", "Pitfall4", "Pitfall6", "Pitfall8", "GPD"))

p_BL <- ggplot(dat2_tr, aes(x = log_BL, y = p_clip, colour = method)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = pred_df,
              aes(x = log_BL, ymin = p_lcl, ymax = p_ucl, fill = method),
              alpha = 0.15, colour = NA, inherit.aes = FALSE) +
  geom_line(data = pred_df, aes(x = log_BL, y = p_fit), linewidth = 1) +
  facet_wrap(~ method, ncol = 5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "log10(body length, mm)", y = "Detectability p̂") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 22))

ggsave(file.path(out_dir, "trait/Detectability_vs_BL_by_method.png"), p_BL, width = 9, height = 6, dpi = 200)


# Motion strategy
trait_sp <- trait0 %>%
  left_join(map, by = "species") %>%
  filter(raw_trait_value == "Wing_development") %>%
  mutate(wing_bin = case_when(
    str_detect(attribute_trait, regex("macropterous", ignore_case = TRUE)) ~ "fully_winged",
    str_detect(attribute_trait, regex("brachypterous|apterous", ignore_case = TRUE)) ~ "not_winged",
    TRUE ~ NA_character_
  ))%>%
  group_by(sp) %>%
  summarise(
    n_records   = sum(!is.na(wing_bin)),
    n_winged    = sum(wing_bin == "fully_winged", na.rm = TRUE),
    n_notwinged = sum(wing_bin == "not_winged",  na.rm = TRUE),
    pct_winged  = n_winged / n_records,
    pct_notwinged = n_notwinged / n_records
  ) %>%
  ungroup() %>%
  rename(species = sp) %>%
  select(species, pct_winged) 

trait_gn <- trait_sp %>%
  mutate(genus = word(species, 1)) %>%
  group_by(genus) %>%
  summarise(pct_winged_gn  = median(pct_winged)) %>%
  ungroup()
  
dat <- tab_p %>% left_join(trait_sp, by = "species") %>%
  mutate(genus = word(species, 1)) 

dat_filled <- dat %>%
  left_join(trait_gn, by = "genus") %>%
  mutate(
    pct_winged = ifelse(is.na(pct_winged), pct_winged_gn, pct_winged)
  ) %>%
  select(-pct_winged_gn)%>%
  mutate(full_wing = case_when(
    pct_winged == 1 ~ "fully_winged",
    pct_winged < 1 ~ "not_fully_winged"
  ))

dat_filled$method <- factor(dat_filled$method, levels = c("Pitfall2", "Pitfall4", "Pitfall6", "Pitfall8", "Pitfall10", "GPD"))

ggplot(subset(dat_filled, !is.na(full_wing) & method != "Pitfall10"), aes(x = full_wing, y = p_hat, colour = method)) +
  geom_violin()+
  geom_jitter(alpha = 0.5, width = 0.15) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(method~., ncol = 6) +
  labs(x = "Wing development", y = "Detectability p̂",
       #title = "Effect of wing development on detectability by method"
       ) +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 22)) ->
  p_wing

ggsave(file.path(out_dir, "trait/detectability_vs_wing_by_method.png"), p_wing, width = 12, height = 6, dpi = 200)

p_traits <- plot_grid(
  p_BL,
  p_wing,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

p_traits

ggsave(
  file.path(out_dir, "trait/detectability_traits_by_method.png"),
  p_traits,
  width = 28,
  height = 12,
  dpi = 200
)

#########################""
# croisement des deux traits

dat_int <- tab_p %>%
  left_join(BL,    by = "species") %>%      # max_BL
  left_join(trait_sp, by = "species") %>%      # pct_winged
  mutate(
    log_BL = log10(mean_BL),
    full_wing = case_when(
      pct_winged == 1 ~ "fully_winged",
      pct_winged <  1 ~ "not_fully_winged",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(log_BL), !is.na(full_wing)) %>%
  mutate(
    p_clip  = pmin(pmax(p_hat, 1e-5), 1 - 1e-5),
    logit_p = qlogis(p_clip),
    se_p    = ifelse(!is.na(lcl) & !is.na(ucl),
                     (ucl - lcl) / (2 * 1.96), NA_real_),
    se_logit = se_p / (p_clip * (1 - p_clip)),
    weight   = ifelse(is.finite(1 / se_logit^2), 1 / se_logit^2, 1),
    full_wing = factor(full_wing,
                       levels = c("not_fully_winged", "fully_winged")),
    method = factor(method,
                    levels = c("Pitfall2","Pitfall4","Pitfall6",
                               "Pitfall8","Pitfall10","GPD"))
  )

library(lme4)

m_int <- lmer(
  logit_p ~ log_BL * full_wing + (1 | species),
  data    = dat_int,
  weights = weight,
  REML    = TRUE
)

summary(m_int)
library(car)
Anova(m_int, type = 3)
#log_BL:full_wingfully_winged: t val = 0.477
# non significatif → effets additifs (taille + ailes indépendants)
m_no_int <- lmer(
  logit_p ~ log_BL + full_wing + (1 | species),
  data = dat_int, weights = weight
)
anova(m_no_int, m_int)
#pas de différnece entre modele avec ou sans interaction

emm <- emmeans(
  m_int,
  ~ log_BL | full_wing,
  at = list(log_BL = seq(min(dat_int$log_BL),
                         max(dat_int$log_BL), length.out = 100)),
  type = "response"
)

emm_df2 <- as.data.frame(emm) %>%
  mutate(
    prob  = plogis(emmean),
    lwr   = plogis(lower.CL),
    upr   = plogis(upper.CL)
  )

p_traits_inter <- ggplot(emm_df2,
       aes(x = log_BL, y = prob, colour = full_wing, fill = full_wing)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "log10(body length, mm)",
       y = "Detectability p̂",
       colour = "Wing development",
       fill   = "Wing development") +
  theme_minimal()

ggsave(
  file.path(out_dir, "trait/detectability_traits_interaction.png"),
  p_traits_inter,
  width = 12,
  height = 7,
  dpi = 200
)
