# ============================================================
# TRAITS → DÉTECTABILITÉ (p̂) ~ méthode * log(masse)
# Sorties: modèles, figures
# ============================================================

tab_p <- readr::read_csv(file.path(out_dir, "p_hat_by_method_unmarked_full.csv"), show_col_types = FALSE)

library(taxize)
trait0 <- read.csv("data/raw-data/1.faune/trait_carabidae.csv", h=T, sep=";")

# Max Body length
trait <- trait0 %>%
  mutate(sp = taxize::gna_parse(species)$canonical_full) %>%
  filter(raw_trait_value == "Body_length") %>%
  mutate(trait_value = as.numeric(trait_value))  %>%
  group_by(sp) %>%
  summarise(max_BL = max(trait_value))%>%
  rename(species = sp)

#Si trait mesuré dans la bdd (à éviter)
trait <- raw %>%
  filter(PROJET == target_project) %>%
  select(LB_NOM, `MASSE (mg)`) %>%
  rename(species = LB_NOM) %>%
  group_by(species) %>%
  summarise(mass = median(`MASSE (mg)`, na.rm = TRUE), .groups="drop")

dat <- tab_p %>% left_join(trait, by = "species")

# Mise à l’échelle logit + poids méta depuis IC95%
dat2 <- dat %>%
  mutate(
    p_clip  = pmin(pmax(p_hat, 1e-5), 1 - 1e-5),
    logit_p = logit(p_clip),
    se_p    = ifelse(!is.na(lcl) & !is.na(ucl), (ucl - lcl) / (2 * 1.96), NA_real_),
    se_logit = ifelse(!is.na(se_p), se_p / (p_clip * (1 - p_clip)), NA_real_),
    weight   = ifelse(is.finite(1/(se_logit^2)), 1/(se_logit^2), 1),
    log_BL = log10(max_BL)
  ) %>%
  filter(!is.na(log_BL), is.finite(log_BL)) %>%
  mutate(species = factor(species))

m_tr_int <- lme4::lmer(logit_p ~ method * log_BL + (1 | species),
                         data = dat2, weights = weight, REML = TRUE)
s <- summary(m_tr_int); capture.output(s, file = file.path(out_dir, "model_traits_detectability_summary.txt"))

trends <- emmeans::emtrends(m_tr_int, ~ method, var = "log_BL")
capture.output(trends, file = file.path(out_dir, "model_traits_emtrends.txt"))
capture.output(pairs(trends), file = file.path(out_dir, "model_traits_emtrends_pairs.txt"))

# Prédictions lissées sur gamme de masses par méthode
newgrid <- dat2 %>%
  group_by(method) %>%
  summarise(m_min = quantile(log_BL, 0.05), m_max = quantile(log_BL, 0.95), .groups="drop") %>%
  rowwise() %>% do(tibble(method = .$method, log_BL = seq(.$m_min, .$m_max, length.out = 100))) %>% ungroup()

pred <- cbind(newgrid,
              predict(m_tr_int, newdata = newgrid, re.form = NA, se.fit = TRUE)) %>%
  as_tibble() %>%
  mutate(p_fit = plogis(fit), p_lcl = plogis(fit - 1.96*se.fit), p_ucl = plogis(fit + 1.96*se.fit))

ggplot(dat2, aes(x = log_BL, y = p_clip, colour = method)) +
  geom_point(alpha = 0.5) +
  geom_line(data = pred, aes(y = p_fit), linewidth = 1) +
  geom_ribbon(data = pred, aes(x = log_BL,ymin = p_lcl, ymax = p_ucl, fill = method), alpha = 0.15, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(method~.) +
  labs(x = "log10(body length, mm)", y = "Detectability p̂ (back-transformed)",
       title = "Effect of body length on detectability by method") +
  theme_minimal() +
  theme(legend.position = "none") ->
  p_BL

ggsave(file.path(out_dir, "detectability_vs_BL_by_method.png"), p_BL, width = 9, height = 6, dpi = 200)

# Motion strategy
trait <- trait0 %>%
  mutate(sp = taxize::gna_parse(species)$canonical_full) %>%
  filter(raw_trait_value == "Wing_development") %>%
  mutate(trait_value = case_when(
    str_detect(trait_value, regex("macropterous", ignore_case = TRUE)) ~ "fully_winged",
    str_detect(trait_value, regex("brachypterous|apterous", ignore_case = TRUE)) ~ "not_winged",
    TRUE ~ NA_character_
  ))%>%
  group_by(sp) %>%
  summarise(
    n_records   = sum(!is.na(trait_value)),
    n_winged    = sum(trait_value == "fully_winged", na.rm = TRUE),
    n_notwinged = sum(trait_value == "not_winged",  na.rm = TRUE),
    pct_winged  = n_winged / n_records,
    pct_notwinged = n_notwinged / n_records
  ) %>%
  ungroup() %>%
  rename(species = sp) %>%
  select(species, pct_winged) %>%
  mutate(full_wing = case_when(
    pct_winged == 1 ~ "fully_winged",
    pct_winged < 1 ~ "not_fully_winged"
  ))

dat <- tab_p %>% left_join(trait, by = "species")


ggplot(dat, aes(x = full_wing, y = p_hat, colour = method)) +
  geom_boxplot()+
  geom_point(alpha = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(method~.) +
  labs(x = "Wing development (% of individuals)", y = "Detectability p̂ (back-transformed)",
       title = "Effect of wing development on detectability by method") +
  theme_minimal() +
  theme(legend.position = "none") ->
  p_wing

ggsave(file.path(out_dir, "detectability_vs_wing_by_method.png"), p_wing, width = 12, height = 6, dpi = 200)
