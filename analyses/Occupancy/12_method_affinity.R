library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(readr)

# --- paramètres ---
top_k <- 12          # nb d'espèces par méthode
thr   <- 0.10        # seuil d'affinité utile (Δ >= 0.10)
cap_p <- 0.98        # option: retirer quasi-plafond pour éviter les "faux tops"
methods_order <- c("GPD","DVAC")

# --- 1) centrage par méthode : Δ = p_hat - médiane_méthode ---
med_by_method <- tab_p %>% group_by(method) %>%
  summarise(p_med = median(p_hat, na.rm = TRUE), .groups = "drop")

tab_delta <- tab_p %>%
  inner_join(med_by_method, by = "method") %>%
  mutate(delta = p_hat - p_med,
         genus = str_extract(species, "^[^\\s]+"),
         method = factor(method, levels = methods_order))

# (option) éviter la saturation : on exclut p_hat > cap_p pour Pitfall8/10
tab_delta_clean <- tab_delta %>%
  filter(!(method %in% c("Pitfall8","Pitfall10") & p_hat > cap_p))

# --- 2) sélection des espèces affines (Δ >= thr) et top-K par méthode ---
top_affines <- tab_delta_clean %>%
  filter(delta >= thr) %>%
  group_by(method) %>%
  arrange(desc(delta), .by_group = TRUE) %>%
  slice_head(n = top_k) %>%
  ungroup() %>%
  mutate(label = paste0(species, " (", sprintf("+%.2f", delta), ")"))

# Si trop peu d'espèces dépassent le seuil, on complète jusqu’à top_k
fill_up <- tab_delta_clean %>%
  group_by(method) %>%
  arrange(desc(delta), .by_group = TRUE) %>%
  slice_head(n = top_k) %>%
  ungroup() %>%
  anti_join(top_affines, by = c("species","method")) %>%
  mutate(label = paste0(species, " (", sprintf("+%.2f", delta), ")"))

top_affines2 <- bind_rows(top_affines, fill_up) %>%
  group_by(method) %>%
  slice_head(n = top_k) %>%
  ungroup()

# ordre des espèces par méthode (du plus affine au moins)
top_affines2 <- top_affines2 %>%
  group_by(method) %>%
  mutate(species_ord = fct_reorder(species, delta, .desc = TRUE)) %>%
  ungroup()%>%
  filter(method %in% methods_order)

# --- 3) graphique : lollipop facetté par méthode ---
p_aff <- ggplot(top_affines2,
                aes(x = species_ord, y = delta)) +
  geom_segment(aes(xend = species_ord, y = 0, yend = delta),
               linewidth = 0.6, color = "grey60") +
  geom_point(size = 2.4, alpha = .9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  coord_flip() +
  facet_wrap(~ method, scales = "free_y", ncol = 3) +
  scale_y_continuous("Δ = p̂ - médiane(méthode)", limits = c(0, NA)) +
  labs(x = NULL,
       title = "Espèces les plus « affines » par méthode") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(face = "bold"))

# export
dir.create(file.path(out_dir, "affinite_species"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "affinite_species/top_affines_lollipop.png"),
       p_aff, width = 12, height = 8, dpi = 300)

# --- 4) (bonus) heatmap espèce × méthode pour un top global ---
top_global <- tab_delta_clean %>%
  arrange(desc(delta)) %>% slice_head(n = 50) %>%
  pull(species) %>% unique()

heat_df <- tab_delta_clean %>%
  filter(species %in% top_global) %>%
  mutate(species = fct_reorder(species, delta, .fun = max, .desc = TRUE),
         method  = factor(method, levels = methods_order))

p_heat <- ggplot(heat_df, aes(method, species, fill = delta)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Δ", option = "C") +
  labs(x = "Méthode", y = "Espèce",
       title = "Affinité relative (Δ) par méthode — top espèces (global)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(out_dir, "affinite_species/affinity_heatmap.png"),
       p_heat, width = 8, height = 12, dpi = 300)
