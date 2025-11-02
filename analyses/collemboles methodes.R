

# --- Chargement des packages ---
library(tidyverse)

# --- Import des données ---
df <- read.csv("data/raw-data/1.faune/coll25.csv", sep=";")

# Aperçu rapide
head(df)

# === PARTIE 1 — Abondance totale et stabilité ===
# --- Étape 1 : calcul de la somme d’abondance par station et par répétition ---
df_sum <- df %>%
  group_by(station, rep) %>%
  summarise(ABONDANCE_TOTALE = sum(ABONDANCE_TOTALE, na.rm = TRUE), .groups = "drop")

# --- Étape 2 : identifier les réplicats CM1 à CM3 et CM1 à CM6 ---
# stabilité = coefficient de variation (CV = sd / mean)
stability <- df_sum %>%
  group_by(station) %>%
  summarise(
    cv_3 = sd(ABONDANCE_TOTALE[rep %in% c("CM1", "CM2", "CM3")], na.rm = TRUE) /
      mean(ABONDANCE_TOTALE[rep %in% c("CM1", "CM2", "CM3")], na.rm = TRUE),
    cv_6 = sd(ABONDANCE_TOTALE[rep %in% paste0("CM", 1:6)], na.rm = TRUE) /
      mean(ABONDANCE_TOTALE[rep %in% paste0("CM", 1:6)], na.rm = TRUE)
  ) %>%
  mutate(
    plus_stable = if_else(cv_3 < cv_6, "3_reps_plus_stable", "6_reps_plus_stable"),
    gain_stabilite = cv_6 - cv_3
  )

# --- Étape 3 : visualisation de la comparaison ---
ggplot(stability, aes(x = cv_3, y = cv_6)) +
  geom_abline(linetype = 2, color = "grey50") +
  geom_point(aes(color = plus_stable), size = 3) +
  theme_minimal(base_size = 13) +
  labs(
    x = "Coefficient de variation (3 réplicats CM1–CM3)",
    y = "Coefficient de variation (6 réplicats CM1–CM6)",
    color = "Signal le plus stable",
    title = "Comparaison de la stabilité de l'abondance moyenne\n par carotte par station (CV)"
  )

# --- Étape 4 : résumé global ---
stability %>%
  count(plus_stable) %>%
  mutate(pct = round(100 * n / sum(n), 1))

# --- Test global (Wilcoxon) ---
wilcox.test(stability$cv_3, stability$cv_6, paired = TRUE)


# === PARTIE 2 — Composition spécifique : PCoA Bray–Curtis ===
# -- Agrégation par station × rép × taxon --
df_wide <- df %>%
  group_by(station, rep, NOM_CITE) %>%
  summarise(abond = sum(ABONDANCE_TOTALE, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = NOM_CITE, values_from = abond, values_fill = 0)

# --- Moyenne par station sur 3 réplicats (CM1–CM3) ---
mat_3 <- df_wide %>%
  filter(rep %in% c("CM1", "CM2", "CM3")) %>%
  group_by(station) %>%
  summarise(across(-rep, mean, na.rm = TRUE)) %>%
  column_to_rownames("station")

# --- Moyenne par station sur 6 réplicats (CM1–CM6) ---
mat_6 <- df_wide %>%
  filter(rep %in% paste0("CM", 1:6)) %>%
  group_by(station) %>%
  summarise(across(-rep, mean, na.rm = TRUE)) %>%
  column_to_rownames("station")

# --- Vérification : mêmes stations et mêmes espèces ---
mat_3 <- mat_3[, intersect(colnames(mat_3), colnames(mat_6))]
mat_6 <- mat_6[, intersect(colnames(mat_3), colnames(mat_6))]

# --- Distances Bray-Curtis ---
dist_3 <- vegdist(mat_3, method = "bray")
dist_6 <- vegdist(mat_6, method = "bray")

# --- PCoA ---
pcoa_3 <- cmdscale(dist_3, k = 2, eig = TRUE)
pcoa_6 <- cmdscale(dist_6, k = 2, eig = TRUE)

# --- Assemblage pour ggplot ---
df_pcoa <- tibble(
  station = rownames(mat_3),
  x3 = pcoa_3$points[,1],
  y3 = pcoa_3$points[,2],
  x6 = pcoa_6$points[,1],
  y6 = pcoa_6$points[,2]
)

# --- Visualisation comparative ---
ggplot(df_pcoa) +
  geom_segment(aes(x = x3, y = y3, xend = x6, yend = y6),
               arrow = arrow(length = unit(0.1, "cm")), color = "grey60") +
  geom_point(aes(x = x3, y = y3), color = "steelblue", size = 3) +
  geom_point(aes(x = x6, y = y6), color = "tomato", size = 3) +
  theme_minimal(base_size = 13) +
  labs(
    title = "PCoA Bray–Curtis : composition des communautés \n(3 vs 6 réplicats)",
    x = "Axe 1", y = "Axe 2"
  )

# --- Concordance entre ordinations : Procrustes ---
pro <- protest(pcoa_3$points, pcoa_6$points, permutations = 999)
summary(pro)

# --- Test PERMANOVA sur matrice combinée (optionnel) ---
# Combine data long avec facteur "Nb_reps"
mat_all <- bind_rows(
  mat_3 %>% as.data.frame() %>% mutate(nb_rep = "3"),
  mat_6 %>% as.data.frame() %>% mutate(nb_rep = "6")
)
adonis2(mat_all[, -ncol(mat_all)] ~ nb_rep, data = mat_all, method = "bray")
