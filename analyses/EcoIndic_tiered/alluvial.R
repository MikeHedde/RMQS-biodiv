
library(tidyverse)
library(ggalluvial)
library(ggnewscale)

#========================
# 1. Import
#========================
nem <- read.csv("data/raw-data/1.faune/nematoda.csv", h = T, sep = ";")
euka <- read.csv("data/raw-data/1.faune/rmqs_2024_euka02_nematoda_reads.csv", h = T, sep = ";")
resp <- read.csv("data/raw-data/3.fonctions/2024_resp.csv", h = T, sep = ";")

# Vérification des noms de colonnes
names(euka)
names(nem)
names(resp)

# Harmonisation du nom de la clé
nem  <- nem  %>% rename(station = STATION)

# Vérification des stations communes
stations_communes <- reduce(
  list(unique(euka$station), unique(nem$station), unique(resp$station)),
  intersect
)


#========================
# 2. Fonction tertiles
#========================
# ntile() répartit les stations en 3 groupes d'effectifs proches
# même en présence de valeurs ex aequo.
make_tertile <- function(x) {
  dplyr::ntile(x, 3)
}

tertile_lab <- c(
  `1` = "T1",
  `2` = "T2",
  `3` = "T3"
)

#========================================
# 3. Fonction pour attribuer les tertiles
#========================================
make_tertile <- function(x) {
  dplyr::ntile(x, 3)
}

tertile_lab <- c(`1` = "T1", `2` = "T2", `3` = "T3")

#========================================
# 4. Richesse en MOTUs par station
#    = nombre de colonnes > 0 par ligne
#========================================
motus_df <- euka %>%
  mutate(across(-station, ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(
    motus_richness = rowSums(across(-station, ~ .x > 0), na.rm = TRUE),
    motus_tertile = factor(
      tertile_lab[as.character(make_tertile(motus_richness))],
      levels = c("T1", "T2", "T3")
    )
  ) %>%
  select(station, motus_richness, motus_tertile)

#========================================
# 5. Extraction de SI dans nematoda.csv
#========================================
si_df <- nem %>%
  filter(indice == "SI") %>%
  mutate(
    SI = suppressWarnings(as.numeric(valeur)),
    si_tertile = factor(
      tertile_lab[as.character(make_tertile(SI))],
      levels = c("T1", "T2", "T3")
    )
  ) %>%
  select(station, SI, si_tertile)

#========================================
# 6. Extraction de tot_resp dans 2024_resp.csv
#========================================
resp_df <- resp %>%
  filter(measure == "tot_resp") %>%
  mutate(
    tot_resp = suppressWarnings(as.numeric(mean)),
    resp_tertile = factor(
      tertile_lab[as.character(make_tertile(tot_resp))],
      levels = c("T1", "T2", "T3")
    )
  ) %>%
  select(station, tot_resp, resp_tertile)

#========================================
# 7. Fusion des trois objets
#========================================
alluv_df <- motus_df %>%
  inner_join(si_df, by = "station") %>%
  inner_join(resp_df, by = "station") %>%
  mutate(
    changed = !(motus_tertile == si_tertile & si_tertile == resp_tertile),
    profile = paste(motus_tertile, si_tertile, resp_tertile, sep = "_")
  ) %>%
  arrange(motus_tertile, si_tertile, resp_tertile, station)

# Vérification
cat("Nombre de stations communes :", nrow(alluv_df), "\n")
print(head(alluv_df))

#========================================
# 8. Palette sobre
#========================================
base_tertile_cols <- c(
  "T1" = "#5F8D7A",  # vert sauge
  "T2" = "#B9854B",  # ocre doux
  "T3" = "#6C5B7B"   # prune grisée
)

#========================================
# 9. Fonction de mélange RGB
#    Ex : T1_T1_T2 = 2/3 T1 + 1/3 T2
#========================================
mix_cols <- function(cols) {
  rgb_mat <- grDevices::col2rgb(cols)
  rgb(
    red   = mean(rgb_mat[1, ]) / 255,
    green = mean(rgb_mat[2, ]) / 255,
    blue  = mean(rgb_mat[3, ]) / 255
  )
}

profiles <- unique(alluv_df$profile)

flow_cols <- setNames(
  vapply(profiles, function(p) {
    terts <- strsplit(p, "_")[[1]]
    mix_cols(base_tertile_cols[terts])
  }, character(1)),
  profiles
)

# Palette globale : colonnes + flux
all_cols <- c(base_tertile_cols, flow_cols)

#========================================
# 10. Mise en forme longue pour ggalluvial
#========================================
plot_df <- bind_rows(
  alluv_df %>%
    transmute(
      station,
      axis = "MOTUs",
      stratum = motus_tertile,
      profile,
      changed
    ),
  alluv_df %>%
    transmute(
      station,
      axis = "SI",
      stratum = si_tertile,
      profile,
      changed
    ),
  alluv_df %>%
    transmute(
      station,
      axis = "Resp",
      stratum = resp_tertile,
      profile,
      changed
    )
) %>%
  mutate(
    axis = factor(axis, levels = c("MOTUs", "SI", "Resp")),
    stratum = factor(stratum, levels = c("T1", "T2", "T3")),
    fill_key = profile,
    alpha_flow = ifelse(changed, 0.90, 0.28)
  )

#========================================
# 11. Alluvial plot
#========================================
p <- ggplot(
  plot_df,
  aes(x = axis, stratum = stratum, alluvium = station, y = 1)
) +
  geom_flow(
    aes(fill = fill_key, alpha = alpha_flow),
    stat = "alluvium",
    width = 0.18,
    knot.pos = 0.45,
    color = scales::alpha("grey20", 0.35),
    linewidth = 0.15,
    show.legend = FALSE
  ) +
  geom_stratum(
    aes(fill = stratum),
    width = 0.18,
    color = "white",
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = all_cols) +
  scale_alpha_identity() +
  labs(
    x = NULL,
    y = NULL #,
    #title = "Transitions de tertiles entre MOTUs, SI et Resp",
    #subtitle = "Couleur des flux = mélange des tertiles traversés ; couleur des colonnes = tertile local"
  ) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 28, face = "bold"),
    plot.subtitle = element_text(size = 14),
    legend.position = "none"
  )

print(p)

#========================================
# 12. Export
#========================================
ggsave(
  filename = "figures/EcoInd_tiered/alluvial_final.png",
  plot = p,
  width = 11,
  height = 7,
  dpi = 300
)


##diagramme ternaire

library(tidyverse)
library(ggtern)

#---------------------------------------
# 1. On part de l'objet déjà construit
#---------------------------------------
# alluv_df contenant :
# station, motus_richness, SI, tot_resp

df <- alluv_df %>%
  select(station, motus_richness, SI, tot_resp)

#---------------------------------------
# 2. Calcul des z-scores
#---------------------------------------
df_z <- df %>%
  mutate(
    z_motus = scale(motus_richness)[,1],
    z_SI    = scale(SI)[,1],
    z_resp  = scale(tot_resp)[,1]
  )

#---------------------------------------
# 3. Rendre les valeurs positives
#---------------------------------------
# On décale par la valeur minimale

shift <- abs(min(df_z %>% select(z_motus, z_SI, z_resp))) + 0.01

df_z <- df_z %>%
  mutate(
    motus_pos = z_motus + shift,
    SI_pos    = z_SI + shift,
    resp_pos  = z_resp + shift
  )

#---------------------------------------
# 4. Normalisation pour coordonnées ternaires
#---------------------------------------

df_tern <- df_z %>%
  rowwise() %>%
  mutate(
    total = motus_pos + SI_pos + resp_pos,
    T = motus_pos / total,
    L = SI_pos / total,
    R = resp_pos / total
  ) %>%
  ungroup() %>%
  left_join(
    alluv_df %>% select(station, profile),
    by = "station"
  )
#---------------------------------------
# 5. Diagramme ternaire
#---------------------------------------

p <- ggtern(
  df_tern,
  aes(
    x = T,
    y = L,
    z = R,
    color = profile
  )
) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  scale_color_manual(values = flow_cols) +
  facet_wrap(profile~.)+
  theme_bw() +
  labs(
    T = "MOTUs richness",
    L = "Structure Index",
    R = "Respiration",
    title = "Position des stations dans l'espace des indicateurs",
    subtitle = "Couleur = trajectoire de tertiles (MOTUs → SI → Resp)"
  ) +
  theme(
    tern.axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold")
  )

print(p)


# Corrélations entre indicateurs

df <- alluv_df %>%
  select(station, motus_richness, SI, tot_resp)

cor_motus_si   <- cor.test(df$motus_richness, df$SI, method = "spearman")
cor_motus_resp <- cor.test(df$motus_richness, df$tot_resp, method = "spearman")
cor_si_resp    <- cor.test(df$SI, df$tot_resp, method = "spearman")

cor_motus_si
cor_motus_resp
cor_si_resp

#Accord de classement entre tertiles

tab_motus_si <- table(alluv_df$motus_tertile, alluv_df$si_tertile)
tab_motus_resp <- table(alluv_df$motus_tertile, alluv_df$resp_tertile)
tab_si_resp <- table(alluv_df$si_tertile, alluv_df$resp_tertile)

tab_motus_si
tab_motus_resp
tab_si_resp

#test d’indépendance simple

chisq.test(tab_motus_si)
chisq.test(tab_motus_resp)
chisq.test(tab_si_resp)

# scatterplots avec rho affiché


ggplot(df, aes(motus_richness, SI)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal()

ggplot(df, aes(motus_richness, tot_resp)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal()

ggplot(df, aes(SI, tot_resp)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal()


# Modèle nul simple
# On mélange les indicateurs aléatoirement entre stations et on regarde la fréquence des spans.
set.seed(1)

n_sim <- 10000

span_null <- replicate(n_sim, {
  
  sim_df <- alluv_df %>%
    mutate(
      si_tertile = sample(si_tertile),
      resp_tertile = sample(resp_tertile)
    )
  
  sim_span <- sim_df %>%
    mutate(
      motus_num = as.numeric(motus_tertile),
      si_num = as.numeric(si_tertile),
      resp_num = as.numeric(resp_tertile),
      span = pmax(motus_num, si_num, resp_num) -
        pmin(motus_num, si_num, resp_num)
    )
  
  mean(sim_span$span == 2)
  
})

mean(span_null)
sd(span_null)

obs <- mean(span_df$span == 2)
p_value <- mean(span_null >= obs)

# Quantifier l’accord entre indicateurs
library(irr)

kappa2(alluv_df[,c("motus_tertile","si_tertile")])
kappa2(alluv_df[,c("motus_tertile","resp_tertile")])
kappa2(alluv_df[,c("si_tertile","resp_tertile")])
