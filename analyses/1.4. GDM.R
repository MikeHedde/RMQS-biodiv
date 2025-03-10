# Analyse de la bêta-diversité des détritivores avec GDM
# Auteur : Mickael Hedde
# Date : 09-03-2025
# Objectif : Modéliser la dissimilarité de composition des communautés de détritivores en fonction des variables environnementales

# Chargement des bibliothèques nécessaires
librarian::shelf(ggplot2, dplyr, tidyr, gdm, forcats, ggpubr, paletteer)

# Chargement et préparation des données environnementales
habitat <- read.csv("data/derived-data/liste_habitat.csv", h = T, sep = ";") %>%
  rename_with(toupper)
sol <- read.csv("data/derived-data/pc_sols.csv", h = T, sep = ",")
clim <- read.csv("data/derived-data/clim_stat.csv", h = T, sep = ",") %>%
  select(-X)
sampling <- read.csv("data/raw-data/2.env/sampling_date.csv", h = T, sep = ";") %>%
  rename_with(toupper)
ndvi <- read.csv("data/derived-data/ndvi.csv", h=T, sep=",") %>%
  select(site, ndvi_m, ndvi_max, ndvi_sd) %>%
  rename(STATION = site)

# Fusion des données environnementales
env <- ndvi %>%
  left_join(habitat) %>%
  left_join(clim) %>%
  left_join(sampling) %>%
  left_join(sol) %>%
  select(-c(X, DATE, HABITAT_DECRIT, HABITAT_CODE, HABITAT_TYPE))

# Sauvegarde des données environnementales consolidées
write.csv(env, "data/derived-data/all_env_variables.csv")

# Chargement des données des différents taxons
## Collemboles
data_coll <- read.csv("data/raw-data/1.faune/collembola.csv", h = T, sep = ";", fileEncoding="latin1")%>%
  filter(RANK == "S") %>%
  select(STATION, LB_NOM, ABONDANCE_TOTALE)%>%
  group_by(STATION, LB_NOM)%>%
  summarise(ab = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(id_cols = STATION, names_from = LB_NOM, values_from = ab, values_fill = 0) %>%
  filter(STATION %in% unique(env$STATION))

## isopodes
data_iso <- read.csv("data/raw-data/1.faune/isopoda.csv", h = T, sep = ";", fileEncoding="latin1")%>%
  filter(RANK == "S") %>%
  select(STATION, LB_NOM, ABONDANCE_TOTALE)%>%
  group_by(STATION, LB_NOM)%>%
  summarise(ab = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(id_cols = STATION, names_from = LB_NOM, values_from = ab, values_fill = 0) %>%
  filter(STATION %in% unique(env$STATION))

## diplopodes
data_diplo <- read.csv("data/raw-data/1.faune/diplopoda.csv", h = T, sep = ";", fileEncoding="latin1")%>%
  filter(RANK == "S") %>%
  select(STATION, LB_NOM, ABONDANCE_TOTALE)%>%
  group_by(STATION, LB_NOM)%>%
  summarise(ab = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(id_cols = STATION, names_from = LB_NOM, values_from = ab, values_fill = 0) %>%
  filter(STATION %in% unique(env$STATION))

# Formatage des données pour le GDM
env2 <- env[env$STATION %in% unique(data_coll$STATION),] %>%
  select(STATION, 
         HAB, ndvi_m, ndvi_sd,                                               #habitat
         clay, sand,                                                         #soil granulo
         CN, COT, calc, CEC,                                                 #C resource & nutrients
         DD_cum, T360_mean, #T30_mean, T30_max, T30_min, T30_sd, T360_sd,    #climat
         pH, Cd_disp, Pb_disp, #Cu_disp,                                      #stress
         LONGITUDE, LATITUDE                                                 #geographic
         ) %>%                              
  mutate_if(is.character, as.numeric)
env2 <- as.data.frame(env2)

gdmdata_coll <- formatsitepair(data_coll, 
                         1, 
                         siteColumn="STATION", 
                         abundance = F,
                         XColumn="LONGITUDE", 
                         YColumn="LATITUDE",
                         predData=env2) %>%
  drop_na()

gdmdata_iso <- formatsitepair(data_iso, 
                               1, 
                               siteColumn="STATION", 
                               abundance = F,
                               XColumn="LONGITUDE", 
                               YColumn="LATITUDE",
                               predData=env2) %>%
  drop_na()

gdmdata_diplo <- formatsitepair(data_diplo, 
                              1, 
                              siteColumn="STATION", 
                              abundance = F,
                              XColumn="LONGITUDE", 
                              YColumn="LATITUDE",
                              predData=env2) %>%
  drop_na()

# Modélisation GDM 
gdm.1 <- gdm(data=gdmdata_coll, geo=TRUE)
gdm.2 <- gdm(data=gdmdata_iso, geo=TRUE)
gdm.3 <- gdm(data=gdmdata_diplo, geo=TRUE)
summary(gdm.1)

# Extraction des coefficients
sumCoef_coll <- as.data.frame(tapply(gdm.1$coefficients, 
                                (seq_along(gdm.1$coefficients) - 1) %/% 3, 
                                sum))%>%
  mutate(predictors=gdm.1$predictors, taxa = "collembola") %>%
  rename("coef" = 1)

sumCoef_iso <- as.data.frame(tapply(gdm.2$coefficients, 
                                     (seq_along(gdm.2$coefficients) - 1) %/% 3, 
                                     sum))%>%
  mutate(predictors=gdm.2$predictors, taxa = "isopoda") %>%
  rename("coef" = 1)

sumCoef_diplo <- as.data.frame(tapply(gdm.3$coefficients, 
                                    (seq_along(gdm.3$coefficients) - 1) %/% 3, 
                                    sum))%>%
  mutate(predictors=gdm.3$predictors, taxa = "diplopoda") %>%
  rename("coef" = 1)

sumCoef <- rbind(sumCoef_coll, sumCoef_iso, sumCoef_diplo) %>%
  mutate(grp = case_when(predictors %in% c("HAB", "ndvi_m", "ndvi_max", "ndvi_sd") ~ "habitat",
                         predictors %in% c("clay", "sand") ~ "soil texture",
                         predictors %in% c("CN", "Ctot", "COT", "N", "Ca", "Mg", "Mn", "K", "Na", "calc", "CEC") ~ "nutrients & C",
                         predictors %in% c("DD_cum", "T30_mean", "T30_max", "T30_min", "T30_sd", "T360_sd", "T360_mean") ~ "climat",
                         predictors %in% c("pH", "Al", "Hum_res", "Cd_disp", "Pb_disp") ~ "stress",
                         predictors %in% c("Geographic") ~ "geographic")) %>%
  mutate(grp = fct_relevel(grp, c("geographic", "habitat", "climat", "stress", "nutrients & C")))

sumCoef_p <- ggplot(sumCoef, aes(x=taxa, y = coef, fill = grp))+
  geom_bar(position="stack", stat="identity", colour = "grey")+
  ylim(0,20)+
  scale_fill_paletteer_d("MexBrewer::Alacena") +
  coord_flip()+
  labs(x = "", y = "Percentage of deviance explained", 
       fill = "Group of \npredictors") +
  theme_bw()+
  theme(legend.position=c(.85,.5),
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        axis.text.y = element_text(size = 6))+
  guides(fill = guide_legend(override.aes = list(size = 0.5)))

# Effet de la distance écologique
predecol_coll <- data.frame(taxa = "Collembola", obs = gdm.1$observed, ecol = gdm.1$ecological, pred = gdm.1$predicted)
predecol_iso <- data.frame(taxa = "Isopoda", obs = gdm.2$observed, ecol = gdm.2$ecological, pred = gdm.2$predicted)
predecol_diplo <- data.frame(taxa = "Diplopoda", obs = gdm.3$observed, ecol = gdm.3$ecological, pred = gdm.3$predicted)

predecol <- tibble(rbind(predecol_coll, predecol_iso, predecol_diplo))


predecol_p <- ggplot(predecol, aes(x = ecol, y = obs))+
  geom_point(aes(colour = taxa, alpha = 0.1))+
  geom_line(aes(x = ecol, y = pred, colour = taxa), lwd=2)+
  scale_colour_paletteer_d("nationalparkcolors::ArcticGates")+
  facet_wrap(~taxa, scales = "free_x")+
  labs(x = "Predicted ecological distance", y = "Composition dissimilarity")+
  theme_bw()+
  theme(legend.position="none",
        strip.background=element_rect(fill="white"))

# Splines
sumCoef_wide <- pivot_wider(sumCoef, id_cols = c(predictors, grp), names_from = taxa, values_from = coef)

gdm.1.y <- as.data.frame(isplineExtract(gdm.1)$y) %>%
  mutate(taxa = "Collembola", param = "ecol_dist")
gdm.1.x <- as.data.frame(isplineExtract(gdm.1)$x)  %>%
  mutate(taxa = "Collembola", param = "predictor")
gdm.2.y <- as.data.frame(isplineExtract(gdm.2)$y)  %>%
  mutate(taxa = "Isopoda", param = "ecol_dist")
gdm.2.x <- as.data.frame(isplineExtract(gdm.2)$x)  %>%
  mutate(taxa = "Isopoda", param = "predictor")
gdm.3.y <- as.data.frame(isplineExtract(gdm.3)$y)  %>%
  mutate(taxa = "Diplopoda", param = "ecol_dist")
gdm.3.x <- as.data.frame(isplineExtract(gdm.3)$x)  %>%
  mutate(taxa = "Diplopoda", param = "predictor")

param.x <- rbind(gdm.1.x, gdm.2.x, gdm.3.x) %>%
  rename_at(vars(-taxa, -param), ~ paste0(., '_x'))
param.y <- rbind(gdm.1.y, gdm.2.y, gdm.3.y) %>%
  rename_at(vars(-taxa, -param), ~ paste0(., '_y'))

splines_COT <- data.frame(x = param.x$COT_x, 
                     y = param.y$COT_y, 
                     taxa = param.x$taxa,
                     param = "COT")

splines_geo <- data.frame(x = param.x$Geographic_x, 
                          y = param.y$Geographic_y, 
                          taxa = param.x$taxa,
                          param = "geo")

splines_ndvi <- data.frame(x = param.x$ndvi_m_x, 
                          y = param.y$ndvi_m_y, 
                          taxa = param.x$taxa,
                          param = "NDVI")

splines_T360mean <- data.frame(x = param.x$T360_mean_x, 
                           y = param.y$T360_mean_y, 
                           taxa = param.x$taxa,
                           param = "T360_m")

splines_gdm <- rbind(splines_COT, splines_ndvi, splines_geo, splines_T360mean)

splines_p <- ggplot(splines_gdm, aes(x = x, y = y, colour = taxa))+
  geom_line(linewidth = 1.5)+
  scale_colour_paletteer_d("nationalparkcolors::ArcticGates")+
  labs(y="Partial ecological distance", x="")+
  facet_wrap(param~., scales = "free_x",
             labeller = as_labeller(c(COT = "Tot Organic C\n(mg. g-1)", 
                                      NDVI = "NDVI",
                                      geo = "Geographic",
                                      T360_m = "Average daily Temperature (°C) \nover last 360 days")),
             strip.position = "bottom", nrow = 1)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

assembled_GDM <- ggarrange(ggarrange(sumCoef_p, predecol_p,  
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1),
          splines_p, nrow = 2, labels = "C", label.y = 0.1)

ggsave("figures/article ASE/GDM detritivore.png", 
       width = 21, height = 12, units = "cm")
