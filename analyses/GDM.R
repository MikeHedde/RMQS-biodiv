librarian::shelf(ggplot2, dplyr, tidyr, gdm, forcats, ggpubr)



#Données environnementales
habitat <- read.csv("data/derived-data/liste_habitat.csv", h = T, sep = ";") %>%
  rename_with(toupper)
  
sol <- read.csv("data/raw-data/pc_sols.csv", h = T, sep = ";") %>%
  select(STATION, PARAM, valeur) %>%
  pivot_wider(id_cols = STATION, names_from = PARAM, values_from = valeur)
clim <- read.csv("data/derived-data/clim_moy.csv", h = T, sep = ",") 
sampling <- read.csv("data/raw-data/sampling_date.csv", h = T, sep = ";") %>%
  rename_with(toupper)

env <- sol %>%
  left_join(habitat) %>%
  left_join(clim) %>%
  left_join(sampling) %>%
  select(-c(X, DATE, HABITAT_DECRIT, HABITAT_CODE, HABITAT_TYPE)) %>%
  filter(!is.na(Al))

#Données collemboles
data_coll <- read.csv("data/raw-data/collembola.csv", h = T, sep = ";", fileEncoding="latin1")%>%
  filter(RANK == "S") %>%
  select(STATION, LB_NOM, ABONDANCE_TOTALE)%>%
  group_by(STATION, LB_NOM)%>%
  summarise(ab = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(id_cols = STATION, names_from = LB_NOM, values_from = ab, values_fill = 0) %>%
  filter(STATION %in% unique(env$STATION))

#Données isopodes
data_iso <- read.csv("data/raw-data/isopoda.csv", h = T, sep = ";", fileEncoding="latin1")%>%
  filter(RANK == "S") %>%
  select(STATION, LB_NOM, ABONDANCE_TOTALE)%>%
  group_by(STATION, LB_NOM)%>%
  summarise(ab = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(id_cols = STATION, names_from = LB_NOM, values_from = ab, values_fill = 0) %>%
  filter(STATION %in% unique(env$STATION))

#Données diplopodes
data_diplo <- read.csv("data/raw-data/diplopoda.csv", h = T, sep = ";", fileEncoding="latin1")%>%
  filter(RANK == "S") %>%
  select(STATION, LB_NOM, ABONDANCE_TOTALE)%>%
  group_by(STATION, LB_NOM)%>%
  summarise(ab = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(id_cols = STATION, names_from = LB_NOM, values_from = ab, values_fill = 0) %>%
  filter(STATION %in% unique(env$STATION))

# Mise au format pour le GDM
env2 <- env[env$STATION %in% unique(data_coll$STATION),] %>%
  select(STATION, 
         HAB,                                       #habitat
                                                    #soil granulo
         CN, Ctot, COT,                             #C resource
         N, Ca, Mg, Mn, K, Na, calc, CEC,           #nutrients
         DD_cum, T2MDEW_moy, T2max_moy, EVLAND_moy, #climat
         pH, Al, Hum_res,                           #stress
         LONGITUDE, LATITUDE,                       #geographic
         ALTITUDE) %>%                              #altitude
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

# GDM 
gdm.1 <- gdm(data=gdmdata_coll, geo=TRUE)
gdm.2 <- gdm(data=gdmdata_iso, geo=TRUE)
gdm.3 <- gdm(data=gdmdata_diplo, geo=TRUE)
summary(gdm.1)

modTest <- gdm.varImp(gdmdata_coll, geo=T, nPerm=50, parallel=T, cores=10, predSelect=T)
barplot(sort(modTest$`Predictor Importance`[,1], decreasing=T))

# extraction des coefficients de chaque GDM
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
  mutate(grp = case_when(predictors == "HAB" ~ "habitat",
                         predictors %in% c("CN", "Ctot", "COT", "N", "Ca", "Mg", "Mn", "K", "Na", "calc", "CEC") ~ "nutrients & C",
                         predictors %in% c("DD_cum", "T2MDEW_moy", "T2max_moy", "EVLAND_moy") ~ "climat",
                         predictors %in% c("pH", "Al", "Hum_res") ~ "stress",
                         predictors == "ALTITUDE" ~ "altitude", 
                         predictors == "Geographic" ~ "geographic")) %>%
  mutate(grp = fct_relevel(grp, c("geographic", "habitat", "climat", "altitude", 
                          "stress", "nutrients & C")))

sumCoef_p <- ggplot(sumCoef, aes(x=taxa, y = coef, fill = grp))+
  geom_bar(position="stack", stat="identity", colour = "grey")+
  coord_flip()+
  labs(x = "", y = "Percentage of deviance explained", 
       fill = "Group of \npredictors") +
  theme_bw()

# Effet de la distance écologique
predecol_coll <- data.frame(taxa = "coll", obs = gdm.1$observed, ecol = gdm.1$ecological, pred = gdm.1$predicted)
predecol_iso <- data.frame(taxa = "iso", obs = gdm.2$observed, ecol = gdm.2$ecological, pred = gdm.2$predicted)
predecol_diplo <- data.frame(taxa = "diplo", obs = gdm.3$observed, ecol = gdm.3$ecological, pred = gdm.3$predicted)

predecol <- tibble(rbind(predecol_coll, predecol_iso, predecol_diplo))

predecol_p <- ggplot(predecol, aes(x = ecol, y = obs))+
  geom_point(aes(colour = taxa, alpha = 0.01))+
  geom_line(aes(x = ecol, y = pred, colour = taxa), lwd=2)+
  facet_wrap(~taxa, scales = "free_x")+
  labs(x = "Predicted ecological distance", y = "Composition dissimilarity")+
  theme_bw()+
  theme(legend.position="none")

# Splines
sumCoef_wide <- pivot_wider(sumCoef, id_cols = c(predictors, grp), names_from = taxa, values_from = coef)

gdm.1.y <- as.data.frame(isplineExtract(gdm.1)$y) %>%
  mutate(taxa = "collembola", param = "ecol_dist")
gdm.1.x <- as.data.frame(isplineExtract(gdm.1)$x)  %>%
  mutate(taxa = "collembola", param = "predictor")
gdm.2.y <- as.data.frame(isplineExtract(gdm.2)$y)  %>%
  mutate(taxa = "isopoda", param = "ecol_dist")
gdm.2.x <- as.data.frame(isplineExtract(gdm.2)$x)  %>%
  mutate(taxa = "isopoda", param = "predictor")
gdm.3.y <- as.data.frame(isplineExtract(gdm.3)$y)  %>%
  mutate(taxa = "diplopoda", param = "ecol_dist")
gdm.3.x <- as.data.frame(isplineExtract(gdm.3)$x)  %>%
  mutate(taxa = "diplopoda", param = "predictor")

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

splines_gdm <- rbind(splines_COT, splines_geo)

splines_p <- ggplot(splines_gdm, aes(x = x, y = y, colour = taxa))+
  geom_line()+
  labs(y="Partial ecological distance", x="")+
  facet_wrap(.~param, scales = "free_x",
             labeller = as_labeller(c(COT = "Tot Organic C\n(mg. g-1)", geo = "Geographic")),
             strip.position = "bottom")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside")


ggarrange(sumCoef_p, predecol_p, splines_p, 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
