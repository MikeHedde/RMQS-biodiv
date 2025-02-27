librarian::shelf(ggplot2, dplyr, tidyr, gdm)



#Données environnementales
habitat <- read.csv("data/derived-data/liste_habitat.csv", h = T, sep = ";") %>%
  rename_with(toupper) %>% 
  mutate(
    HAB = case_when(
      HABITAT_TYPE=="FERME" ~ 0,
      HABITAT_TYPE=="OUVERT" ~ 1
    ))
  
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

#Données bio
data_coll <- read.csv("data/raw-data/collembola.csv", h = T, sep = ";", fileEncoding="latin1")%>%
  filter(RANK == "S") %>%
  select(STATION, LB_NOM, ABONDANCE_TOTALE)%>%
  group_by(STATION, LB_NOM)%>%
  summarise(ab = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(id_cols = STATION, names_from = LB_NOM, values_from = ab, values_fill = 0) %>%
  filter(STATION %in% unique(env$STATION))


# Mise au format pour le GDM
env2 <- env[env$STATION %in% unique(data_coll$STATION),] %>%
  select(STATION, N, CN, Ca, Mg, Mn, K, Na, Ctot, HAB, calc, COT, CEC, pH, DD_cum, 
         T2MDEW_moy, T2max_moy, EVLAND_moy, LONGITUDE, LATITUDE, ALTITUDE) %>%
  mutate_if(is.character, as.numeric)
env2 <- as.data.frame(env2)
gdmdata <- formatsitepair(data_coll, 
                         1, 
                         siteColumn="STATION", 
                         abundance = T,
                         XColumn="LONGITUDE", 
                         YColumn="LATITUDE",
                         predData=env2) %>%
  drop_na()

# GDM 
gdm.1 <- gdm(data=gdmdata, geo=TRUE)
summary(gdm.1)
gdm.crossval <- gdm.crossvalidation(gdmdata,train.proportion=0.9, n.crossvalid.tests=1,
                    geo=T, splines=NULL, knots=NULL)
plot(gdm.1, plot.layout=c(3,6))
gdm.1.splineDat <- isplineExtract(gdm.1)

