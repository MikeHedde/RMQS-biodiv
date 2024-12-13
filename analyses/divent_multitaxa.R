librarian::shelf(ggplot2, dplyr)


data_car <- read.csv("data/raw-data/carabidae.csv", h = T, sep = ";", fileEncoding="latin1")
data_ara <- read.csv("data/raw-data/araneae.csv", h = T, sep = ";", fileEncoding="latin1")
data_iso <- read.csv("data/raw-data/isopoda.csv", h = T, sep = ";", fileEncoding="latin1")
data_diplo <- read.csv("data/raw-data/diplopoda.csv", h = T, sep = ";", fileEncoding="latin1")
data_coll <- read.csv("data/raw-data/collembola.csv", h = T, sep = ";", fileEncoding="latin1")

list_hab <- read.csv("data/derived-data/liste_habitat.csv", h = T, sep = ";") %>%
  rename_with(toupper) 

source("analyses/0. fonctions maison/DIVERSITY.R")
divent_car_ba <- DIVERSITY(data = data_car[data_car$METHODE == "Pitfall",])
divent_ara_ba <- DIVERSITY(data = data_ara[data_ara$METHODE == "pitfall",])
divent_iso <- DIVERSITY(data = data_iso)
divent_iso_tm <- DIVERSITY(data = data_iso[data_iso$METHODE == "trimanuel",])
divent_iso_ba <- DIVERSITY(data = data_iso[data_iso$METHODE == "Pitfall",])
divent_diplo <- DIVERSITY(data = data_diplo)
divent_diplo_tm <- DIVERSITY(data = data_diplo[data_diplo$METHODE == "trimanuel",])
divent_diplo_ba <- DIVERSITY(data = data_diplo[data_diplo$METHODE == "Pitfall",])
divent_coll <- DIVERSITY(data = data_coll)

# deficit de couverture multitaxa
def_cov <- rbind(cbind(taxa = "carabidaeBA", divent_car_ba$def_cov), 
                 cbind(taxa = "araneaeBA", divent_ara_ba$def_cov), 
                 cbind(taxa = "collembola", divent_coll$def_cov)) %>%
  left_join(divent_ara$station) %>%
  left_join(list_hab)

# Représentation du déficit global
ggplot(def_cov, aes(x= taxa, y = def_cov)) +
  geom_hline(aes(linetype = "5%"), yintercept = 5, 
             color = "#1b9e77", linewidth = 1) +
  geom_hline(aes(linetype = "33%"), yintercept = 33, 
                 color = "#d95f02") +
  geom_boxplot(outlier.shape = NA, colour = "violet")+
  lims(y=c(0,100))+
  geom_jitter(width = 0.15, alpha = 0.2, colour = "violet")+
  labs(y = "Deficit de couverture\n(% nombre d'espèces)", linetype = "")+
  theme_bw()

# Représentation du déficit par habitat
ggplot(def_cov, aes(x = taxa, y = def_cov)) +
  geom_hline(aes(linetype = "5%"), yintercept = 5, 
             color = "#1b9e77", linewidth = 1) +
  geom_hline(aes(linetype = "33%"), yintercept = 33, 
             color = "#d95f02") +
  geom_boxplot(outlier.shape = NA, colour = "violet")+
  geom_jitter(width = 0.15, alpha = 0.2, colour = "violet")+
  labs(y = "Deficit de couverture\n(% nombre d'espèces)", linetype = "")+
  lims(y=c(0,100))+
  facet_grid(.~HABITAT_CODE)+
  theme_bw()


# partition de diversité multitaxa
div_part <- rbind(cbind(taxa = "carabidaeBA", divent_car_ba$div_part), 
                 cbind(taxa = "araneaeBA", divent_ara_ba$div_part), 
                 cbind(taxa = "collembola", divent_coll$div_part)) %>%
  select(taxa, scale, estimator, diversity) %>%
  add_row(taxa = "carabidaeBA", scale = "gamma_obs", estimator = "", 
          diversity = length(unique(data_car$LB_NOM[data_car$METHODE == "Pitfall"]))) %>%
  add_row(taxa = "araneaeBA", scale = "gamma_obs", estimator = "", 
          diversity = length(unique(data_ara$LB_NOM[data_ara$METHODE == "pitfall"]))) %>%
  add_row(taxa = "collembola", scale = "gamma_obs", estimator = "", 
          diversity = length(unique(data_coll$LB_NOM))) %>%
  tibble
  

ggplot(div_part, aes(x = scale, y = log10(diversity)))+
  geom_bar(stat="identity")+
  facet_grid(taxa~., scales = "free") +
  labs(y = "Diversité \n(nombre d'espèces, log-transfromé)", x = "Echelles")+
  theme_bw()

# Completude
cov <- rbind(cbind(taxa = "carabidaeBA", divent_car_ba$cov), 
             cbind(taxa = "araneaeBA", divent_ara_ba$cov), 
             cbind(taxa = "collembola", divent_coll$cov)) %>%
  left_join(list_hab) %>%
  filter(!is.na(HABITAT_CODE))

ggplot(cov, aes(x = HABITAT_CODE, y = coverage)) +
  geom_boxplot()+
  labs(x = "habitats", y = "complétude (%)")+
  facet_grid( ~ taxa, scales = "free")+
  theme_bw()

# Species accumulation curves
acc  <- rbind(cbind(taxa = "carabidaeBA", divent_car_ba$acc), 
                     cbind(taxa = "araneaeBA", divent_ara_ba$acc), 
                     cbind(taxa = "collembola", divent_coll$acc)) %>%
  left_join(divent_ara$station) %>%
  left_join(list_hab) %>% 
  select(site, taxa, level, diversity, STATION, HABITAT_CODE) %>%
  left_join(cov[,1:3]) %>%
  group_by(taxa, STATION) %>%
  filter(level < weight)

# représentation
ggplot(acc, aes(x = level, y = diversity, group = STATION)) +
  geom_path(aes(colour = STATION))+
  facet_grid(taxa ~HABITAT_CODE, scales = "free")+
  theme_bw()

# profil de diversité
Hill_number <- rbind(cbind(taxa = "carabidaeBA", divent_car_ba$hill_prof), 
             cbind(taxa = "araneaeBA", divent_ara_ba$hill_prof), 
             cbind(taxa = "collembola", divent_coll$hill_prof)) %>%
  left_join(divent_ara$station) %>%
  left_join(list_hab) %>%
  filter(order %in% c(0,1,2))

ggplot(Hill_number, aes(x = HABITAT_CODE, y = diversity)) +
  geom_boxplot(outlier.shape = NA, colour = "violet")+
  geom_jitter(width = 0.15, alpha = 0.2, colour = "violet")+
  labs(x = "Nombre de Hill", y = "Diversité")+
  facet_wrap(as.factor(order)~taxa, scales = "free")+
  theme_bw()

