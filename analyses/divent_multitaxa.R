librarian::shelf(ggplot2, dplyr, forcats, ade4)


data_car <- read.csv("data/raw-data/carabidae.csv", h = T, sep = ";", fileEncoding="latin1")
data_ara <- read.csv("data/raw-data/araneae.csv", h = T, sep = ";", fileEncoding="latin1")
data_iso <- read.csv("data/raw-data/isopoda.csv", h = T, sep = ";", fileEncoding="latin1")
data_diplo <- read.csv("data/raw-data/diplopoda.csv", h = T, sep = ";", fileEncoding="latin1")
data_coll <- read.csv("data/raw-data/collembola.csv", h = T, sep = ";", fileEncoding="latin1")
data_form <- read.csv("data/raw-data/formicidae.csv", h = T, sep = ";", fileEncoding="latin1")


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
divent_form_tm <- DIVERSITY(data = data_form[data_form$METHODE == "trimanuel",])
divent_form_ba <- DIVERSITY(data = data_form[data_form$METHODE == "Pitfall",])

# deficit de couverture multitaxa
def_cov <- rbind(cbind(taxa = "carabidaeBA", divent_car_ba$def_cov), 
                 cbind(taxa = "araneaeBA", divent_ara_ba$def_cov), 
                 cbind(taxa = "collembola", divent_coll$def_cov),
                 cbind(taxa = "formicidaeBA", divent_form_ba$def_cov)) %>%
  left_join(divent_ara$station) %>%
  left_join(list_hab)

    # fusion avec les degrés-jours
    clim_moy <- read.csv("data/derived-data/clim_moy.csv", h = T, sep = ",")
    pca_clim <- dudi.pca(clim_moy[-48, 3:6], scannf = F)
    pca_climli <- cbind(STATION = clim_moy$STATION[-48], pca_clim$li)
    threshold = 20
    
    def_cov <- left_join(def_cov, clim_moy) %>%
      left_join(pca_climli) %>%
      mutate(cov_grp = cut(def_cov,
                     breaks = c(0, threshold, 100),
                     labels = c("low", "high")))

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
        facet_grid(.~HABITAT_TYPE)+
        theme_bw()
      
      # Représentation du déficit par degre-jour cumulés
      ggplot(def_cov, aes(x = DD_cum, y = def_cov)) +
        geom_hline(aes(linetype = "5%"), yintercept = 5, 
                   color = "#1b9e77", linewidth = 1) +
        geom_hline(aes(linetype = "33%"), yintercept = 33, 
                   color = "#d95f02") +
        geom_point(size = 3, colour = "#6C3BAA", alpha = 0.5)+
        geom_smooth(se=FALSE, method = "lm", formula = y ~ poly(x, 1), lty = "dashed", color = "#6C3BAA")+
        labs(y = "Deficit de couverture\n(% nombre d'espèces)", 
             x = "Degrés-jours (>10 °C) cumulés\nau jour du prélèvement", linetype = "")+
        lims(y=c(0,100))+
        facet_grid(HABITAT_TYPE~taxa)+
        theme_bw()


      ggplot(def_cov, aes(x=Axis1, y = Axis2))+
        geom_point(aes(alpha = def_cov, size = def_cov, 
                       colour = HABITAT_TYPE))+
        facet_wrap(taxa~HABITAT_TYPE)+
        theme_bw()
      
      
# partition de diversité multitaxa
div_part <- rbind(cbind(taxa = "carabidaeBA", divent_car_ba$div_part), 
                 cbind(taxa = "araneaeBA", divent_ara_ba$div_part), 
                 cbind(taxa = "collembola", divent_coll$div_part), 
                 cbind(taxa = "formicidaeBA", divent_form_ba$div_part)) %>%
  select(taxa, scale, estimator, diversity) %>%
  add_row(taxa = "carabidaeBA", scale = "gamma_obs", estimator = "", 
          diversity = length(unique(data_car$LB_NOM[data_car$METHODE == "Pitfall"]))) %>%
  add_row(taxa = "araneaeBA", scale = "gamma_obs", estimator = "", 
          diversity = length(unique(data_ara$LB_NOM[data_ara$METHODE == "pitfall"]))) %>%
  add_row(taxa = "collembola", scale = "gamma_obs", estimator = "", 
          diversity = length(unique(data_coll$LB_NOM))) %>%
  add_row(taxa = "formicidaeBA", scale = "gamma_obs", estimator = "", 
          diversity = length(unique(data_form$LB_NOM))) %>%  
  tibble
  

ggplot(div_part, aes(x = scale, y = log10(diversity)))+
  geom_bar(stat="identity")+
  facet_grid(taxa~., scales = "free") +
  labs(y = "Diversité \n(nombre d'espèces, log-transfromé)", x = "Echelles")+
  theme_bw()

# Completude
cov <- rbind(cbind(taxa = "carabidaeBA", divent_car_ba$cov), 
             cbind(taxa = "araneaeBA", divent_ara_ba$cov), 
             cbind(taxa = "collembola", divent_coll$cov), 
             cbind(taxa = "formicidaeB1", divent_form_ba$cov)) %>%
  left_join(list_hab) %>%
  filter(!is.na(HABITAT_CODE))

        # fusion avec les degrés-jours
        cumul_dd <- read.csv("data/derived-data/cumul_dd.csv", h = T, sep = ",")
        cov <- left_join(cov, cumul_dd)

        ggplot(cov, aes(x = HABITAT_TYPE, y = coverage)) +
          geom_boxplot(outlier.shape = NA)+
          geom_jitter(colour = "violet", width = 0.25, size = 2)+
          labs(x = "habitats", y = "complétude (%)")+
          #facet_grid( ~ taxa, scales = "free")+
          theme_bw()
        
        # Représentation du déficit par degre-jour cumulés
        ggplot(cov, aes(x = cumulative_dd, y = coverage)) +
          geom_point(colour = "violet")+
          #geom_smooth(se=FALSE, method = "lm", formula = y ~ poly(x, 1), lty = "dashed", color = "#DEEBF7")+
          geom_smooth(se=T, method = "lm", formula = y ~ poly(x, 2), lty = "dashed", color = "#9ECAE1")+
          #geom_smooth(se=FALSE, method = "lm", formula = y ~ poly(x, 3), lty = "dashed", color = "#4292C6")+
          labs(y = "complétude (%)", 
               x = "Degrés-jours (>10 °C) cumulés\nau jour du prélèvement", linetype = "")+
          #lims(y=c(0,100))+
          #facet_grid(HABITAT_TYPE~.)+
          theme_bw()

# Species accumulation curves
acc  <- rbind(cbind(taxa = "carabidaeBA", divent_car_ba$acc), 
                     cbind(taxa = "araneaeBA", divent_ara_ba$acc), 
                     cbind(taxa = "collembola", divent_coll$acc), 
                     cbind(taxa = "formicidaeBA", divent_form_ba$acc)) %>%
  left_join(divent_ara$station) %>%
  left_join(list_hab) %>% 
  select(site, taxa, level, diversity, STATION, HABITAT_CODE, HABITAT_TYPE) %>%
  left_join(cov[,1:4], relationship = "many-to-many") %>%
  group_by(taxa, STATION) %>%
  filter(level < weight) %>%
  left_join(cumul_dd) 

      # représentation
      ggplot(acc, aes(x = level, y = diversity, group = STATION)) +
        geom_path(aes(colour = STATION))+
        facet_grid(taxa ~ HABITAT_CODE, scales = "free")+
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


# Effectifs
eff  <- acc %>%
  select(taxa, level, STATION, HABITAT_TYPE, weight, cumulative_dd)%>%
  group_by(STATION, taxa) %>%
  filter(level == max(level))

ggplot(eff, aes(x=as.factor(STATION), y = weight, colour = HABITAT_TYPE))+
  geom_point()+
  facet_grid(taxa~., scales = "free_y")+
  labs(x = "Stations",
       y = "Effectifs \n(nombre d'individus)")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

ggplot(subset(eff, weight>6), aes(cumulative_dd, weight, 
                fill=HABITAT_TYPE, colour = HABITAT_TYPE))+
  geom_point()+
  geom_smooth(se=F, method = "lm", formula = y ~ poly(x, 2), 
              lty = "dashed")+
  scale_y_log10()+
  scale_x_log10()+
  labs(x = "Degrés-jours (>10 °C) cumulés\nau jour du prélèvement",
       y = "Effectifs \n(nombre d'individus)")+
  facet_grid(taxa~.)+
  theme_bw()
  