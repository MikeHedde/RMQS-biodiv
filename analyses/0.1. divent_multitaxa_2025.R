
librarian::shelf(dplyr, tidyr, ggplot2, vegan, forcats, divent, ggdist)

data_car <- read.csv("data/raw-data/lab files/databases update/update_carabidae.csv", h = T, sep = ";", fileEncoding="latin1")
data_ara <- read.csv("data/raw-data/lab files/databases update/update_araneae.csv", h = T, sep = ";", fileEncoding="latin1")

X <- data_car

PB3 <- X %>%
  filter(RANG == "ES",
         PROJET == "RMQS_2025",
         METHODE == "Pitfall") %>%
  mutate(METHODE = "Pitfall 1-3") %>%
  filter(as.numeric(substr(ID_ECHANTILLON, nchar(ID_ECHANTILLON), nchar(ID_ECHANTILLON))) %in% 1:3) %>%
  bind_rows(X)

PB6 <- X %>%
  filter(RANG == "ES",
         PROJET == "RMQS_2025",
         METHODE == "Pitfall") %>%
  mutate(METHODE = "Pitfall 1-6") %>%
  filter(as.numeric(substr(ID_ECHANTILLON, nchar(ID_ECHANTILLON), nchar(ID_ECHANTILLON))) %in% 1:6) %>%
  bind_rows(PB3)

mymat0 <- PB6 %>%
  filter(RANG == "ES",
         PROJET == "RMQS_2025") %>%
  group_by(NOM_STATION, METHODE, LB_NOM) %>%
  summarise(abtot = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(names_from = LB_NOM, values_from = abtot, values_fill = 0) %>%
  mutate_if(is.numeric, ~replace_na(., 0))


# Étape 2. Calcul du taux d'accumulation d'espèces
abd <- as_abundances(mymat0[,3:ncol(mymat0)])

abd_more_than_1_sp <- abd %>%
  filter(weight>1)

acc0 <- accum_hill(
  abd_more_than_1_sp, 
  q = 0,
  n_simulations = 10)  # l'augmentation du nombre de simulations diminue la variabilité de l'enveloppe

# Species accumulation curves
site_names <- tibble("site" = abd$site,
      "methode" = mymat0$METHODE,
      "station" = mymat0$NOM_STATION)

acc <- left_join(acc0, site_names) %>%
  filter(estimator != "Sample") %>%
  mutate(methode = fct_recode(methode,
                              "Pitfall 1-10" = "Pitfall"))
acc$methode <- factor(acc$methode, levels=c('GPD', "Pitfall 1-10", "Pitfall 1-6", "Pitfall 1-3",
                                            "Tri manuel", "DVAC"))

acc_sub <- acc %>%
  filter(station %in% c("RMQS_25_SAINT-GERMAIN_637"))

      # représentation
      ggplot(acc, aes(x = level, y = diversity, group = methode)) +
        geom_path(aes(colour = methode, linewidth = estimator))+
        labs(x="Nb of individuals (log10)", y = "Species Number (log10)", title =  "Araneae")+
        facet_wrap(station~.)+
        scale_x_log10()+        scale_y_log10()+
        #lims(x = c(0,300))+
        theme_bw()+
        scale_color_manual(values = c("GPD" = "black",
                                      "Pitfall 1-10" = "mediumpurple4", 
                                      "Pitfall 1-6" = "mediumpurple3",
                                      "Pitfall 1-3" = "mediumpurple1",
                                      "Tri manuel" = "aquamarine",
                                      "DVAC" = "chartreuse"))+
        scale_discrete_manual("linewidth", values = c(1.5, 0, .5))


      
acc_sample <- left_join(acc0, site_names) %>%
  filter(estimator == "Sample") %>%
  mutate(methode = fct_recode(methode,
                              "Pitfall 1-10" = "Pitfall"))
acc_sample$methode <- factor(acc_sample$methode, levels=c('GPD', "Pitfall 1-10", "Pitfall 1-6", "Pitfall 1-3",
                                                          "DVAC", "Tri manuel"))                                            

ggplot(acc_sample, aes(x = methode, y = diversity, fill = methode, alpha = 0.5))+
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA) +
  stat_slab(linewidth = 3)+
  theme_bw()
