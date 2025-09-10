
librarian::shelf(dplyr, tidyr, ggplot2, vegan, forcats, divent)

data_car <- read.csv("data/raw-data/lab files/databases update/carabidae_250906.csv", h = T, sep = ";", fileEncoding="latin1")
data_ara <- read.csv("data/raw-data/lab files/databases update/araneae_250906.csv", h = T, sep = ";", fileEncoding="latin1")

X <- data_car

mymat0 <- X %>%
  #filter(RANG == "ES") %>%
  group_by(STATION, METHODE, LB_NOM) %>%
  summarise(abtot = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(names_from = LB_NOM, values_from = abtot, values_fill = 0) %>%
  mutate_if(is.numeric, ~replace_na(., 0))

# Étape 2. Calcul du taux d'accumulation d'espèces
abd <- as_abundances(mymat0[,3:ncol(mymat0)])

abd_more_than_1_sp <- abd %>%
  filter(weight>1)

acc <- accum_hill(
  abd_more_than_1_sp, 
  q = 0,
  n_simulations = 10)  # l'augmentation du nombre de simulations diminue la variabilité de l'enveloppe

# Species accumulation curves
site_names <- tibble("site" = abd$site,
      "methode" = mymat0$METHODE,
      "station" = mymat0$STATION)

acc <- left_join(acc, site_names) %>%
  filter(estimator == "Interpolation")

      # représentation
      ggplot(acc, aes(x = level, y = diversity, group = station)) +
        geom_path(aes(colour = station))+
        labs(x="Nb of individuals", y = "Species Number", title =  "Carabidae")+
        facet_grid(methode ~ .)+
        lims(y=c(0,20))+
        theme_bw()+
        guides(colour="none")
      

  