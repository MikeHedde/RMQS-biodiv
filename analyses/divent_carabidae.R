librarian::shelf(dplyr, tidyr, vegan, nls2, ggplot2, pacman, tidytable, iNEXT.4steps, purrr)
install.packages('divent', repos = c('https://ericmarcon.r-universe.dev', 'https://cloud.r-project.org'))

data <- read.csv("data/raw-data/carabidae.csv", h = T, sep = ";", fileEncoding="latin1")
list_hab <- read.csv("data/derived-data/liste_habitat.csv", h = T, sep = ";") %>%
  rename_with(toupper) 



# Étape 1. Création d'une matrice d'abondance-activité

mymat <- data %>%
  mutate(site = paste(STATION, "_", REPETITION)) %>%
  group_by(site, LB_NOM) %>%
  summarise(abtot = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(names_from = LB_NOM, values_from = abtot, values_fill = 0)

mymat <- as_species_distribution(mymat)
abd_sum(mymat)

