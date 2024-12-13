# wrapper pour calculer plusieurs metriques de diversité avec divent (Marcon)

DIVERSITY <- function(data){
# Étape 1. Création d'une matrice d'abondance-activité
librarian::shelf(dplyr, tidyr, ggplot2, vegan, forcats, divent)

mymat0 <- data %>%
  group_by(STATION, LB_NOM) %>%
  summarise(abtot = sum(ABONDANCE_TOTALE)) %>%
  pivot_wider(names_from = LB_NOM, values_from = abtot, values_fill = 0) %>%
  mutate_if(is.numeric, ~replace_na(., 0))

# Étape 2. Calcul du taux de couverture
abd <- as_abundances(mymat0[,2:ncol(mymat0)])
obs_SR <- specnumber(mymat0[,2:ncol(mymat0)])
cov <- cbind(coverage(abd), STATION = mymat0$STATION)

# Etape 3. Courbes d'accumulation d'espèces
abd_more_than_1_sp <- abd %>%
  filter(weight>1)
acc <- accum_hill(
  abd_more_than_1_sp, 
  q = 0,
  n_simulations = 10)  # l'augmentation du nombre de simulations diminue la variabilité de l'enveloppe

# Etape 4. Profils de diversité
hill_prof <- profile_hill(abd)

# metacommunity
mc <- metacommunity(abd)
div_part0 <- div_part(abd, q = 0) 
div_part <- filter(div_part0, !scale == "community")

# Deficit de couverture
def_cov <- div_part0 %>%
  filter(scale == "community") %>%
  select(site, estimator, diversity) %>%
  mutate(obs_div = obs_SR) %>%
  mutate(def_cov = (diversity-obs_div)/diversity*100)

# station
station <- cov %>%
  select(site, STATION)

results <- list(station = station, cov = cov, acc = acc, hill_prof = hill_prof, 
                div_part = div_part, def_cov = def_cov)
return(results)
}