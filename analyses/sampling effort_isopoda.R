librarian::shelf(dplyr, tidyr, vegan, nls2, ggplot2, pacman, tidytable, iNEXT.4steps)


data <- read.csv("data/raw-data/isopoda.csv", h = T, sep = ";", fileEncoding="latin1")
list_hab <- read.csv("data/derived-data/liste_habitat.csv", h = T, sep = ";") %>%
  rename_with(toupper) 


# Étape 1. Attribution des habitats standardisés
data <- data %>%
  left_join(list_hab) %>%
  filter(ORDRE == "Isopoda")

# Étape 2. Création d'une matrice de présence/absence
# Regroupement des données par les colonnes spécifiées et calcul de l'abondance totale
# Transformation en matrice de présence/absence (1 = présence, 0 = absence)
mymat <- data %>%
  group_by(HABITAT, METHODE, STATION, REPETITION, LB_NOM) %>%
  summarise(abtot = sum(ABONDANCE_TOTALE), .groups = 'drop') %>%
  mutate(presence = ifelse(abtot > 0, 1, 0)) %>%
  select(-abtot) %>%
  # Ajout des combinaisons manquantes pour répétitions
  complete(METHODE, STATION, REPETITION := 1:6, fill = list(presence = 0)) %>%
  # Transformation au format large avec chaque espèce comme colonne
  pivot_wider(names_from = LB_NOM, values_from = presence, values_fill = 0) %>%
  # Remplissage des informations manquantes sur les habitats
  group_by(STATION) %>%
  fill(HABITAT, .direction = "downup") %>%
  ungroup() %>%
  mutate(met_stat = paste(METHODE, STATION, sep = "_"))

# Étape 2 : Calcul des courbes d'accumulation
# Initialisation d'une liste pour stocker les résultats des courbes d'accumulation
accumulation_data <- list()

# Boucle pour chaque combinaison de méthode et de site
for (met_stat in unique(mymat$met_stat)) {
    # Sous-ensemble des données pour une méthode et un site donnés
    site_data <- mymat %>%
      filter(met_stat == met_stat) %>%
      select(-c(STATION, REPETITION, METHODE, HABITAT, met_stat))
    
    # Calcul de la courbe d'accumulation
    site_accumulation <- specaccum(site_data, method = "random")
    
    # Stockage des résultats dans un dataframe
    accumulation_data[[met_stat]] <- data.frame(
      #methode = method,
      #station = station,
      sites = site_accumulation$sites,
      richness = site_accumulation$richness
    )
  }

# Conversion des résultats en un dataframe unique
accumulation_df <- bind_rows(accumulation_data, .id = "column_label") %>%
  separate(column_label, into = c("METHODE", "STATION"), sep = "_") %>%
  rename_with(toupper) %>%
  rownames_to_column() %>%
  mutate(STATION = as.integer(STATION)) %>%
  left_join(list_hab)


# Fonction pour la conversion en data.frame avec rownames
convert_to_data_frame_with_rownames <- function(tibble, rowname_column) {
  # Convertir en data.frame
  df <- as.data.frame(tibble)
  
  # Vérifier que la colonne existe et l'utiliser comme rownames
  if (rowname_column %in% names(df)) {
    rownames(df) <- make.unique(as.character(df[[rowname_column]]))
    df[[rowname_column]] <- NULL  # Retirer la colonne pour éviter la duplication
  } else {
    stop("convert_to_data_frame_with_rownames : La colonne spécifiée pour les rownames n'existe pas")
  }
  
  return(df)
}


# Transformation des données et calcul du sampling coverage
mymat2 <- mymat %>%
  select(-c(HABITAT, METHODE, STATION)) %>% 
  relocate(met_stat) %>%
  # Retirer les colonnes inutiles
  pivot_longer(cols = c(3:ncol(.)),                  # Long format (de la 3ème à dernière colonne)
               names_to = "esp", 
               values_to = "presence") %>%
  filter(presence > 0) %>%                           # Filtrer pour les présences > 0
  group_by(met_stat, REPETITION) %>%                  # Regrouper par station et répétition
  filter(length(unique(.)) > 1) %>%                  # Garder les groupes avec plus d'une ligne
  ungroup() %>%                                      # Dégroupement
  pivot_wider(names_from = REPETITION,               # Transformer au format large
              values_from = presence, 
              values_fill = 0,
              values_fn = max) %>%
  tidytable::group_split(met_stat, .keep = FALSE, 
                         .named = TRUE) %>%                    # Diviser par station
  keep(~ nrow(.) >= 3) %>%                           # Garder les dataframes avec au moins 2 lignes
  keep(~ sum(rowSums(.[, -1]) > 1) > 0)              # Vérifier si plus d'une colonne a des "1"


mymat_dataframes <- lapply(mymat2, function(tbl) convert_to_data_frame_with_rownames(tbl, "esp"))
lapply(mymat_dataframes, class)

iNEXT_results <- iNEXT4steps(data = mymat_dataframes, datatype = "incidence_raw")
iNEXT_results2 <- iNEXT(mymat_dataframes[[21]], q = 0, datatype = "incidence_raw", endpoint = NULL)

SamplCov <- Completeness(data = mymat_dataframes, datatype = "incidence_raw")
SamplCov <- SamplCov %>%
  rename(STATION = Assemblage) %>%
  separate(STATION, into = c("METHODE", "STATION"), sep = "_") %>%
  mutate(STATION = as.integer(STATION)) %>%
  left_join(list_hab)

# Relation abondance-richesse  
mymat3 <- data %>%
  group_by(METHODE, STATION, REPETITION, LB_NOM) %>%
  summarise(ABONDANCE_TOTALE = sum(ABONDANCE_TOTALE), .groups = 'drop') %>%
  complete(METHODE, STATION, REPETITION = 1:6,  # Toutes les répétitions possibles
           fill = list(PRESENCE = 0)) %>%  # Valeurs par défaut pour PRESENCE
  pivot_wider(names_from = LB_NOM, values_from = ABONDANCE_TOTALE, values_fill = 0) %>%
  left_join(select(list_hab, c(STATION, HABITAT))) %>% 
  group_by(STATION) %>%
  fill(HABITAT, .direction = "downup") %>%  # Remplit vers le bas puis le haut
  ungroup()

numeric_columns <- mymat %>%
  select(-STATION, -REPETITION) %>%
  select(where(is.numeric))%>%
  colnames()

mymat3 <- mymat3 %>%
  mutate(
    AB = rowSums(select(., all_of(numeric_columns)), na.rm = TRUE),  # Abondance totale
    SR = specnumber(select(., all_of(numeric_columns)))) %>%              # Richesse spécifique
  select(STATION, REPETITION, HABITAT, AB, SR)



# Etape 3. Visualisation des résultats
# Graphique des courbes d'accumulation
ggplot(subset(accumulation_df, !is.na(HABITAT)), 
       aes(x = SITES, y = RICHNESS, 
           group = STATION, colour = HABITAT)) +
  geom_line() +
  lims(x = c(1,6), y = c(0,5))+
  facet_wrap(HABITAT ~.) +
  labs(x = "Nombre d'échantillons",
       y = "Richesse spécifique") +
  theme_bw()+
  theme(legend.position = "none")

# sampling coverage
ggplot(SamplCov, aes(x = Order.q, y = Estimate.SC, 
                     group = STATION, colour = HABITAT))+
  geom_line(linewidth = 1)+
  labs(title = "Complétude d'échantillonnage par station",
       x = "valeur de q",
       y = "Sample coverage (%)",
       color = "Type d'habitat") +
  lims(y = c(0, 1))+
  facet_wrap(~HABITAT)+
  theme_bw()+
  theme(legend.position = "none")

# sampling coverage moy pour q=0
ggplot(subset(SamplCov, Order.q %in% c(0, 1, 2)), aes(x = HABITAT, y = Estimate.SC, 
                                                      colour = HABITAT))+
  geom_boxplot()+
  labs(title = "Complétude d'échantillonnage par habitat",
       x = "valeur de q",
       y = "Sample coverage (%)",
       color = "Type d'habitat") +
  facet_grid(Order.q~.)+
  lims(y = c(0, 1))+
  theme_bw()+
  theme(legend.position = "none")

# ab - richesse
ggplot(subset(mymat3, !is.na(HABITAT)), 
       aes(x= log1p(AB), y = SR, colour = HABITAT))+
  geom_point()+
  lims(y=c(0,3))+
  labs(y = "Richesse spécifique \n(log)", 
       x = "Nombre d'individus collectés \n(log)")+
  geom_smooth(method = lm, se = T)+
  facet_wrap(HABITAT~.)+
  theme_bw()
