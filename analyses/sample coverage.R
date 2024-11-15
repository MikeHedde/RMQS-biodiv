# Chargement des packages nécessaires
install.packages("dplyr")
install.packages("vegan")
install.packages("ggplot2")
install.packages("iNEXT")
install.packages("ggiNEXT")
library(iNEXT)
library(dplyr)
library(vegan)
library(ggplot2)
library(tibble)

# Chargement des données
sandbox <- read.csv("data/raw-data/sandbox.csv", h = T, sep = ";", encoding = "Latin-1")

# Fonction pour la conversion en data.frame avec rownames
convert_to_data_frame_with_rownames <- function(tibble, rowname_column) {
  # Convertir en data.frame
  df <- as.data.frame(tibble)
  
  # Vérifier que la colonne existe et l'utiliser comme rownames
  if (rowname_column %in% names(df)) {
    rownames(df) <- make.unique(as.character(df[[rowname_column]]))
    df[[rowname_column]] <- NULL  # Retirer la colonne pour éviter la duplication
  } else {
    stop("La colonne spécifiée pour les rownames n'existe pas.")
  }
  
  return(df)
}

# Transformation des données et calcul du sampling coverage
mymat <- sandbox %>%
  filter(METHODE == "Pitfall") %>%
  filter(!ORDRE %in% c("VIDE", "MANQUANT")) %>%
  pivot_wider(id_cols = c(STATION, ORDRE), names_from = REPETITION, values_from = ABONDANCE_TOTALE, 
              values_fn = sum, values_fill = 0) %>%
  group_split(STATION, .keep = F)

mymat_dataframes <- lapply(mymat, function(tbl) convert_to_data_frame_with_rownames(tbl, "ORDRE"))
lapply(mymat_dataframes, class)

sample_coverage_results <- iNEXT(mymat_dataframes, q = 0, datatype = "incidence_raw", endpoint = NULL)


# 4. Représentation
ggplot(sample_coverage_results$DataInfo, aes(y = SC))+
  geom_boxplot()+
  labs(y = "Sample coverage (%)")+
  theme_minimal()+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) 
  
