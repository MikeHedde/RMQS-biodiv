# Simplification des types d'habitats en regroupant les catégories similaires

simplif_habitat <- function(data, habitat_column){
  data <- data %>%
    mutate(HABITAT = case_when(
      !!sym(habitat_column) %in% c("Monoculture", "Culture", "grande culture", "culture", "Culture de coriande", "Verger", "Culture de blé") ~ "culture",
      !!sym(habitat_column) %in% c("Chênaie claire", "Hetraie", "Forêt mixte", "Chenaie","Chênaie/ Érablaie", "Forêt") ~ "foret",
      !!sym(habitat_column) %in% c("Pelouse écorchée à salicorne", "Friche humide", "Lande à bruyère / fougeraie", "Pelouse alpine", "Parc urbain") ~ "autre_ouvert", 
      !!sym(habitat_column) %in% c("Prairie de fauche", "Prairie", "Prairie de pâture") ~ "prairie",
      TRUE ~ NA_character_
    ))
}
