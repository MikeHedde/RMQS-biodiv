library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggpubr)
library(rstatix)
library(broom)


setwd("C:/Users/Thoma/OneDrive/Bureau/R stage/SCRIPT moi")
masse <- read.csv("C:/Users/Thoma/OneDrive/Bureau/R stage/SCRIPT moi/clean.csv", h = TRUE, sep = ";", fileEncoding = "ISO-8859-1")
masse <- read.csv("data/raw-data/stage_Thomas_241126.csv", h=T, sep=";", dec=".", fileEncoding = "ISO-8859-1")


library(dplyr)

# Liste des mots à exclure
mots_a_exclure <- c("endommagé", "a peser", "non", "abdomen creux", "corps coupé en 2 ", "abdomen absent","repeser", "absence d'élytres ", "vide")  # Ajoutez d'autres mots ici

# Création de l'objet mass sans les lignes contenannt les mots a exclures car ils indiquent une masse qui n'est pas représentative
mass <- masse %>%
  filter(!grepl(paste(mots_a_exclure, collapse = "|"), remarque, ignore.case = TRUE))







# Convertir la colonne 'MASSE (mg)' en numérique --> permet de remplacer toutes les , par .
mass$MASSE <- as.numeric(gsub(",", ".", mass$MASSE))




# mass par esp
mass_esp <- mass %>%
  select(NOM_VALIDE, MASSE) 

nb_mass <- mass_esp %>%
  group_by(NOM_VALIDE)%>%
  summarise(nb = n())%>%
  filter(nb>1)

mass_esp2 <- mass_esp%>%
  filter(NOM_VALIDE %in% nb_mass$NOM_VALIDE)

ggplot(mass_esp2, aes(x = reorder(NOM_VALIDE, MASSE, median), y = MASSE))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.5, width = 0.25, colour = "lightgreen")+
  labs(y = "Masse (mg)", x ="")+
  scale_y_continuous(trans='log10')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # fond cadrillage ? 

###############################################################################################""
 sp_list = c("Poecilus cupreus (Linnaeus, 1758)", 
            "Harpalus distinguendus (Duftschmid, 1812)",
            "Amara aenea (De Geer, 1774)", 
            "Carabus auratus Linnaeus, 1761",
            "Calathus fuscipes (Goeze, 1777)",
            "Abax parallelepipedus (Piller & Mitterpacher, 1783)",
            "Carabus monilis Fabricius, 1792",
            "Carabus nemoralis")
#Lignes a garder uniquement si je veux faire une liste a un momment 



# Filtrer les données selon la liste d'espèces et utiliser "HABITAT" pour la couleur et la comparaison
SP <- mass %>%
  filter(NOM_VALIDE %in% sp_list)  # Garder seulement les espèces de sp_list


#Version plus simple : 
#Boxplot  milieux ouvert / fermé de la liste 


#carabe_mass_individual <- mass %>%
  filter(NOM_VALIDE %in% sp_list) %>%   # Garder les espèces de sp_list (carabes)
  mutate(MASSE = as.numeric(MASSE))     # Convertir MASSE en numérique si nécessaire

# Créer le boxplot avec les masses individuelles par habitat de la liste 
#ggplot(carabe_mass_individual, aes(x = HABITAT, y = MASSE, fill = HABITAT)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.5, width = 0.25, colour = "lightgreen") + # Ajouter des points pour voir les valeurs individuelles
  scale_y_continuous(trans = 'log10')
  labs(y = "Masse individuelle des carabes (mg)", x = "Type de milieu") +
  theme_bw()
#########################################################################################""

# Créer le boxplot avec les masses individuelles par habitat (sans filtrer par sp_list)
ggplot(mass, aes(x = HABITAT, y = MASSE, fill = HABITAT)) +
  geom_violin()+
  geom_boxplot(outlier.shape = NA, alpha = 0.1) + #masque valeurs abérantes
  geom_jitter(alpha = 0.25, width = 0.25, colour = "lightgreen") + # Ajouter des points pour voir les valeurs individuelles; alpha =opacité width= dispersition evite superposition des points 
  scale_y_continuous(trans = 'log10') + # Appliquer la transformation logarithmique sur l'axe des y
  labs(y = "Masse individuelle des carabes (mg)", x = "Type de milieu") + #titre des axes
  theme_bw()#ajoute un fond = meilleur visibilité du graph 


#Constat probable différence avec milieu fermé abritant en moyenne + d'espèces "grosse" mais besoin test stats
# car boites se chevauchent

#MEME GRAPH QUE PRECEDANT MAIS PAR LOCALIE :
#Changement d'echelle afin d'aténuer les différences d'individus 

# Créer un nouveau tableau avec la somme des masses par localité
mass_summed <- mass %>%
  group_by(NOM_LOCALITE, HABITAT) %>%
  summarise(SOMME_MASSE = sum(MASSE, na.rm = TRUE)) # na.rm pour ignorer les NA
library(ggplot2)

# Créer le graphique avec les masses additioné par localité 
ggplot(mass_summed, aes(x = HABITAT, y = SOMME_MASSE, fill = HABITAT)) +
  #geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.1) + # Masquer les valeurs aberrantes
  geom_jitter(alpha = 0.25, width = 0.25, colour = "lightgreen") + # Ajouter des points pour les localités
  scale_y_continuous(trans = 'log10') + # Transformation logarithmique pour l'axe des y
  labs(y = "Somme des masses des carabes (mg)", x = "Type de milieu") + # Titre des axes
  theme_bw() # Fond plus visible



#Voir si données suivent une loie normale :
shapiro.test(mass$MASSE)
#p value =2.2e-16 < 0.05 donc les données ne suivent pas une loie normale 
#L'équivalent non paramétrique de l'ANOVA a 1 facteur est : le test de Kruskal-Wallis mais k<3
# donc test Wilcoxon : 

## NON, ce ne sont pas les données qui doivent suivre
## une loi normale. Ce sont les résidus de la relation
#1. Transformer les colonnes `male` et `femelle` en une colonne `sex`
mass_long <- mass %>%
  filter(!is.na(MASSE)) %>%  # Enlever les lignes sans masse
  mutate(sex = case_when(
    male == 1 ~ "male",  # Si la colonne 'male' vaut 1, on attribue "male" à la nouvelle colonne 'sex'
    femelle == 1 ~ "female",  # Si la colonne 'femelle' vaut 1, on attribue "female" à la nouvelle colonne 'sex'
    TRUE ~ NA_character_  # Sinon, on assigne NA
  )) %>%
  filter(!is.na(sex))  # Supprimer les lignes sans information de sexe


#2. Calculer la masse moyenne par station, habitat, espèce, et sexe
rensch_wide <- mass_long %>%
  group_by(NOM_LOCALITE, HABITAT, NOM_VALIDE, sex) %>%
  summarise(mean_mass = mean(MASSE, na.rm = TRUE)) %>% #Calcule moyenne de chaque groupe
  filter(all(c("male", "female") %in% sex)) %>%  # Garder les localité avec les deux sexes présents ; c() permet de créer un vecteur 
  pivot_wider(names_from = sex, values_from = mean_mass) %>%  # Étendre pour avoir une colonne `male` et `female`
  mutate(log_female = log(female), log_male = log(male))

ggscatter(rensch_wide, x = "log_male", y = "log_female",
  color = "HABITAT", add = "reg.line")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = HABITAT)
  )
# Milieu ouvert présente un R2 à 0.93, plus faible que en fermé. 
# Certainement du à des outliers qu'il va falloir identifier et corriger

#3. Normalité des résidus et homogénéité des variances

  # Calculer le modèle, la covariable passe en premier
  model <- lm(log_female ~ log_male + HABITAT, data = rensch_wide)
  # Inspecter les paramètres de diagnostic du modèle
  model.metrics <- augment(model) %>%
    select(-.hat, -.sigma, -.fitted) # Supprimer les détails
  head(model.metrics, 3)

      # Évaluer la normalité des résidus à l'aide du test de Shapiro-Wilk
      shapiro_test(model.metrics$.resid)
      #Le test de Shapiro-Wilk est  significatif (p > 0,05), 
      # résidus ne suivent pas une distribution normale

            # Homogénéité des variances
            model.metrics %>% levene_test(.resid ~ HABITAT)
            # Le test de Levene n’était pas significatif (p > 0,05), 
            # homogénéité des variances résiduelles pour tous les groupes

                #Valeurs aberrantes
                model.metrics %>% 
                  filter(abs(.std.resid) > 2) %>%
                  as.data.frame()
            

#4. ANCOVA (si normalité !)
    res.lm <- lm(log_female ~ log_male + HABITAT, data = rensch_wide)
    tidy(summary(res.lm))
    Anova(res.lm)
    #peu de dimorphisme sexuel (pente fermé = 0.989, soit quasi 1)
    #pas de différence de dimorphisme sexuel significatif (p=0.099) entre milieux
            
# Filtrer les données pour les deux milieux
mass_ouvert <- mass %>% filter(HABITAT == "ouvert") %>% pull(MASSE)
mass_ferme <- mass %>% filter(HABITAT == "ferme") %>% pull(MASSE)

# Exécuter le test de Mann-Whitney
wilcox.test(mass_ouvert, mass_ferme)
#On a p value = 2.2e-16 <0.05 donc il y a une différence significative entre les masses des carabidae de milieux ouverts et fermes 
# Le milieu a donc bien un effet sur la masse. Les carabidae des milieux fermes ont une masse plus importante que ceux des milieux ouvert



##################################

#Le milieu a t-il une influance sur le dimorphisme sexuel de masse chez les carabidae ? 

# Transformer les colonnes `male` et `femelle` en une colonne `sex`
mass_long <- mass %>%
  filter(!is.na(MASSE)) %>%  # Enlever les lignes sans masse
  mutate(sex = case_when(
    male == 1 ~ "male",  # Si la colonne 'male' vaut 1, on attribue "male" à la nouvelle colonne 'sex'
    femelle == 1 ~ "female",  # Si la colonne 'femelle' vaut 1, on attribue "female" à la nouvelle colonne 'sex'
    TRUE ~ NA_character_  # Sinon, on assigne NA
  )) %>%
  filter(!is.na(sex))  # Supprimer les lignes sans information de sexe


# Calculer la masse moyenne par station, habitat, espèce, et sexe
rensch_wide <- mass_long %>%
  group_by(NOM_LOCALITE, HABITAT, NOM_VALIDE, sex) %>%
  summarise(mean_mass = mean(MASSE, na.rm = TRUE)) %>% #Calcule moyenne de chaque groupe
  filter(all(c("male", "female") %in% sex)) %>%  # Garder les localité avec les deux sexes présents ; c() permet de créer un vecteur 
  pivot_wider(names_from = sex, values_from = mean_mass)  # Étendre pour avoir une colonne `male` et `female`

# Créer le graphique avec la règle de Rensch
ggplot(rensch_wide, aes(x = log(male), y = log(female), colour = HABITAT)) +
  geom_point() +  # Points colorés par type de milieu
  geom_smooth(method = "glm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE) +  # Régression quadratique
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Bissectrice
  labs(x = "Log(Male Mass)", y = "Log(Female Mass)", title = "Rensch's Rule Representation") +
  theme_minimal()
#On constate que les femelle seraient plus grosse que les malles, chez les petites sp, dans les milieux fermés
# probablement : pente < 1 indique une allométrie conforme a la règle de Resch 
#En milieu ouvert, la pente semble +/- = 1 isométrie les males et les femelles augmentent de la meme proportion ?  

#Aalyse stat: 




# Créer une nouvelle table avec les masses log-transformées +moyenne des logs pour les males et les femelles

rensch_data <- mass %>%
  filter(!is.na(MASSE), male == 1 | femelle == 1) %>%
  mutate(
    log_male = ifelse(male == 1, log(MASSE), NA),
    log_female = ifelse(femelle == 1, log(MASSE), NA)
  ) %>%
  group_by(HABITAT, NOM_VALIDE) %>%
  summarise(
    log_male = mean(log_male, na.rm = TRUE),
    log_female = mean(log_female, na.rm = TRUE)
  ) %>%
  filter(!is.na(log_male) & !is.na(log_female))

# Modèle ANCOVA
model1 <- lm(log_female ~ log_male + HABITAT, data = rensch_data)
summary(ancova_model)
#On a p value significative que sur l'effet de la masse des males sur celui des femelles ( j'ai du male a formulé et comprendre cette partie )

#Vérification du model : 

#On teste la normalité des résidus : 
shapiro.test(model1$residuals)
#On a pvalue = 1.032e-09 donc < 0.05 donc les résidus ne sont pas normaux et on ne peut pas appliquer une ANCOVA



