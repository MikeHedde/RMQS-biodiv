# Workflow taxonomic uncertainty — Collembola

Ce dossier remplace le gros script monolithique par un workflow relançable par étapes.
Les figures ne sont plus produites dans les scripts d'analyse : elles sont uniquement dans `figures.R`.

## Structure

- `R/00_config.R` : paramètres, chemins, options, système de cache RDS.
- `R/functions_taxonomic_uncertainty.R` : fonctions partagées.
- `01_import_clean.R` : import, nettoyage, audit taxonomique initial.
- `02_prepare_taxonomy_pools.R` : TAXREF mainland, espèces observées, pools de confusion.
- `03_prepare_environment_drivers.R` : variables environnementales et drivers focaux.
- `04_build_scenarios.R` : scénarios déterministes et stochastiques.
- `05_compute_alpha_metrics.R` : alpha diversité, couverture, gamma par scénario.
- `06_compute_stability.R` : stabilité des inférences vs baseline.
- `07_compute_drivers.R` : PERMANOVA et modèles de drivers alpha.
- `08_rare_and_inventory_audits.R` : contribution des espèces rares et inventaires par scénario.
- `figures.R` : figures publication-ready uniquement.
- `run_all.R` : lance tout dans l'ordre.

## Utilisation standard

Place ce dossier au même niveau que ton dossier `data/`, ou définis explicitement le dossier projet :

```r
options(taxo.project_dir = "C:/chemin/vers/ton/projet")
source("run_all.R")
```

Les sorties CSV restent écrites dans :

```r
outputs_taxonomic_uncertainty_v9/
```

Les caches intermédiaires sont écrits dans :

```r
outputs_taxonomic_uncertainty_v9/_rds_cache/
```

Les figures publication-ready sont écrites dans :

```r
outputs_taxonomic_uncertainty_v9/figures_publication_ready/
```

## Relancer seulement une étape

Par défaut, si le cache RDS existe, le script recharge l'étape au lieu de la recalculer.
Pour forcer une seule étape :

```r
options(taxo.force_steps = "06_compute_stability")
source("06_compute_stability.R")
source("figures.R")
```

Pour forcer plusieurs étapes :

```r
options(taxo.force_steps = c(
  "04_build_scenarios",
  "05_compute_alpha_metrics",
  "06_compute_stability",
  "07_compute_drivers",
  "08_rare_and_inventory_audits"
))
source("run_all.R")
```

Pour tout recalculer :

```r
options(taxo.force_rebuild = TRUE)
source("run_all.R")
```

## Logique des dépendances

Chaque script recharge les caches dont il dépend. Par exemple, `06_compute_stability.R` recharge `04_build_scenarios.rds`, et `figures.R` recharge les sorties nécessaires sans relancer les calculs lourds.

## Note importante

Le script d'origine contenait aussi un bloc `10. FIGURES` avec des figures exploratoires. Il n'a pas été conservé dans le workflow d'analyse, pour respecter la consigne : les figures sont uniquement dans `figures.R`.
