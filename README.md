# Devoir maison n°1 : Régressions Régularisées

## Description

Ce projet s'inscrit dans le cadre du cours **Analyse et Modélisation Multivariées** du Master 2 SSD de l'Université de Montpellier. L'objectif principal est de modéliser la densité globale de peuplement arboré de 27 espèces d'arbres dans le bassin du fleuve Congo à partir de variables géographiques, climatiques et de photosynthèse.

## Objectifs

Comparer différentes techniques de régression pour prédire la densité arborée :
- Régression sur Composantes Principales (RCP) globale
- Régression sur Composantes Principales thématiques (Géographie + Photosynthèse)
- Régression PLS (Partial Least Squares)
- Régression Ridge
- Régression LASSO

## Données

Les données proviennent du projet **CoForTips** et comprennent 1000 parcelles forestières du bassin du fleuve Congo avec :

### Variables explicatives
- **Variables géographiques** : latitude, longitude, altitude
- **Variables climatiques** : pluviométries annuelle et mensuelles
- **Variables qualitatives** : geology (type de sol)
- **Indices de photosynthèse** : 23 indices EVI

### Variable à prédire
- **Densité arborée** : somme des abondances des 27 espèces (gen1 à gen27) divisée par la surface de la parcelle

## Packages R

### Packages utilisés
- `FactoMineR` : Analyses en Composantes Principales
- `factoextra` : visualisations graphiques
- `pls` : régression PLS
- `glmnet` : régressions Ridge et LASSO

### Installation
```r
install.packages(c("FactoMineR", "factoextra", "pls", "glmnet"))
```

## Structure du projet

```
.
├── README.md
├── Analyse_multivariée.Rmd    # fichier de code qui effectue toutes les analyses
├── codes_archives/            # codes de test pour les analyses
│   ├── TP1_Regressions_regularisees.R
│   ├── Test Q4.Rmd
│   └── DM1.R  
├── Consignes.pdf              # consignes données pour ce devoir
├── Datagenus.csv              # données pour l'analyse                
└── figures/                   # dossier contenant les graphiques générés
    └── ...
```

## Utilisation

1. Cloner le dépôt :
```bash
git clone [https://github.com/EmilieAig/Analyse_multivariee]
cd Analyse_multivariee
```

3. Exécuter les scripts dans l'ordre (Q1 à Q6)

## Auteurs

AIGOIN Emilie
THOMAS Anne-Laure
STETSUN Kateryna
Master 2 SSD - Université de Montpellier  
Année universitaire 2024-2025

## Références

- Cours d'Analyse Multivariée - X. Bry
- Projet CoForTips - Données du bassin du fleuve Congo
- Documentation R : [CRAN](https://cran.r-project.org/)

