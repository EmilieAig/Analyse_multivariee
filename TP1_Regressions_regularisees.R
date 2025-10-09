################################################################################
## TP1 Régressions régularisées
################################################################################

# Charger les packages nécessaires
library(dplyr)
library(ade4)
library(corrplot)
library(caret)
library(car)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(pls)
library(glmnet)

# Charger les données
data <- read.csv("Datagenus.csv", sep = '')

# Explorer les données
head(data)
names(data)
str(data)

################################################################################

# Question 1

data$surface <- gsub(",", ".", data$surface)
data$surface <- as.numeric(data$surface)

# 1.a Calcul de la variable dépendante : treedensity
# Sélection des colonnes gen1 à gen27
abundances <- data[, grep("^gen[0-9]+$", names(data))]

# Calcul de la densité de peuplement
data$treedensity <- rowSums(abundances) / data$surface

summary(data$treedensity)

# Visualisation avec un histogramme
hist(data$treedensity,
     main = "Distribution de la densité de peuplement",
     xlab = "Densité des arbres",
     ylab = "Fréquence",
     col = "#CC79A7")

#1.b Construction du tableau des variables explicatives

# Variables géographiques quantitatives
geo_vars <- c("lat", "lon", "altitude", "pluvio_yr", paste0("pluvio_", 1:12))
geo_data <- data[, geo_vars]
geo_data <- data.frame(lapply(geo_data, function(x) gsub(",", ".", x)))
geo_data <- data.frame(lapply(geo_data, function(x) as.numeric(as.character(x))))

# Carrés des variables géographiques
geo_squared <- geo_data^2
colnames(geo_squared) <- paste0(geo_vars, "_sq")

# Indicateurs de la variable qualitative geology
geology_dummies <- acm.disjonctif(data["geology"])

# Variables EVI (supposons qu'elles commencent par "EVI")
evi_vars <- grep("^evi_", names(data), value = TRUE)
evi_data <- data[, evi_vars]
evi_data <- data.frame(lapply(evi_data, function(x) as.numeric(gsub(",", ".", x))))

# Construction du tableau explicatif
X <- cbind(geo_data, geo_squared, geology_dummies, evi_data)

# Vérification
dim(X)
summary(X)
str(X)

# Inventaire théorique des multi-colinéarités

# Matrice de corrélation
cor_matrix <- cor(X, use = "pairwise.complete.obs")

# Visualisation graphique
corrplot(cor_matrix, method = "color", tl.cex = 0.6,
         title = "Corrélation entre variables explicatives")

# Détection des corrélations fortes (> 0.9)
high_corr <- which(abs(cor_matrix) > 0.9 & abs(cor_matrix) < 1, arr.ind = TRUE)
high_corr_vars <- unique(c(rownames(high_corr), colnames(high_corr)))

cat("Variables fortement corrélées (corr > 0.9) :\n")
print(high_corr_vars)

# Analyse théorique :
# 1. Les carrés des variables géographiques sont naturellement corrélés aux variables originales.
# 2. Les indicatrices de 'geology' sont redondantes (leur somme = 1). On supprime la première pour éviter l'alias.
# 3. Certaines EVI peuvent être très corrélées entre elles. On pourrait éventuellement garder une seule variable par groupe fortement corrélé.

# Suppression de la première colonne d'indicatrices de 'geology' pour éviter la redondance
geology_dummies <- geology_dummies[, -1]

# Reconstruction du tableau explicatif sans redondance
X <- cbind(geo_data, geo_squared, geology_dummies, evi_data)

# Vérification finale
dim(X)
summary(X)

################################################################################

# Question 2 : Régression PLS

# Standardisation des variables explicatives
X_scaled <- scale(X)

# 2.a — ACP réduite globale
acp <- PCA(X_scaled, graph = FALSE)

# Visualisation de la variance expliquée
fviz_screeplot(acp,
               addlabels = TRUE,
               ylim = c(0, 50),
               main ="Variance expliquée",
               xlab = "Dimensions",
               ylab = "Pourcentages de variance expliquée")

# Sélection des composantes non bruitées
# Critère de Kaiser : garder les composantes avec une valeur propre > 1
# selected_dims <- which(eig_vals[, "eigenvalue"] > 1)
# Critère du seuil cumulatif : garder les composantes jusqu’à 80% de variance cumulée
# selected_dims <- which(eig_vals[, "cumulative percentage of variance"] <= 80) 
# A savoir que je peux combiner les deux si il le faut : selected_dims <- which(eig_vals[, "eigenvalue"] > 1 & eig_vals[, "cumulative percentage of variance"] <= 80)
eig_vals <- acp$eig
eig_vals
selected_dims <- which(eig_vals[, "cumulative percentage of variance"] <= 80)

cat("Composantes retenues :\n")
print(selected_dims)

# Extraction des composantes principales retenues
PCs <- acp$ind$coord[, selected_dims]

# 2.b — Régression linéaire sur les composantes retenues
pc_df <- data.frame(PCs)
pc_df$treedensity <- data$treedensity

model_pc <- lm(treedensity ~ ., data = pc_df)
summary(model_pc)

# Élimination des composantes non significatives (p > 0.05)
significant <- summary(model_pc)$coefficients[-1, 4] < 0.05
selected_PC_names <- names(significant)[significant]

# Nouveau modèle avec composantes significatives
model_pc_reduced <- lm(treedensity ~ ., data = pc_df[, c(selected_PC_names, "treedensity")])
summary(model_pc_reduced)

# R² du modèle réduit
cat("R² du modèle réduit :", summary(model_pc_reduced)$r.squared, "\n")

# Graphe Y vs Ŷ avec courbe de régression 
Y_hat <- predict(model_pc_reduced)

# Couleur violette semi-transparente (alpha = 0.5)
point_col <- rgb(0, 0, 1, alpha = 0.5) # Plus l’alpha est bas, plus les points cumulés deviennent foncés

# Crée le graphique
plot(data$treedensity, Y_hat,
     xlab = "Densité observée (Y)",
     ylab = "Densité prédite (Ŷ)",
     main = "Régression sur composantes principales",
     col = point_col,
     pch = 20)

# Ajouter la droite y = x (droite identité) : rouge
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

# Droite de régression empirique : bleue pointillée
reg_line <- lm(Y_hat ~ data$treedensity)
abline(reg_line, col = "blue", lty = 3, lwd = 2)

# Récupérer les limites du graphique
xlim <- par("usr")[1:2]
ylim <- par("usr")[3:4]

# Légende compacte, collée au bord droit
legend(x = xlim[2] * 0.995,     # très proche du bord droit
       y = ylim[2],             # tout en haut du graphique
       legend = c("Points observés", "Droite identité (Y=Ŷ)", "Droite de régression"),
       col = c(point_col, "red", "blue"),
       pch = c(16, NA, NA),       
       lty = c(NA, 2, 3),         
       lwd = c(NA, 2, 2),
       bty = "o",
       box.lwd = 0.8,
       bg = rgb(1,1,1,0.8),
       cex = 0.8,
       xjust = 1,
       yjust = 1,
       x.intersp = 0.6,
       y.intersp = 0.8,
       seg.len = 1)

# Calculer la corrélation pour commentaire
correlation_YYhat <- cor(data$treedensity, Y_hat)
cat("Corrélation entre Y et Ŷ :", round(correlation_YYhat, 4), "\n")

# 2.c — Retrouver les coefficients des variables originelles
# Coefficients = combinaison linéaire des loadings et des coefficients du modèle
loadings <- acp$var$coord[, selected_dims]
beta_pc <- coef(model_pc)[-1]
beta_orig <- loadings %*% beta_pc
beta_orig <- as.vector(beta_orig)
names(beta_orig) <- colnames(X_scaled)

cat("Coefficients des variables originelles dans Ŷ :\n")
print(round(beta_orig, 4))

# 2.d — Correction de la linéarité avec transformation log(Y)
# Transformation log1p (log(1 + Y)) pour préserver les zéros
data$log_treedensity <- log1p(data$treedensity)

# Nouveau modèle avec log(Y)
pc_df$log_treedensity <- data$log_treedensity
model_log <- lm(log_treedensity ~ ., data = pc_df[, c(selected_PC_names, "log_treedensity")])
summary(model_log)

# Graphe log(Y) vs Ŷ
Y_hat_log <- predict(model_log)

# Couleur verte semi-transparente pour les points
point_col_log <- rgb(0, 0.6, 0, 0.5)

plot(data$log_treedensity, Y_hat_log,
     xlab = "log(Densité + 1) observé", 
     ylab = "log(Densité + 1) prédit",
     main = "Régression avec transformation log",
     pch = 20,
     col = point_col_log,
     cex.main = 0.9)

# Droite identité (Y = Ŷ) : rouge
abline(a = 0, b = 1, col = "red",lty = 3, lwd = 2)

# Droite de régression empirique : vert foncé pointillée
reg_log_line <- lm(Y_hat_log ~ data$log_treedensity)
abline(reg_log_line, col = rgb(0, 0.3, 0, 0.8), lty = 2, lwd = 2)

# Limites du graphique
xlim <- par("usr")[1:2]
ylim <- par("usr")[3:4]

# Légende compacte, collée au bord droit
legend(x = xlim[1] + 0.02 * diff(xlim), y = ylim[2],
       legend = c("Points observés", "Droite identité (Y=Ŷ)", "Droite de régression"),
       col = c(point_col_log, "red", rgb(0, 0.3, 0, 0.8)),
       pch = c(20, NA, NA),
       lty = c(NA, 3, 2),
       lwd = c(NA, 2, 2),
       bty = "o",
       box.lwd = 0.8,
       bg = rgb(1,1,1,0.8),
       cex = 0.8,
       xjust = 0,
       yjust = 1,
       x.intersp = 0.6,
       y.intersp = 0.8,
       seg.len = 1)

# Calculer les corrélations pour comparer
cor_original <- cor(data$treedensity, Y_hat)
cor_log <- cor(data$log_treedensity, Y_hat_log)

cat("\n=== CORRÉLATIONS Y vs Ŷ ===\n")
cat("Modèle original :", round(cor_original, 4), "\n")
cat("Modèle log      :", round(cor_log, 4), "\n")

beta_orig_globale <- beta_orig  # on garde les coefficients de la RCP globale

################################################################################

# Question 3 : Régression sur ACP séparées (Photosynthèse & Géographie)

# Standardisation des deux blocs
geo_scaled <- scale(cbind(geo_data, geo_squared, geology_dummies))
evi_scaled <- scale(evi_data)

# 3.a — ACP séparée sur le bloc Géographie
acp_geo <- PCA(geo_scaled, graph = FALSE)
fviz_screeplot(acp_geo,
               addlabels = TRUE,
               main = "Variance expliquée",
               xlab = "Dimensions",
               ylab = "Pourcentages de variance expliquée")

# Sélection des composantes géographiques (valeurs propres > 1)
eig_geo <- acp_geo$eig
eig_geo

# A voir si on peut pas faire avec cumul parce que error
geo_dims <- which(eig_geo[, "cumulative percentage of variance"] <= 80)
geo_dims
geo_PCs <- acp_geo$ind$coord[, geo_dims]

# ACP séparée sur le bloc Photosynthèse (EVI)
acp_evi <- PCA(evi_scaled, graph = FALSE)
fviz_screeplot(acp_evi, addlabels = TRUE, main = "ACP Photosynthèse")

# Sélection des composantes EVI (valeurs propres > 1)
eig_evi <- acp_evi$eig
evi_dims <- which(eig_evi[, "eigenvalue"] > 1)
evi_PCs <- acp_evi$ind$coord[, evi_dims]

# Fusion des composantes retenues
PC_all <- cbind(geo_PCs, evi_PCs)
pc_df <- data.frame(PC_all)
pc_df$treedensity <- data$treedensity

# 3.b — Régression sur toutes les composantes
model_all <- lm(treedensity ~ ., data = pc_df)
summary(model_all)

# Sélection des composantes significatives (p < 0.05)
significant <- summary(model_all)$coefficients[-1, 4] < 0.05
selected_PC_names <- names(significant)[significant]

# Modèle réduit
model_reduced <- lm(treedensity ~ ., data = pc_df[, c(selected_PC_names, "treedensity")])
summary(model_reduced)

# R² du modèle réduit
cat("R² du modèle réduit :", summary(model_reduced)$r.squared, "\n")

# Graphe Y vs Ŷ
Y_hat <- predict(model_reduced)
plot(data$treedensity, Y_hat,
     xlab = "Y observé (treedensity)",
     ylab = "Y prédit (Ŷ)",
     main = "Régression sur ACP thématique",
     col = "#0072B2", pch = 16)
abline(0, 1, lty = 2)

# 3.c — Retrouver les coefficients des variables originelles
# Coefficients = combinaison linéaire des loadings et des coefficients du modèle
load_geo <- acp_geo$var$coord[, geo_dims]
load_evi <- acp_evi$var$coord[, evi_dims]

beta_geo <- coef(model_all)[1 + seq_along(geo_dims)]
beta_evi <- coef(model_all)[(1 + length(geo_dims)):(length(geo_dims) + length(evi_dims))]

beta_orig_geo <- load_geo %*% beta_geo
beta_orig_evi <- load_evi %*% beta_evi

beta_orig <- c(beta_orig_geo, beta_orig_evi)
names(beta_orig) <- c(colnames(geo_scaled), colnames(evi_scaled))

cat("Coefficients des variables originelles dans Ŷ :\n")
print(round(beta_orig, 4))

# 3.d — Transformation log(Y) pour corriger la linéarité
data$log_treedensity <- log1p(data$treedensity)
pc_df$log_treedensity <- data$log_treedensity

model_log <- lm(log_treedensity ~ ., data = pc_df[, c(selected_PC_names, "log_treedensity")])
summary(model_log)

# Graphe log(Y) vs Ŷ
Y_hat_log <- predict(model_log)
plot(data$log_treedensity, Y_hat_log,
     xlab = "log(1 + Y) observé",
     ylab = "log(1 + Ŷ) prédit",
     main = "Régression log-transformée sur ACP thématique",
     col = "#CC79A7", pch = 16)
abline(0, 1, lty = 2)

beta_orig_themes <- beta_orig  # on garde les coefficients de la RCP thématique

################################################################################

# Question 4 : Régression PLS

# Préparation des données : transformation log si nécessaire
Y_pls <- data$log_treedensity
X_pls <- scale(X)  # variables explicatives standardisées

# 4.a — Régression PLS avec validation croisée (K-fold, k=10)
set.seed(123)
pls_model <- plsr(Y_pls ~ X_pls, ncomp = 20, validation = "CV") 

# Visualiser l'erreur de prédiction (PRESS) pour choisir le nombre optimal de composantes
validationplot(pls_model, val.type = "MSEP", main="Erreur de validation croisée (MSEP)")

# Nombre optimal de composantes
optimal_ncomp <- which.min(pls_model$validation$PRESS)
cat("Nombre optimal de composantes PLS :", optimal_ncomp, "\n")

# 4.b — Interprétation des composantes retenues
# Scores et loadings
pls_scores <- pls_model$scores[, 1:optimal_ncomp]
pls_loadings <- pls_model$loadings[, 1:optimal_ncomp]

# 4.c — Retrouver les coefficients des variables originales
pls_coef <- coef(pls_model, ncomp = optimal_ncomp)
pls_coef <- as.vector(pls_coef)
names(pls_coef) <- colnames(X_pls)
cat("Coefficients PLS des variables originales :\n")
print(round(pls_coef, 4))

# 4.d — Prédiction et visualisation
Y_hat_pls <- predict(pls_model, ncomp = optimal_ncomp)
plot(Y_pls, Y_hat_pls,
     xlab = "log(1 + Y) observé",
     ylab = "log(1 + Ŷ) prédit (PLS)",
     main = "Régression PLS",
     col = "#D55E00", pch = 16)
abline(0, 1, lty = 2)

################################################################################

# Question 5 : Régressions pénalisées (Ridge et LASSO)

# Préparation des données : standardisation et transformation log si nécessaire
Y_glm <- data$log_treedensity
X_glm <- scale(X)  # variables explicatives standardisées

# Conversion en matrice (glmnet nécessite des matrices)
X_matrix <- as.matrix(X_glm)

# -----------------------------
# 5.a — Régression Ridge
# -----------------------------
set.seed(123)
# alpha = 0 pour Ridge
ridge_cv <- cv.glmnet(X_matrix, Y_glm, alpha = 0, nfolds = 10)
plot(ridge_cv)
cat("Lambda optimal Ridge :", ridge_cv$lambda.min, "\n")

# Modèle Ridge final
ridge_model <- glmnet(X_matrix, Y_glm, alpha = 0, lambda = ridge_cv$lambda.min)
ridge_coef <- as.vector(coef(ridge_model))[-1]  # enlever l'intercept
names(ridge_coef) <- colnames(X_matrix)
cat("Coefficients Ridge :\n")
print(round(ridge_coef, 4))

# Prédiction et visualisation
Y_hat_ridge <- predict(ridge_model, X_matrix)
plot(Y_glm, Y_hat_ridge,
     xlab = "log(1 + Y) observé",
     ylab = "log(1 + Ŷ) prédit (Ridge)",
     main = "Ridge Regression",
     col = "#009E73", pch = 16)
abline(0, 1, lty = 2)

# -----------------------------
# 5.b — Régression LASSO
# -----------------------------
set.seed(123)
# alpha = 1 pour LASSO
lasso_cv <- cv.glmnet(X_matrix, Y_glm, alpha = 1, nfolds = 10)
plot(lasso_cv)
cat("Lambda optimal LASSO :", lasso_cv$lambda.min, "\n")

# Modèle LASSO final
lasso_model <- glmnet(X_matrix, Y_glm, alpha = 1, lambda = lasso_cv$lambda.min)
lasso_coef <- as.vector(coef(lasso_model))[-1]  # enlever l'intercept
names(lasso_coef) <- colnames(X_matrix)
cat("Coefficients LASSO :\n")
print(round(lasso_coef, 4))

# Prédiction et visualisation
Y_hat_lasso <- predict(lasso_model, X_matrix)
plot(Y_glm, Y_hat_lasso,
     xlab = "log(1 + Y) observé",
     ylab = "log(1 + Ŷ) prédit (LASSO)",
     main = "LASSO Regression",
     col = "#E69F00", pch = 16)
abline(0, 1, lty = 2)

################################################################################

# Question 6 : Synthèse et Tableau comparatif des coefficients

# Harmonisation des longueurs
all_vars <- colnames(X)

# Compléter les vecteurs manquants avec des zéros si besoin
pad <- function(v) {
  out <- rep(0, length(all_vars))
  names(out) <- all_vars
  common <- intersect(names(v), all_vars)
  out[common] <- v[common]
  return(out)
}

coeff_matrix <- data.frame(
  Variable = all_vars,
  RCP_globale = pad(beta_orig_globale),
  RCP_themes  = pad(beta_orig_themes),
  PLS         = pad(pls_coef),
  Ridge       = pad(ridge_coef),
  LASSO       = pad(lasso_coef)
)

# Ajouter R² pour chaque méthode
R2_values <- data.frame(
  Methode = c("RCP_globale", "RCP_themes", "PLS", "Ridge", "LASSO"),
  R2 = c(
    summary(model_pc_reduced)$r.squared,
    summary(model_reduced)$r.squared,
    cor(Y_pls, Y_hat_pls)^2,   # PLS
    cor(Y_glm, Y_hat_ridge)^2, # Ridge
    cor(Y_glm, Y_hat_lasso)^2  # LASSO
  )
)

# Affichage
print(R2_values)

# Visualisation des coefficients
library(reshape2)
coeff_long <- melt(coeff_matrix, id.vars = "Variable",
                   variable.name = "Methode", value.name = "Coefficient")

ggplot(coeff_long, aes(x = Variable, y = Coefficient, fill = Methode)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Comparaison des coefficients selon la méthode",
       x = "Variables explicatives",
       y = "Coefficient") +
  scale_fill_brewer(palette = "Set2")

