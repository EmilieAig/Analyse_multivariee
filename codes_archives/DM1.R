library(dplyr)
library(ade4)
library(pls)
library(ggplot2)
library(FactoMineR)
library(factoextra)


# Lire les données avec tab comme séparateur
df1 <- read.csv("Datagenus.csv", sep = "\t", dec = ",", header = TRUE, stringsAsFactors = FALSE)

# Transform all character columns into numerical values
df1 <- df1 %>% mutate_if(is.character,~ as.numeric(.))

# Enlever la colonne forest
df1 <- df1 %>% select(-forest)

# 1.a)

df1$sum_sp <- rowSums(df1[, paste0("gen", 1:27)], na.rm = TRUE)
df1$treedensity <- df1$sum_sp / df1$surface

head(df1$treedensity)

# 1.b)

# geo
geo_vars <- c("lat", "lon", "altitude", "pluvio_yr", paste0("pluvio_", 1:12))

# carre
df1[paste0(geo_vars, "_2")] <- df1[, geo_vars]^2

# indicators for geology
df1$geology <- as.factor(df1$geology)
geology_dummies <- acm.disjonctif(df1["geology"])

evi_vars <- paste0("evi_", 1:23)

# Объединяем все переменные объяснения
X <- cbind(
  df1[, geo_vars],                       # исходные географические
  df1[, paste0(geo_vars, "_2")],         # квадраты географических
  geology_dummies,                        # dummy-переменные geology
  df1[, evi_vars]                        # EVI индексы
)

#var_explicative (dans Emilie code) = my X 
head(X)

# 2.a)

acp <- prcomp(X, scale. = TRUE, center = TRUE)

# Variance expliquée
var_expl <- acp$sdev^2 / sum(acp$sdev^2) * 100
var_expl_cum <- cumsum(var_expl)

fviz_eig(acp, addlabels = TRUE, ylim = c(0, max(var_expl[1:10])*1.2)) +
  ggtitle("Variance expliquée par composante") +
  ylab("% de variance expliquée") +
  xlab("Composantes principales")

# Просмотр накопленной дисперсии
print(data.frame(
  PC = 1:10,
  Var_pct = round(var_expl[1:10], 2),
  Var_cum = round(var_expl_cum[1:10], 2)
))

# Сохраняем первые 5 компонент
scores <- as.data.frame(acp$x[, 1:5])
colnames(scores) <- paste0("PC", 1:5)

# 2.b)

# Modèle linéaire : densité en fonction des composantes principales

# Добавляем зависимую переменную в таблицу компонент
scores$treedensity <- df1$treedensity

model_full  <- lm(treedensity ~ ., data = scores)
#same result as with
#model_full  <- lm(df1$treedensity ~ PC1 + PC2 + PC3 + PC4 + PC5, data = scores)
summary(model_full)

summary(model_full)$coefficients

model <- model_full
summary(model)

# R² models
R2 <- summary(model)$r.squared
cat("R² du modèle =", round(R2, 3), "\n")

# Graphique Y observé vs Y prédit
pred <- predict(model)
plot(scores$treedensity, pred,
     main = "Y observé vs Y prédit",
     xlab = "Y observé",
     ylab = "Y prédit",
     pch = 19, col = "blue",
     xlim = range(scores$treedensity),
     ylim = range(pred))
# Диагональ (идеальная линия)
abline(0, 1, col = "red", lwd = 2)
# Линия регрессии
abline(lm(pred ~ scores$treedensity), col = "darkgreen", lty = 2, lwd = 2)
# Подпись с R²
legend("topleft",
       legend = paste("R² =", round(R2, 3)),
       bty = "n", text.col = "black", cex = 1.1)


#graph like Emilie
plot(scores$treedensity, pred, 
     xlab = "Densité observée (Y)", 
     ylab = "Densité prédite (Ŷ)",
     main = "Valeurs observées vs prédites",
     pch = 20,
     col = rgb(0, 0, 1, 0.5))
# Ajouter la droite y = x (droite identité)
abline(a = 0, b = 1, col = "red", lwd = 2)
# Ajouter une droite de régression pour voir l'ajustement
abline(lm(pred ~ scores$treedensity), col = "blue", lty = 2, lwd = 2)
legend("topright", 
       legend = c("Droite identité (y=x)", "Droite de régression"),
       col = c("red", "blue"), 
       lty = c(1, 2), 
       lwd = 2,
       cex = 0.8)

# 2.c)
# Коэффициенты модели на главных компонентах (без константы)
beta_PC <- coef(model)[-1]

# Матрица загрузок (rotation) из PCA
rotation_matrix <- acp$rotation[, 1:5]

# Восстанавливаем коэффициенты оригинальных переменных (в стандартизированной форме)
beta_std <- rotation_matrix %*% beta_PC

# Коррекция стандартизации (возвращаемся к масштабу исходных переменных) - не надо
#sd_X <- apply(X, 2, sd)
#beta_orig <- as.vector(beta_std / sd_X)
beta_orig <- as.vector(beta_std)

# Создаём таблицу с именами переменных
beta_table <- data.frame(
  Variable = colnames(X),
  Coefficient = beta_orig
)

# Сортировка по вкладу
beta_table <- beta_table[order(abs(beta_table$Coefficient), decreasing = TRUE), ]

# Просмотр 10 самых влиятельных переменных
head(beta_table, 10)

# 2.d) Transformation log(Y + 1)
# Убираем Y, чтобы не попасть в ловушку
scores_noY <- scores[, !names(scores) %in% "treedensity"]
# Добавляем логарифмированную зависимую переменную
scores_noY$logY <- log1p(df1$treedensity)

# Строим модель
model_log <- lm(logY ~ ., data = scores_noY)
summary(model_log)

# comparison R²
R2_lin <- summary(model)$r.squared
R2_log <- summary(model_log)$r.squared
cat("R² linéaire :", round(R2_lin, 5), "\n")
cat("R² log(Y+1) :", round(R2_log, 5), "\n")

pred_log <- predict(model_log)
plot(df1$logY, pred_log,
     main = paste0("log(Y+1)"),
     xlab = "log(Y+1) observé", ylab = "log(Y+1) prédit",
     pch = 20, col = rgb(0,0.4,0,0.5))
abline(0, 1, col = "red", lwd = 2)
abline(lm(pred_log ~ df1$logY), col = "darkgreen", lty = 2, lwd = 2)
legend("topleft", 
       legend = c("Droite identité (y=y^)", "Droite de régression"),
       col = c("red", "blue"), 
       lty = c(1, 2), 
       lwd = 2,
       cex = 0.8)

head(X)
head(df1$treedensity)

# 3. Partition des variables explicatives

# Variables EVI
evi_vars <- paste0("evi_", 1:23)

# Partie Photosynthèse (EVI)
X_evi <- X[, evi_vars]
# Partie Géographie (tout le reste)
X_geo <- X[, !(colnames(X) %in% evi_vars)]

# Проверим размеры
dim(X_geo)
dim(X_evi)

cat("Variables géographiques :", ncol(X_geo), "\n")
cat("Variables de photosynthèse (EVI) :", ncol(X_evi), "\n")

# 3.a)
# ---- ACP pour Géographie ----
acp_geo <- prcomp(X_geo, scale. = TRUE, center = TRUE)

# Variance expliquée
var_expl_geo <- acp_geo$sdev^2 / sum(acp_geo$sdev^2) * 100
var_expl_geo_cum <- cumsum(var_expl_geo)

# first 10
fviz_eig(acp_geo, addlabels = TRUE, ylim = c(0, max(var_expl_geo[1:10])*1.2)) +
  ggtitle("Variance expliquée par composantes géographiques") +
  ylab("% de variance expliquée") +
  xlab("Composantes principales")

# disters
print(data.frame(
  PC = 1:10,
  Var_pct = round(var_expl_geo[1:10], 3),
  Var_cum = round(var_expl_geo_cum[1:10], 3)
))

# Сохраняем первые n_компонент (например, 5)
scores_geo <- as.data.frame(acp_geo$x[, 1:3])
colnames(scores_geo) <- paste0("Geo_PC", 1:3)

# ----  ACP pour EVI ----
acp_evi <- prcomp(X_evi, scale. = TRUE, center = TRUE)

var_expl_evi <- acp_evi$sdev^2 / sum(acp_evi$sdev^2) * 100
var_expl_evi_cum <- cumsum(var_expl_evi)

fviz_eig(acp_evi, addlabels = TRUE, ylim = c(0, max(var_expl_evi[1:10])*1.2)) +
  ggtitle("Variance expliquée par composantes EVI") +
  ylab("% de variance expliquée") +
  xlab("Composantes principales")

print(data.frame(
  PC = 1:10,
  Var_pct = round(var_expl_evi[1:10], 3),
  Var_cum = round(var_expl_evi_cum[1:10], 3)
))

# Сохраняем первые m_компонент (например, 5)
scores_evi <- as.data.frame(acp_evi$x[, 1:5])
colnames(scores_evi) <- paste0("EVI_PC", 1:5)

# ---- Объединяем компоненты для последующей регрессии ----
scores_combined <- cbind(scores_geo, scores_evi)

head(scores_combined)

# 3.b) -----тут началась какая-то хуйня------------

# Полная модель на выбранных компонентах
model_combined <- lm(treedensity ~ ., data = scores_combined)
summary(model_combined)

# Смотрим, какие компоненты значимы (p < 0.05)
summary(model_combined)$coefficients

# Объединяем компоненты
scores_combined <- cbind(scores_geo_sel, scores_evi_sel)

# Добавляем зависимую переменную
scores_combined$logY <- log1p(df1$treedensity)  # если мы используем трансформацию из 2.d

# Предсказанные значения
pred_combined <- predict(model_combined)

# R²
R2_combined <- summary(model_combined)$r.squared
cat("R² du modèle combiné =", round(R2_combined, 3), "\n")

# График
plot(scores_combined$logY, pred_combined,
     main = paste0("Y observé vs Y prédit (ACP combinée)\nR² = ", round(R2_combined,3)),
     xlab = "log(Y+1) observé",
     ylab = "log(Y+1) prédit",
     pch = 19, col = rgb(0,0.4,0,0.5))
abline(0, 1, col = "red", lwd = 2)                      # линия y=x
abline(lm(pred_combined ~ scores_combined$logY), 
       col = "darkgreen", lty = 2, lwd = 2)            # линия регрессии
legend("topleft", 
       legend = c("Droite identité (y=y^)", "Droite de régression"),
       col = c("red", "darkgreen"), lty = c(1,2), lwd = 2, cex = 0.8)

# 3.c) Projection sur ces composantes principales

#X_geo_acp <- acp_geo$x[, 1:nb_geo]
#X_evi_acp <- acp_evi$x[, 1:nb_evi]

X_geo_acp <- acp_geo$x[, 1:n_comp_geo]
X_evi_acp <- acp_evi$x[, 1:n_comp_evi]

# Réunion des composantes principales
X_final <- cbind(X_geo_acp, X_evi_acp)

# Transformation de Y : log(Y + 1)
Y_log <- log(df1$treedensity + 1)

# 3.d) Régression linéaire sur les composantes principales

#modele2 <- lm(df1$treedensity ~ X_final)

modele2 <- lm(Y_log ~ X_final)
summary(modele2)

# Vérification des résidus
plot(modele2$residuals, main = "Résidus du modèle ACP séparé")
hist(modele2$residuals, main = "Histogramme des résidus")

# 3.e) Vérification de la corrélation entre les composantes principales
#cor_matrix <- cor(X_final)
#round(cor_matrix, 2)

# Afficher sous forme de heatmap pour visualiser les corrélations
#corrplot::corrplot(cor_matrix, method = "color", tl.cex = 0.7, number.cex = 0.6)

# === 3.e) Graphique Y observé vs Y prédit ===

# === 3.f) Vérification des corrélations entre composantes principales ===


# 4. Régression PLS:

# а
