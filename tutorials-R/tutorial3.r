# TP d'Introduction à l'Économétrie
# Modèle Linéaire Simple
#
# Ce TP reprend les concepts fondamentaux du modèle linéaire simple présentés dans le cours.
# Nous allons explorer les méthodes d'estimation (MCO et maximum de vraisemblance), 
# les hypothèses du modèle, les tests statistiques et l'interprétation des résultats.

## 1. Importation des bibliothèques nécessaires
library(ggplot2)
library(dplyr)
library(tidyr)
library(car)
library(lmtest)
library(stats)

# Configuration pour reproduire les résultats
set.seed(0)

## 2. Génération de données simulées
# Nous allons simuler un modèle linéaire simple de la forme : Y_i = α + β X_i + u_i
#
# Où :
# - α est la constante (l'ordonnée à l'origine)
# - β est le coefficient de la variable explicative
# - u_i est le terme d'erreur qui suit une loi normale

# Paramètres du modèle
alpha_vrai <- 10  # Vraie valeur de alpha
beta_vrai <- 2    # Vraie valeur de beta
sigma_u <- 3      # Écart-type du terme d'erreur
n <- 100          # Nombre d'observations

# Génération de la variable explicative X
X <- runif(n, 0, 10)

# Génération du terme d'erreur u (suivant une loi normale centrée de variance sigma_u²)
u <- rnorm(n, 0, sigma_u)

# Génération de la variable dépendante Y selon le modèle
Y <- alpha_vrai + beta_vrai * X + u

# Création d'un data frame pour stocker les données
df <- data.frame(
  X = X,
  Y = Y,
  u = u  # On garde les erreurs pour vérification ultérieure
)

# Affichage des premières lignes du data frame
head(df)

## 3. Visualisation des données
ggplot(df, aes(x = X, y = Y)) +
  geom_point(alpha = 0.7) +
  geom_abline(intercept = alpha_vrai, slope = beta_vrai, color = "red", 
              linetype = "solid", size = 1) +
  labs(title = "Relation entre X et Y avec le vrai modèle",
       subtitle = paste("Y =", alpha_vrai, "+", beta_vrai, "X")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 14))

## 4. Estimation par la méthode des Moindres Carrés Ordinaires (MCO)
#
# Rappel des formules des estimateurs MCO :
#
# β̂ = Σ(X_i - X̄)(Y_i - Ȳ) / Σ(X_i - X̄)²
#
# α̂ = Ȳ - β̂X̄
#
# Nous allons implémenter ces formules manuellement puis comparer avec l'estimateur de R.

# Calcul manuel des estimateurs MCO
X_mean <- mean(X)
Y_mean <- mean(Y)

# Calcul de beta
numerator <- sum((X - X_mean) * (Y - Y_mean))
denominator <- sum((X - X_mean)^2)
beta_hat <- numerator / denominator

# Calcul de alpha
alpha_hat <- Y_mean - beta_hat * X_mean

# Affichage des résultats
cat(sprintf("Estimateur MCO manuel pour α: %.4f\n", alpha_hat))
cat(sprintf("Estimateur MCO manuel pour β: %.4f\n", beta_hat))
cat(sprintf("Valeurs réelles: α = %s, β = %s\n", alpha_vrai, beta_vrai))

# Utilisation de la fonction lm de R pour l'estimation MCO
model <- lm(Y ~ X, data = df)
summary(model)

## 5. Calcul et interprétation des résidus

# Calcul des valeurs prédites
Y_pred <- alpha_hat + beta_hat * X

# Calcul des résidus
residus <- Y - Y_pred

# Ajout des résidus et valeurs prédites au data frame
df$Y_pred <- Y_pred
df$residus <- residus

# Estimation de la variance des erreurs
sigma_u_hat <- sqrt(sum(residus^2) / (n - 2))

cat(sprintf("Estimation de l'écart-type des erreurs: %.4f\n", sigma_u_hat))
cat(sprintf("Vraie valeur de l'écart-type des erreurs: %.4f\n", sigma_u))

# Visualisation des résidus
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Premier graphique: résidus en fonction de X
plot(X, residus, main = "Résidus en fonction de X", 
     xlab = "X", ylab = "Résidus", pch = 19, col = "blue", alpha = 0.7)
abline(h = 0, col = "red", lwd = 2)
grid()

# Deuxième graphique: histogramme des résidus
hist(residus, main = "Distribution des résidus", 
     xlab = "Résidus", col = "lightblue", border = "white", breaks = 15)
lines(density(residus), col = "red", lwd = 2)

# Alternative avec ggplot2
p1 <- ggplot(df, aes(x = X, y = residus)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red", size = 1) +
  labs(title = "Résidus en fonction de X") +
  theme_minimal()

p2 <- ggplot(df, aes(x = residus)) +
  geom_histogram(aes(y = ..density..), fill = "lightblue", 
                bins = 15, alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  labs(title = "Distribution des résidus") +
  theme_minimal()

# Pour afficher ces graphiques côte à côte, utilisez la bibliothèque gridExtra
# library(gridExtra)
# grid.arrange(p1, p2, ncol = 2)

## 6. Vérification des hypothèses du modèle linéaire simple
#
# Rappel des hypothèses :
# 1. La distribution de l'erreur u est indépendante de X
# 2. L'erreur est centrée et de variance constante (homoscédasticité)
# 3. Les coefficients α et β sont constants
# 4. Les résidus suivent une loi normale
# 5. Non-autocorrélation des résidus

# Diagnostics du modèle
par(mfrow = c(2, 2))
plot(model)

# Test de normalité des résidus (Shapiro-Wilk)
shapiro_test <- shapiro.test(residus)
cat("Test de Shapiro-Wilk pour la normalité des résidus:\n")
cat(sprintf("Statistique W: %.4f\n", shapiro_test$statistic))
cat(sprintf("p-value: %.4f\n", shapiro_test$p.value))
cat(sprintf("Conclusion: %s\n", 
            ifelse(shapiro_test$p.value > 0.05, "Résidus normaux", "Résidus non normaux")))

# Test d'homoscédasticité (Breusch-Pagan)
bp_test <- bptest(model)
cat("\nTest de Breusch-Pagan pour l'homoscédasticité:\n")
cat(sprintf("Statistique BP: %.4f\n", bp_test$statistic))
cat(sprintf("p-value: %.4f\n", bp_test$p.value))
cat(sprintf("Conclusion: %s\n", 
            ifelse(bp_test$p.value > 0.05, "Homoscédasticité", "Hétéroscédasticité")))

## 7. Analyse de la variance et qualité de l'ajustement (R²)
#
# Rappel de la décomposition de la variance totale :
#
# Σ(Y_i - Ȳ)² = Σ(Ŷ_i - Ȳ)² + Σ(Y_i - Ŷ_i)²
#
# SCT = SCE + SCR
#
# Et le coefficient de détermination :
#
# R² = SCE/SCT = 1 - SCR/SCT

# Calcul de la somme des carrés totale (SCT)
SCT <- sum((Y - Y_mean)^2)

# Calcul de la somme des carrés expliquée (SCE)
SCE <- sum((Y_pred - Y_mean)^2)

# Calcul de la somme des carrés résiduelle (SCR)
SCR <- sum((Y - Y_pred)^2)

# Calcul du R²
R2 <- SCE / SCT
# Vérification: R² = 1 - SCR/SCT
R2_verif <- 1 - SCR / SCT

cat(sprintf("Somme des carrés totale (SCT): %.4f\n", SCT))
cat(sprintf("Somme des carrés expliquée (SCE): %.4f\n", SCE))
cat(sprintf("Somme des carrés résiduelle (SCR): %.4f\n", SCR))
cat(sprintf("Vérification: SCT = SCE + SCR -> %.4f = %.4f + %.4f -> %.4f\n", 
            SCT, SCE, SCR, SCE + SCR))
cat(sprintf("\nR² calculé manuellement: %.4f\n", R2))
cat(sprintf("R² calculé par vérification: %.4f\n", R2_verif))
cat(sprintf("R² fourni par la fonction lm: %.4f\n", summary(model)$r.squared))

# Visualisation de la décomposition de la variance
ggplot(df, aes(x = X, y = Y)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, color = "red", size = 1.5,
              se = FALSE) +
  geom_hline(yintercept = Y_mean, color = "green", linetype = "dashed",
            alpha = 0.5) +
  geom_vline(xintercept = X_mean, color = "black", linetype = "dashed",
            alpha = 0.5) +
  # Ajouter quelques lignes pour illustrer la décomposition de la variance
  geom_segment(data = df[seq(1, n, 10), ], 
              aes(x = X, y = Y, xend = X, yend = Y_pred),
              color = "blue", alpha = 0.3) +
  geom_segment(data = df[seq(1, n, 10), ], 
              aes(x = X, y = Y_pred, xend = X, yend = Y_mean),
              color = "green", alpha = 0.3) +
  labs(title = paste("Régression linéaire avec décomposition de la variance (R² =", 
                      round(R2, 4), ")"),
       subtitle = paste("Y =", round(alpha_hat, 4), "+", round(beta_hat, 4), "X")) +
  theme_minimal()

## 8. Tests statistiques sur les coefficients
#
# Nous allons tester les hypothèses : H0: β = 0 vs H1: β ≠ 0
#
# Rappel : sous les hypothèses du modèle linéaire simple :
#
# (β̂ - β) / sqrt(s² * C_β) ~ T(n-2)
#
# Où C_β = 1 / (n * Σ(X_i - X̄)²)

# Calcul de l'écart-type de beta
C_beta <- 1 / (n * sum((X - X_mean)^2))
se_beta <- sqrt(sigma_u_hat^2 * C_beta)

# Statistique de test pour H0: beta = 0
t_stat <- beta_hat / se_beta

# Valeur critique et p-value
dof <- n - 2  # degrés de liberté
p_value <- 2 * (1 - pt(abs(t_stat), dof))  # test bilatéral

# Intervalle de confiance à 95%
t_crit <- qt(0.975, dof)  # quantile à 97.5% (test bilatéral à 5%)
ci_lower <- beta_hat - t_crit * se_beta
ci_upper <- beta_hat + t_crit * se_beta

cat("Test de significativité de β:\n")
cat("H0: β = 0 vs H1: β ≠ 0\n")
cat(sprintf("Estimateur de β: %.4f\n", beta_hat))
cat(sprintf("Écart-type de β: %.4f\n", se_beta))
cat(sprintf("Statistique t: %.4f\n", t_stat))
cat(sprintf("p-value: %.6f\n", p_value))
cat(sprintf("Conclusion: %s au seuil de 5%%\n", 
            ifelse(p_value < 0.05, "Rejet de H0", "Non-rejet de H0")))
cat(sprintf("\nIntervalle de confiance à 95%% pour β: [%.4f, %.4f]\n", 
            ci_lower, ci_upper))
cat(sprintf("La vraie valeur de β (%s) est-elle dans l'intervalle? %s\n", 
            beta_vrai, ifelse(ci_lower <= beta_vrai & beta_vrai <= ci_upper, "Oui", "Non")))

## 9. Prévisions avec le modèle estimé
#
# Maintenant que nous avons estimé notre modèle, nous pouvons l'utiliser pour faire des prévisions.

# Nouvelle valeur de X pour la prévision
X_new <- 15  # Une valeur en dehors de la plage des données

# Prévision ponctuelle
Y_pred_new <- alpha_hat + beta_hat * X_new

# Intervalle de prévision à 95%
# Formule: Y_pred ± t_{n-2, 0.975} * s * sqrt(1 + 1/n + (X_new - X̄)²/Σ(X_i - X̄)²)
prediction_variance <- sigma_u_hat^2 * (1 + 1/n + ((X_new - X_mean)^2) / sum((X - X_mean)^2))
prediction_se <- sqrt(prediction_variance)
prediction_margin <- t_crit * prediction_se

pred_lower <- Y_pred_new - prediction_margin
pred_upper <- Y_pred_new + prediction_margin

cat(sprintf("Prévision pour X = %s:\n", X_new))
cat(sprintf("Valeur prédite: %.4f\n", Y_pred_new))
cat(sprintf("Intervalle de prévision à 95%%: [%.4f, %.4f]\n", 
            pred_lower, pred_upper))

# Visualisation avec predict() pour calculer l'intervalle de prévision
new_data <- data.frame(X = seq(0, 16, length.out = 100))
pred_interval <- predict(model, newdata = new_data, interval = "prediction", level = 0.95)
new_data$fit <- pred_interval[, 1]
new_data$lwr <- pred_interval[, 2]
new_data$upr <- pred_interval[, 3]

ggplot() +
  geom_point(data = df, aes(x = X, y = Y), alpha = 0.7) +
  geom_line(data = new_data, aes(x = X, y = fit), color = "red", size = 1.2) +
  geom_ribbon(data = new_data, aes(x = X, ymin = lwr, ymax = upr), 
              fill = "gray", alpha = 0.2) +
  geom_point(aes(x = X_new, y = Y_pred_new), color = "green", size = 4, shape = 4) +
  geom_errorbar(aes(x = X_new, ymin = pred_lower, ymax = pred_upper), 
                color = "green", width = 0.5, linetype = "dashed", alpha = 0.7) +
  labs(title = "Régression linéaire avec prévision",
       subtitle = paste("Y =", round(alpha_hat, 4), "+", round(beta_hat, 4), "X")) +
  theme_minimal()

## 10. Différents types de modèles et interprétation des coefficients
#
# Comme indiqué dans le cours, il existe différentes formes de modèles linéaires selon la transformation 
# appliquée aux variables. Nous allons explorer quatre types de modèles :
#
# 1. lin-lin: Y = α + β X + u
# 2. log-lin: log(Y) = α + β X + u
# 3. lin-log: Y = α + β log(X) + u
# 4. log-log: log(Y) = α + β log(X) + u

# Création de nouvelles variables avec des distributions positives
set.seed(0)
X_pos <- rlnorm(n, 0, 0.5)  # Distribution log-normale pour assurer X > 0
beta_log_lin <- 0.05  # Pour log-lin
beta_lin_log <- 5     # Pour lin-log
beta_log_log <- 0.8   # Pour log-log

# Génération des Y pour les différents modèles
Y_lin_lin <- 10 + 2 * X_pos + rnorm(n, 0, 1)
Y_log_lin <- exp(2 + beta_log_lin * X_pos + rnorm(n, 0, 0.1))
Y_lin_log <- 5 + beta_lin_log * log(X_pos) + rnorm(n, 0, 1)
Y_log_log <- exp(1 + beta_log_log * log(X_pos) + rnorm(n, 0, 0.1))

# Création d'un data frame pour ces données
df_models <- data.frame(
  X = X_pos,
  Y_lin_lin = Y_lin_lin,
  Y_log_lin = Y_log_lin,
  Y_lin_log = Y_lin_log,
  Y_log_log = Y_log_log,
  log_X = log(X_pos),
  log_Y_log_lin = log(Y_log_lin),
  log_Y_log_log = log(Y_log_log)
)

# Affichage des premières lignes du data frame
head(df_models)

# Estimation et interprétation des différents modèles
model_lin_lin <- lm(Y_lin_lin ~ X, data = df_models)
model_log_lin <- lm(log(Y_log_lin) ~ X, data = df_models)
model_lin_log <- lm(Y_lin_log ~ log(X), data = df_models)
model_log_log <- lm(log(Y_log_log) ~ log(X), data = df_models)

# Affichage des résultats
cat("Résumé du modèle lin-lin:\n")
print(summary(model_lin_lin))

cat("\nRésumé du modèle log-lin:\n")
print(summary(model_log_lin))

cat("\nRésumé du modèle lin-log:\n")
print(summary(model_lin_log))

cat("\nRésumé du modèle log-log:\n")
print(summary(model_log_log))

# Interprétation des coefficients
cat("\nInterprétation des coefficients:\n")
cat("1. Modèle lin-lin: Une augmentation de 1 unité de X est associée à une variation de", 
    round(coef(model_lin_lin)[2], 4), "unités de Y.\n")

cat("2. Modèle log-lin: Une augmentation de 1 unité de X est associée à une variation de", 
    round(100 * coef(model_log_lin)[2], 4), "% de Y.\n")

cat("3. Modèle lin-log: Une augmentation de 1% de X est associée à une variation de", 
    round(coef(model_lin_log)[2] / 100, 4), "unités de Y.\n")

cat("4. Modèle log-log: Une augmentation de 1% de X est associée à une variation de", 
    round(coef(model_log_log)[2], 4), "% de Y (élasticité).\n")

# Visualisation des différents modèles
par(mfrow = c(2, 2))

# 1. lin-lin
plot(df_models$X, df_models$Y_lin_lin, main = "Modèle lin-lin", 
     xlab = "X", ylab = "Y", pch = 19)
abline(model_lin_lin, col = "red", lwd = 2)

# 2. log-lin
plot(df_models$X, log(df_models$Y_log_lin), main = "Modèle log-lin", 
     xlab = "X", ylab = "log(Y)", pch = 19)
abline(model_log_lin, col = "red", lwd = 2)

# 3. lin-log
plot(log(df_models$X), df_models$Y_lin_log, main = "Modèle lin-log", 
     xlab = "log(X)", ylab = "Y", pch = 19)
abline(model_lin_log, col = "red", lwd = 2)

# 4. log-log
plot(log(df_models$X), log(df_models$Y_log_log), main = "Modèle log-log", 
     xlab = "log(X)", ylab = "log(Y)", pch = 19)
abline(model_log_log, col = "red", lwd = 2)

