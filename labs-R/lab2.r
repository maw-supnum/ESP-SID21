# TP Estimation Statistique Avancée
#
# Dans ce TP, nous allons explorer différentes méthodes d'estimation statistique :
# * Estimation ponctuelle (méthode des moments et maximum de vraisemblance)
# * Estimation par intervalle de confiance pour une moyenne (variance connue et inconnue)

# Importation des bibliothèques nécessaires
library(ggplot2)
library(stats)
library(optimx)

# Configuration pour la reproductibilité
set.seed(0)

# --------------------------------------------------------
# 1. Estimation ponctuelle
#
# L'estimation ponctuelle consiste à estimer un paramètre inconnu d'une distribution 
# par une valeur unique calculée à partir d'un échantillon. Nous allons explorer 
# deux méthodes principales :
# * Méthode des moments
# * Maximum de vraisemblance
# --------------------------------------------------------

# --------------------------------------------------------
# 1.1 Méthode des moments
#
# La méthode des moments consiste à égaler les moments théoriques de la distribution 
# aux moments empiriques calculés à partir de l'échantillon.
#
# Prenons l'exemple d'une distribution exponentielle où nous voulons estimer le paramètre λ.
# --------------------------------------------------------

# Génération d'un échantillon suivant une loi exponentielle
lambda_true <- 0.5  # Valeur réelle du paramètre λ
taille_echantillon <- 1000
echantillon_exp <- rexp(n = taille_echantillon, rate = lambda_true)

# Visualisation de l'échantillon
hist(echantillon_exp, breaks = 30, freq = FALSE, 
     main = "Distribution exponentielle", 
     xlab = "x", ylab = "Densité", 
     col = rgb(0, 0, 1, 0.7))

# Superposition de la densité théorique
curve(lambda_true * exp(-lambda_true * x), 
      from = 0, to = 15, 
      col = "red", lwd = 2, add = TRUE)

legend("topright", 
       legend = c("Échantillon", paste0("Densité théorique (λ=", lambda_true, ")")), 
       col = c(rgb(0, 0, 1, 0.7), "red"), 
       lty = c(NA, 1), lwd = c(NA, 2), pch = c(15, NA), 
       pt.cex = 2)

# Estimation de λ par la méthode des moments
# Pour une loi exponentielle, E[X] = 1/λ
# Donc λ = 1/E[X]
lambda_moments <- 1 / mean(echantillon_exp)

cat("Valeur réelle de λ :", lambda_true, "\n")
cat("Estimation de λ par la méthode des moments :", sprintf("%.6f", lambda_moments), "\n")
cat("Erreur relative :", sprintf("%.2f%%", 100 * abs(lambda_moments - lambda_true) / lambda_true), "\n")

# Illustrons maintenant la méthode des moments pour une distribution normale où nous voulons 
# estimer μ et σ².

# Génération d'un échantillon suivant une loi normale
mu_true <- 3
sigma_true <- 1.5
echantillon_norm <- rnorm(n = taille_echantillon, mean = mu_true, sd = sigma_true)

# Visualisation de l'échantillon
hist(echantillon_norm, breaks = 30, freq = FALSE, 
     main = "Distribution normale", 
     xlab = "x", ylab = "Densité", 
     col = rgb(0, 0, 1, 0.7))

# Superposition de la densité théorique
curve(dnorm(x, mean = mu_true, sd = sigma_true), 
      from = mu_true - 4*sigma_true, to = mu_true + 4*sigma_true, 
      col = "red", lwd = 2, add = TRUE)

legend("topright", 
       legend = c("Échantillon", paste0("Densité théorique (μ=", mu_true, ", σ=", sigma_true, ")")), 
       col = c(rgb(0, 0, 1, 0.7), "red"), 
       lty = c(NA, 1), lwd = c(NA, 2), pch = c(15, NA), 
       pt.cex = 2)

# Estimation par la méthode des moments
# Pour une loi normale, E[X] = μ et Var(X) = σ²
mu_moments <- mean(echantillon_norm)
sigma2_moments <- var(echantillon_norm) * (taille_echantillon - 1) / taille_echantillon  # Variance biaisée
sigma_moments <- sqrt(sigma2_moments)

cat("Valeurs réelles : μ =", mu_true, ", σ =", sigma_true, ", σ² =", sigma_true^2, "\n")
cat("Estimations par la méthode des moments :\n")
cat("  μ =", sprintf("%.6f", mu_moments), "\n")
cat("  σ =", sprintf("%.6f", sigma_moments), "\n")
cat("  σ² =", sprintf("%.6f", sigma2_moments), "\n")
cat("Erreurs relatives : μ =", sprintf("%.2f%%", 100 * abs(mu_moments - mu_true) / mu_true), 
    ", σ =", sprintf("%.2f%%", 100 * abs(sigma_moments - sigma_true) / sigma_true), "\n")

# --------------------------------------------------------
# 1.2 Maximum de vraisemblance
#
# Le principe du maximum de vraisemblance consiste à choisir comme estimateur 
# la valeur du paramètre qui maximise la probabilité d'observer l'échantillon donné.
#
# Commençons par l'exemple de la loi exponentielle.
# --------------------------------------------------------

# Fonction de log-vraisemblance négative pour une loi exponentielle
log_vraisemblance_exp <- function(lambda_val, data) {
  # log L(λ) = n log(λ) - λ Σx_i
  n <- length(data)
  # On minimise -log L
  return(-1 * (n * log(lambda_val) - lambda_val * sum(data)))
}

# Estimation de λ par maximum de vraisemblance
resultat <- optimize(f = log_vraisemblance_exp, 
                     interval = c(1e-10, 10), 
                     data = echantillon_exp)
lambda_mle <- resultat$minimum

# Comparaison avec la méthode des moments
cat("Valeur réelle de λ :", lambda_true, "\n")
cat("Estimation de λ par maximum de vraisemblance :", sprintf("%.6f", lambda_mle), "\n")
cat("Estimation de λ par la méthode des moments :", sprintf("%.6f", lambda_moments), "\n")
cat("Erreur relative (MLE) :", sprintf("%.2f%%", 100 * abs(lambda_mle - lambda_true) / lambda_true), "\n")

# Pour une loi exponentielle, l'estimateur du maximum de vraisemblance est en fait identique 
# à l'estimateur des moments : λ = 1/x̄.
#
# Vérifions :

# Calcul direct de l'estimateur du maximum de vraisemblance
lambda_mle_direct <- 1 / mean(echantillon_exp)
cat("MLE direct :", sprintf("%.6f", lambda_mle_direct), "\n")
cat("MLE par optimisation :", sprintf("%.6f", lambda_mle), "\n")
cat("Différence :", sprintf("%.10f", abs(lambda_mle_direct - lambda_mle)), "\n")

# Maintenant, utilisons le maximum de vraisemblance pour estimer les paramètres d'une distribution normale.

# Fonction de log-vraisemblance négative pour une loi normale
log_vraisemblance_norm <- function(params, data) {
  mu <- params[1]
  sigma <- params[2]
  n <- length(data)
  # log L(μ,σ) = -n/2 log(2π) - n log(σ) - 1/(2σ²) Σ(x_i - μ)²
  return(-1 * (-n/2 * log(2*pi) - n * log(sigma) - 1/(2*sigma^2) * sum((data - mu)^2)))
}

# Estimation par maximum de vraisemblance
resultat <- optim(par = c(0, 1), 
                  fn = log_vraisemblance_norm, 
                  data = echantillon_norm, 
                  method = "L-BFGS-B", 
                  lower = c(-Inf, 1e-10), 
                  upper = c(Inf, Inf))
mu_mle <- resultat$par[1]
sigma_mle <- resultat$par[2]

# Comparaison avec la méthode des moments
cat("Valeurs réelles : μ =", mu_true, ", σ =", sigma_true, "\n")
cat("Estimations par maximum de vraisemblance : μ =", sprintf("%.6f", mu_mle), 
    ", σ =", sprintf("%.6f", sigma_mle), "\n")
cat("Estimations par la méthode des moments : μ =", sprintf("%.6f", mu_moments), 
    ", σ =", sprintf("%.6f", sigma_moments), "\n")

# Pour une loi normale, les estimateurs du maximum de vraisemblance pour μ et σ sont :
# - μ̂ = 1/n Σx_i (la moyenne empirique)
# - σ̂² = 1/n Σ(x_i - μ̂)² (la variance biaisée)
#
# Vérifions :

# Calcul direct des estimateurs du maximum de vraisemblance
mu_mle_direct <- mean(echantillon_norm)
# Notez que c'est la variance biaisée (divisée par n et non par n-1)
sigma_mle_direct <- sqrt(sum((echantillon_norm - mu_mle_direct)^2) / length(echantillon_norm))

cat("MLE direct : μ =", sprintf("%.6f", mu_mle_direct), 
    ", σ =", sprintf("%.6f", sigma_mle_direct), "\n")
cat("MLE par optimisation : μ =", sprintf("%.6f", mu_mle), 
    ", σ =", sprintf("%.6f", sigma_mle), "\n")

# Visualisons la fonction de vraisemblance pour la loi normale :

# Création d'une grille de valeurs pour μ et σ
mu_range <- seq(mu_true - 1, mu_true + 1, length.out = 100)
sigma_range <- seq(sigma_true - 0.5, sigma_true + 0.5, length.out = 100)
grid_values <- expand.grid(mu = mu_range, sigma = sigma_range)

# Calcul de la log-vraisemblance pour chaque point de la grille
log_vrais <- numeric(nrow(grid_values))
for(i in 1:nrow(grid_values)) {
  log_vrais[i] <- -log_vraisemblance_norm(c(grid_values$mu[i], grid_values$sigma[i]), echantillon_norm)
}

# Visualisation en 2D (courbes de niveau)
contour_data <- data.frame(grid_values, log_vrais = log_vrais)
ggplot(contour_data, aes(x = mu, y = sigma, z = log_vrais)) +
  geom_contour_filled() +
  geom_point(aes(x = mu_mle, y = sigma_mle), color = "red", size = 3) +
  labs(title = "Courbes de niveau de la Log-Vraisemblance",
       x = "μ", y = "σ", fill = "Log-Vraisemblance") +
  theme_minimal()

# --------------------------------------------------------
# 2. Estimation par intervalle de confiance
#
# Un intervalle de confiance fournit une plage de valeurs plausibles pour un paramètre 
# inconnu, avec un certain niveau de confiance (généralement 95%).
# --------------------------------------------------------

# --------------------------------------------------------
# 2.1 Intervalle de confiance pour la moyenne d'une loi normale (variance connue)
#
# Supposons que nous connaissons la variance σ² de la population et que nous voulons 
# estimer la moyenne μ.
# --------------------------------------------------------

# Génération d'un nouvel échantillon
mu_true <- 10
sigma_true <- 2
taille_echantillon <- 50  # Échantillon plus petit pour mieux illustrer
echantillon <- rnorm(n = taille_echantillon, mean = mu_true, sd = sigma_true)

# Calcul de la moyenne échantillonnale
moyenne_echantillon <- mean(echantillon)

# Niveau de confiance
alpha <- 0.05  # 95% de confiance
z_alpha_2 <- qnorm(1 - alpha/2)  # Quantile de la loi normale centrée réduite

# Erreur standard (écart-type de la moyenne échantillonnale)
erreur_standard <- sigma_true / sqrt(taille_echantillon)

# Intervalle de confiance
ic_inf <- moyenne_echantillon - z_alpha_2 * erreur_standard
ic_sup <- moyenne_echantillon + z_alpha_2 * erreur_standard

cat("Moyenne réelle :", mu_true, "\n")
cat("Moyenne échantillonnale :", sprintf("%.4f", moyenne_echantillon), "\n")
cat("Intervalle de confiance à 95% : [", sprintf("%.4f", ic_inf), 
    ", ", sprintf("%.4f", ic_sup), "]\n", sep = "")
cat("Largeur de l'intervalle :", sprintf("%.4f", ic_sup - ic_inf), "\n")
cat("La moyenne réelle est dans l'intervalle :", ic_inf <= mu_true && mu_true <= ic_sup, "\n")

# Visualisons l'intervalle de confiance :

# Densité de la loi normale de la moyenne échantillonnale
x <- seq(moyenne_echantillon - 4*erreur_standard, 
         moyenne_echantillon + 4*erreur_standard, 
         length.out = 1000)
y <- dnorm(x, mean = moyenne_echantillon, sd = erreur_standard)

# Création du graphique
plot(x, y, type = "l", col = "blue", lwd = 2, 
     main = "Intervalle de confiance pour la moyenne (variance connue)",
     xlab = "Valeur", ylab = "Densité")

# Coloration de l'intervalle de confiance
x_ic <- seq(ic_inf, ic_sup, length.out = 100)
y_ic <- dnorm(x_ic, mean = moyenne_echantillon, sd = erreur_standard)
polygon(c(x_ic, rev(x_ic)), c(y_ic, rep(0, length(y_ic))), 
        col = rgb(0.5, 0.8, 1, 0.4), border = NA)

# Ajout des lignes verticales
abline(v = mu_true, col = "red", lty = 2)
abline(v = moyenne_echantillon, col = "green", lty = 1)
abline(v = ic_inf, col = "blue", lty = 3)
abline(v = ic_sup, col = "blue", lty = 3)

# Ajout de la légende
legend("topright", 
       legend = c("Distribution de la moyenne échantillonnale", 
                  "Intervalle de confiance à 95%", 
                  "Moyenne réelle", 
                  "Moyenne échantillonnale", 
                  "Bornes de l'intervalle"),
       col = c("blue", rgb(0.5, 0.8, 1, 0.4), "red", "green", "blue"),
       lty = c(1, NA, 2, 1, 3), 
       lwd = c(2, NA, 1, 1, 1),
       pch = c(NA, 15, NA, NA, NA),
       pt.cex = 2)

# Illustrons l'effet de la taille de l'échantillon sur la largeur de l'intervalle de confiance :

# Différentes tailles d'échantillon
tailles <- c(10, 30, 50, 100, 200, 500, 1000)
largeurs <- numeric(length(tailles))

for(i in 1:length(tailles)) {
  n <- tailles[i]
  # Erreur standard
  erreur_std <- sigma_true / sqrt(n)
  # Largeur de l'intervalle de confiance
  largeur <- 2 * z_alpha_2 * erreur_std
  largeurs[i] <- largeur
}

# Création du graphique
plot(tailles, largeurs, type = "b", col = "blue", pch = 16, lwd = 2, 
     log = "x", # échelle logarithmique pour l'axe x
     main = "Largeur de l'intervalle de confiance en fonction de la taille de l'échantillon",
     xlab = "Taille de l'échantillon", ylab = "Largeur de l'intervalle")
grid()

# --------------------------------------------------------
# 2.2 Intervalle de confiance pour la moyenne d'une loi normale (variance inconnue)
#
# Lorsque la variance σ² est inconnue, nous l'estimons à partir de l'échantillon 
# et utilisons la loi de Student (t) au lieu de la loi normale.
# --------------------------------------------------------

# Utilisation du même échantillon que précédemment, mais nous "oublions" la variance réelle
# et l'estimons à partir de l'échantillon

# Estimation de la variance (non biaisée)
variance_estimee <- var(echantillon)
ecart_type_estime <- sqrt(variance_estimee)

# Erreur standard estimée
erreur_standard_estimee <- ecart_type_estime / sqrt(taille_echantillon)

# Quantile de la loi de Student
t_alpha_2 <- qt(1 - alpha/2, df = taille_echantillon - 1)  # df = degrés de liberté = n-1

# Intervalle de confiance selon la loi de Student
ic_inf_t <- moyenne_echantillon - t_alpha_2 * erreur_standard_estimee
ic_sup_t <- moyenne_echantillon + t_alpha_2 * erreur_standard_estimee

cat("Moyenne réelle :", mu_true, "\n")
cat("Moyenne échantillonnale :", sprintf("%.4f", moyenne_echantillon), "\n")
cat("Écart-type réel :", sigma_true, "\n")
cat("Écart-type estimé :", sprintf("%.4f", ecart_type_estime), "\n")
cat("Intervalle de confiance à 95% (loi normale) : [", sprintf("%.4f", ic_inf), 
    ", ", sprintf("%.4f", ic_sup), "]\n", sep = "")
cat("Intervalle de confiance à 95% (loi de Student) : [", sprintf("%.4f", ic_inf_t), 
    ", ", sprintf("%.4f", ic_sup_t), "]\n", sep = "")

# Comparons les deux approches (variance connue vs. inconnue) :

# Densité de la loi normale
x <- seq(moyenne_echantillon - 4*erreur_standard, 
         moyenne_echantillon + 4*erreur_standard, 
         length.out = 1000)
y_norm <- dnorm(x, mean = moyenne_echantillon, sd = erreur_standard)

# Densité de la loi de Student (ajustée pour la comparaison)
y_t <- dt((x - moyenne_echantillon) / erreur_standard_estimee, 
          df = taille_echantillon - 1) / erreur_standard_estimee

# Création du graphique
plot(x, y_norm, type = "l", col = "blue", lwd = 2, 
     main = "Comparaison des intervalles de confiance (variance connue vs. inconnue)",
     xlab = "Valeur", ylab = "Densité", 
     ylim = c(0, max(c(y_norm, y_t))))
lines(x, y_t, col = "red", lwd = 2)

# Coloration des intervalles de confiance
rect(ic_inf, 0, ic_sup, max(y_norm), 
     col = rgb(0, 0, 1, 0.2), border = NA)
rect(ic_inf_t, 0, ic_sup_t, max(y_t), 
     col = rgb(1, 0, 0, 0.2), border = NA)

# Ajout des lignes verticales
abline(v = mu_true, col = "black", lty = 2)
abline(v = moyenne_echantillon, col = "green", lty = 1)

# Ajout de la légende
legend("topright", 
       legend = c("Loi normale (variance connue)", 
                  "Loi de Student (variance inconnue)",
                  "IC Loi normale",
                  "IC Loi de Student",
                  "Moyenne réelle",
                  "Moyenne échantillonnale"),
       col = c("blue", "red", rgb(0, 0, 1, 0.2), rgb(1, 0, 0, 0.2), "black", "green"),
       lty = c(1, 1, NA, NA, 2, 1), 
       lwd = c(2, 2, NA, NA, 1, 1),
       pch = c(NA, NA, 15, 15, NA, NA),
       pt.cex = 2)

# Voyons l'impact de la taille de l'échantillon sur la différence entre les deux approches :

# Différentes tailles d'échantillon
tailles <- c(5, 10, 20, 30, 50, 100, 500)
largeurs_norm <- numeric(length(tailles))
largeurs_t <- numeric(length(tailles))

for(i in 1:length(tailles)) {
  n <- tailles[i]
  # Erreur standard
  erreur_std <- sigma_true / sqrt(n)
  # Largeur de l'intervalle de confiance (variance connue)
  largeur_norm <- 2 * z_alpha_2 * erreur_std
  largeurs_norm[i] <- largeur_norm
  
  # Quantile de Student avec n-1 degrés de liberté
  t_val <- qt(1 - alpha/2, df = n - 1)
  # Erreur standard estimée
  erreur_std_estimee <- ecart_type_estime / sqrt(n)
  # Largeur de l'intervalle de confiance (variance inconnue)
  largeur_t <- 2 * t_val * erreur_std_estimee
  largeurs_t[i] <- largeur_t
}

# Création du graphique
plot(tailles, largeurs_norm, type = "b", col = "blue", pch = 16, lwd = 2, 
     log = "x", # échelle logarithmique pour l'axe x
     ylim = range(c(largeurs_norm, largeurs_t)),
     main = "Largeur de l'intervalle de confiance en fonction de la taille de l'échantillon",
     xlab = "Taille de l'échantillon", ylab = "Largeur de l'intervalle")
lines(tailles, largeurs_t, type = "b", col = "red", pch = 16, lwd = 2)
legend("topright", 
       legend = c("Variance connue", "Variance inconnue"),
       col = c("blue", "red"),
       lty = 1, pch = 16, lwd = 2)
grid()

