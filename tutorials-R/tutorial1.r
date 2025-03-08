# TP Estimation Statistique Avancée (Suite)
#
# Dans cette partie du TP, nous allons explorer deux méthodes d'estimation par intervalle de confiance supplémentaires :
# * Estimation par intervalle de confiance pour une proportion
# * Estimation par intervalle de confiance pour une variance (moyenne connue et inconnue)

# Importation des bibliothèques nécessaires et configuration
set.seed(0) # Pour la reproductibilité

# -----------------------------------------------------------------------------
# 1. Estimation par intervalle de confiance pour une proportion
# -----------------------------------------------------------------------------
# Lorsque l'on veut estimer une proportion p dans une population, on utilise généralement la statistique 
# p_hat = X/n où X est le nombre de succès dans un échantillon de taille n. 
# Sous certaines conditions (n assez grand), cette statistique suit approximativement une loi normale.

# Génération d'un échantillon binomial (succès/échec)
p_true <- 0.3  # Proportion réelle
taille_echantillon <- 100
echantillon_bin <- rbinom(taille_echantillon, 1, p_true)

# Estimation de la proportion
p_hat <- mean(echantillon_bin)

cat("Proportion réelle :", p_true, "\n")
cat("Proportion estimée :", p_hat, "\n")
cat("Nombre de succès :", p_hat * taille_echantillon, "\n")

# -----------------------------------------------------------------------------
# 1.1 Intervalle de confiance avec l'approximation normale
# -----------------------------------------------------------------------------
# Pour n assez grand (généralement n*p ≥ 5 et n*(1-p) ≥ 5), on peut utiliser l'approximation normale.
# L'erreur standard de p_hat est donnée par sqrt(p(1-p)/n).

# Niveau de confiance
alpha <- 0.05  # 95% de confiance
z_alpha_2 <- qnorm(1 - alpha/2)  # Quantile de la loi normale centrée réduite

# Erreur standard estimée
erreur_standard <- sqrt(p_hat * (1 - p_hat) / taille_echantillon)

# Intervalle de confiance
ic_inf <- p_hat - z_alpha_2 * erreur_standard
ic_sup <- p_hat + z_alpha_2 * erreur_standard

cat("Proportion réelle :", p_true, "\n")
cat("Proportion estimée :", round(p_hat, 4), "\n")
cat("Intervalle de confiance à 95% : [", max(0, round(ic_inf, 4)), ", ", min(1, round(ic_sup, 4)), "]\n")
cat("Largeur de l'intervalle :", round(ic_sup - ic_inf, 4), "\n")
cat("La proportion réelle est dans l'intervalle :", ic_inf <= p_true && p_true <= ic_sup, "\n")

# Visualisons l'intervalle de confiance
visualiser_IC_proportion <- function() {
  x <- seq(max(0, p_hat - 4*erreur_standard), min(1, p_hat + 4*erreur_standard), length.out = 1000)
  y <- dnorm(x, p_hat, erreur_standard)
  
  plot(x, y, type = "l", col = "blue", lwd = 2, 
       main = "Intervalle de confiance pour une proportion",
       xlab = "Valeur", ylab = "Densité",
       xlim = c(max(0, p_hat - 4*erreur_standard), min(1, p_hat + 4*erreur_standard)))
  
  # Coloration de l'intervalle de confiance
  x_ic <- seq(max(0, ic_inf), min(1, ic_sup), length.out = 100)
  y_ic <- dnorm(x_ic, p_hat, erreur_standard)
  polygon(c(x_ic, rev(x_ic)), c(rep(0, length(x_ic)), rev(y_ic)), 
          col = rgb(0.53, 0.81, 0.92, 0.4), border = NA)
  
  # Ajout des lignes verticales
  abline(v = p_true, col = "red", lty = 2)
  abline(v = p_hat, col = "green", lty = 1)
  abline(v = max(0, ic_inf), col = "blue", lty = 3)
  abline(v = min(1, ic_sup), col = "blue", lty = 3)
  
  # Légende
  legend("topright", 
         legend = c("Distribution de la proportion échantillonnale", 
                    "Intervalle de confiance à 95%", 
                    "Proportion réelle",
                    "Proportion échantillonnale",
                    "Bornes de l'intervalle"),
         col = c("blue", rgb(0.53, 0.81, 0.92), "red", "green", "blue"),
         lty = c(1, 1, 2, 1, 3), lwd = c(2, 10, 1, 1, 1),
         cex = 0.8)
}

# Pour exécuter et visualiser le graphique, décommentez la ligne suivante:
visualiser_IC_proportion()

# -----------------------------------------------------------------------------
# 1.2 Correction de continuité (Intervalle de Wilson)
# -----------------------------------------------------------------------------
# L'intervalle de Wilson est une méthode plus précise pour estimer l'intervalle de confiance d'une proportion,
# notamment pour les petits échantillons ou les proportions proches de 0 ou 1.

# Calcul de l'intervalle de Wilson
# Nombre de succès
X <- p_hat * taille_echantillon

# Terme central
centre <- (X + z_alpha_2^2/2) / (taille_echantillon + z_alpha_2^2)

# Terme de marge
marge <- z_alpha_2 * sqrt((X * (taille_echantillon - X) / taille_echantillon + z_alpha_2^2/4) / 
                         (taille_echantillon + z_alpha_2^2))

# Limites de l'intervalle
ic_inf_wilson <- max(0, centre - marge)
ic_sup_wilson <- min(1, centre + marge)

cat("Proportion réelle :", p_true, "\n")
cat("Proportion estimée :", round(p_hat, 4), "\n")
cat("Intervalle de confiance normal à 95% : [", max(0, round(ic_inf, 4)), ", ", min(1, round(ic_sup, 4)), "]\n")
cat("Intervalle de Wilson à 95% : [", round(ic_inf_wilson, 4), ", ", round(ic_sup_wilson, 4), "]\n")
cat("La proportion réelle est dans l'intervalle normal :", ic_inf <= p_true && p_true <= ic_sup, "\n")
cat("La proportion réelle est dans l'intervalle de Wilson :", ic_inf_wilson <= p_true && p_true <= ic_sup_wilson, "\n")

# Comparons les deux méthodes pour différentes tailles d'échantillon
comparer_methodes <- function() {
  # Différentes tailles d'échantillon
  tailles <- c(10, 20, 30, 50, 100, 200, 500, 1000)
  largeurs_norm <- numeric(length(tailles))
  largeurs_wilson <- numeric(length(tailles))
  couverture_norm <- numeric(length(tailles))
  couverture_wilson <- numeric(length(tailles))
  
  n_simulations <- 1000  # Nombre de simulations pour chaque taille
  
  for (i in 1:length(tailles)) {
    n <- tailles[i]
    compteur_norm <- 0
    compteur_wilson <- 0
    largeur_norm_total <- 0
    largeur_wilson_total <- 0
    
    for (j in 1:n_simulations) {
      # Génération d'un échantillon
      echantillon <- rbinom(n, 1, p_true)
      p_est <- mean(echantillon)
      X <- p_est * n
      
      # Méthode normale
      err_std <- sqrt(p_est * (1 - p_est) / n)
      ic_inf_n <- max(0, p_est - z_alpha_2 * err_std)
      ic_sup_n <- min(1, p_est + z_alpha_2 * err_std)
      
      # Méthode de Wilson
      centre <- (X + z_alpha_2^2/2) / (n + z_alpha_2^2)
      marge <- z_alpha_2 * sqrt((X * (n - X) / n + z_alpha_2^2/4) / (n + z_alpha_2^2))
      ic_inf_w <- max(0, centre - marge)
      ic_sup_w <- min(1, centre + marge)
      
      # Vérification de la couverture
      if (ic_inf_n <= p_true && p_true <= ic_sup_n) {
        compteur_norm <- compteur_norm + 1
      }
      if (ic_inf_w <= p_true && p_true <= ic_sup_w) {
        compteur_wilson <- compteur_wilson + 1
      }
      
      # Calcul des largeurs
      largeur_norm_total <- largeur_norm_total + (ic_sup_n - ic_inf_n)
      largeur_wilson_total <- largeur_wilson_total + (ic_sup_w - ic_inf_w)
    }
    
    # Moyenne des résultats
    couverture_norm[i] <- compteur_norm / n_simulations
    couverture_wilson[i] <- compteur_wilson / n_simulations
    largeurs_norm[i] <- largeur_norm_total / n_simulations
    largeurs_wilson[i] <- largeur_wilson_total / n_simulations
  }
  
  # Configuration pour l'affichage en 1x2
  par(mfrow = c(1, 2))
  
  # Taux de couverture
  plot(tailles, couverture_norm, type = "b", col = "blue", pch = 16, lwd = 2, 
       xlab = "Taille de l'échantillon", ylab = "Taux de couverture",
       main = "Taux de couverture des intervalles de confiance",
       log = "x", ylim = c(0.7, 1))
  points(tailles, couverture_wilson, type = "b", col = "red", pch = 16, lwd = 2)
  abline(h = 1-alpha, col = "black", lty = 2)
  legend("bottomright", 
         legend = c("Méthode normale", "Méthode de Wilson", "Niveau nominal (95%)"),
         col = c("blue", "red", "black"),
         lty = c(1, 1, 2), pch = c(16, 16, NA), lwd = c(2, 2, 1),
         cex = 0.8)
  
  # Largeur des intervalles
  plot(tailles, largeurs_norm, type = "b", col = "blue", pch = 16, lwd = 2, 
       xlab = "Taille de l'échantillon", ylab = "Largeur moyenne",
       main = "Largeur moyenne des intervalles de confiance",
       log = "x")
  points(tailles, largeurs_wilson, type = "b", col = "red", pch = 16, lwd = 2)
  legend("topright", 
         legend = c("Méthode normale", "Méthode de Wilson"),
         col = c("blue", "red"),
         lty = c(1, 1), pch = c(16, 16), lwd = c(2, 2),
         cex = 0.8)
  
  # Réinitialiser la configuration
  par(mfrow = c(1, 1))
}

# Pour exécuter et visualiser les graphiques comparatifs, décommentez la ligne suivante:
comparer_methodes()

# -----------------------------------------------------------------------------
# 2. Estimation par intervalle de confiance pour la variance d'une loi normale
# -----------------------------------------------------------------------------
# L'estimation de la variance sigma² d'une loi normale repose sur la propriété suivante:
# si X1, X2, ..., Xn est un échantillon aléatoire d'une distribution normale de moyenne µ et de variance sigma²,
# alors la statistique (n-1)S²/sigma² suit une loi du χ² à (n-1) degrés de liberté,
# où S² est la variance échantillonnale non biaisée.

# -----------------------------------------------------------------------------
# 2.1 Intervalle de confiance pour la variance (moyenne connue)
# -----------------------------------------------------------------------------
# Lorsque la moyenne µ est connue, la statistique sum((Xi - µ)²)/sigma² suit une loi du χ² à n degrés de liberté.

# Génération d'un échantillon suivant une loi normale
mu_true <- 10
sigma_true <- 2
taille_echantillon <- 50
echantillon <- rnorm(taille_echantillon, mean = mu_true, sd = sigma_true)

# Calcul de la variance biaisée (avec mu connu)
variance_biaisee_mu_connu <- sum((echantillon - mu_true)^2) / taille_echantillon

# Niveau de confiance
alpha <- 0.05  # 95% de confiance

# Quantiles de la loi du chi-2 à n degrés de liberté
chi2_inf <- qchisq(alpha/2, df = taille_echantillon)
chi2_sup <- qchisq(1 - alpha/2, df = taille_echantillon)

# Intervalle de confiance pour sigma²
ic_inf_var <- taille_echantillon * variance_biaisee_mu_connu / chi2_sup
ic_sup_var <- taille_echantillon * variance_biaisee_mu_connu / chi2_inf

# Intervalle de confiance pour sigma
ic_inf_sigma <- sqrt(ic_inf_var)
ic_sup_sigma <- sqrt(ic_sup_var)

cat("Variance réelle :", sigma_true^2, "\n")
cat("Variance estimée (moyenne connue) :", round(variance_biaisee_mu_connu, 4), "\n")
cat("Intervalle de confiance à 95% pour la variance : [", round(ic_inf_var, 4), ", ", round(ic_sup_var, 4), "]\n")
cat("Écart-type réel :", sigma_true, "\n")
cat("Intervalle de confiance à 95% pour l'écart-type : [", round(ic_inf_sigma, 4), ", ", round(ic_sup_sigma, 4), "]\n")
cat("La variance réelle est dans l'intervalle :", ic_inf_var <= sigma_true^2 && sigma_true^2 <= ic_sup_var, "\n")

# Visualisons l'intervalle de confiance pour la variance
visualiser_IC_variance <- function() {
  x <- seq(0.1, chi2_sup + 10, length.out = 1000)
  stat <- taille_echantillon * variance_biaisee_mu_connu / sigma_true^2
  y <- dchisq(x, df = taille_echantillon)
  
  plot(x, y, type = "l", col = "blue", lwd = 2, 
       main = "Distribution du Chi² pour l'intervalle de confiance de la variance (moyenne connue)",
       xlab = "Valeur de la statistique (n*sigma²_hat)/(sigma²)",
       ylab = "Densité",
       xlim = c(0, chi2_sup + 10))
  
  # Coloration des zones de rejet
  x_reject_left <- seq(0.1, chi2_inf, length.out = 100)
  y_reject_left <- dchisq(x_reject_left, df = taille_echantillon)
  polygon(c(x_reject_left, rev(x_reject_left)), c(rep(0, length(x_reject_left)), rev(y_reject_left)), 
          col = rgb(1, 0, 0, 0.3), border = NA)
  
  x_reject_right <- seq(chi2_sup, chi2_sup + 10, length.out = 100)
  y_reject_right <- dchisq(x_reject_right, df = taille_echantillon)
  polygon(c(x_reject_right, rev(x_reject_right)), c(rep(0, length(x_reject_right)), rev(y_reject_right)), 
          col = rgb(1, 0, 0, 0.3), border = NA)
  
  # Zone d'acceptation
  x_accept <- seq(chi2_inf, chi2_sup, length.out = 100)
  y_accept <- dchisq(x_accept, df = taille_echantillon)
  polygon(c(x_accept, rev(x_accept)), c(rep(0, length(x_accept)), rev(y_accept)), 
          col = rgb(0, 1, 0, 0.3), border = NA)
  
  # Statistique observée
  abline(v = stat, col = "black", lty = 2)
  abline(v = chi2_inf, col = "red", lty = 3)
  abline(v = chi2_sup, col = "red", lty = 3)
  
  # Légende
  legend("topright", 
         legend = c("Distribution du Chi²", 
                    "Zones de rejet (2.5%)", 
                    "Zone d'acceptation (95%)",
                    "Statistique observée"),
         col = c("blue", rgb(1, 0, 0, 0.7), rgb(0, 1, 0, 0.7), "black"),
         lty = c(1, 1, 1, 2), lwd = c(2, 10, 10, 1),
         cex = 0.8)
}

# Pour exécuter et visualiser le graphique, décommentez la ligne suivante:
visualiser_IC_variance()

# -----------------------------------------------------------------------------
# 2.2 Intervalle de confiance pour la variance (moyenne inconnue)
# -----------------------------------------------------------------------------
# Lorsque la moyenne µ est inconnue et doit être estimée, on utilise la variance échantillonnale
# S² = (1/(n-1))*sum((Xi - X̄)²) et la statistique (n-1)S²/sigma² qui suit une loi du χ² à (n-1) degrés de liberté.

# Calcul de la variance non biaisée (avec mu inconnu)
variance_non_biaisee <- var(echantillon) # var utilise n-1 au dénominateur par défaut en R

# Quantiles de la loi du chi-2 à (n-1) degrés de liberté
chi2_inf_n1 <- qchisq(alpha/2, df = taille_echantillon-1)
chi2_sup_n1 <- qchisq(1 - alpha/2, df = taille_echantillon-1)

# Intervalle de confiance pour sigma²
ic_inf_var_n1 <- (taille_echantillon - 1) * variance_non_biaisee / chi2_sup_n1
ic_sup_var_n1 <- (taille_echantillon - 1) * variance_non_biaisee / chi2_inf_n1

# Intervalle de confiance pour sigma
ic_inf_sigma_n1 <- sqrt(ic_inf_var_n1)
ic_sup_sigma_n1 <- sqrt(ic_sup_var_n1)

cat("Variance réelle :", sigma_true^2, "\n")
cat("Variance estimée (moyenne connue) :", round(variance_biaisee_mu_connu, 4), "\n")
cat("Variance estimée (moyenne inconnue) :", round(variance_non_biaisee, 4), "\n")
cat("Intervalle de confiance à 95% pour la variance (moyenne connue) : [", round(ic_inf_var, 4), ", ", round(ic_sup_var, 4), "]\n")
cat("Intervalle de confiance à 95% pour la variance (moyenne inconnue) : [", round(ic_inf_var_n1, 4), ", ", round(ic_sup_var_n1, 4), "]\n")
cat("Écart-type réel :", sigma_true, "\n")
cat("Intervalle de confiance à 95% pour l'écart-type (moyenne inconnue) : [", round(ic_inf_sigma_n1, 4), ", ", round(ic_sup_sigma_n1, 4), "]\n")
cat("La variance réelle est dans l'intervalle (moyenne connue) :", ic_inf_var <= sigma_true^2 && sigma_true^2 <= ic_sup_var, "\n")
cat("La variance réelle est dans l'intervalle (moyenne inconnue) :", ic_inf_var_n1 <= sigma_true^2 && sigma_true^2 <= ic_sup_var_n1, "\n")

# Comparons les intervalles de confiance pour la variance avec moyenne connue et inconnue
comparer_IC_variances <- function() {
  # Création d'une grille de valeurs pour sigma²
  sigma2_range <- seq(0.5, 8, length.out = 1000)
  
  # Fonction de densité de probabilité pour S² (moyenne connue)
  pdf_var_mu_connu <- function(sigma2) {
    # Transformation de variable
    # Si X ~ chi2(n), alors Y = aX suit une loi gamma
    dgamma(variance_biaisee_mu_connu * taille_echantillon / sigma2, 
           shape = taille_echantillon/2, scale = 2/taille_echantillon) * 
      (variance_biaisee_mu_connu * taille_echantillon / sigma2^2)
  }
  
  # Fonction de densité de probabilité pour S² (moyenne inconnue)
  pdf_var_mu_inconnu <- function(sigma2) {
    dgamma(variance_non_biaisee * (taille_echantillon-1) / sigma2, 
           shape = (taille_echantillon-1)/2, scale = 2/(taille_echantillon-1)) * 
      (variance_non_biaisee * (taille_echantillon-1) / sigma2^2)
  }
  
  # Calcul des densités
  y_connu <- sapply(sigma2_range, pdf_var_mu_connu)
  y_inconnu <- sapply(sigma2_range, pdf_var_mu_inconnu)
  
  # Tracé des densités
  plot(sigma2_range, y_connu, type = "l", col = "blue", lwd = 2,
       main = "Comparaison des intervalles de confiance pour la variance",
       xlab = expression(sigma^2),
       ylab = "Densité de probabilité",
       xlim = range(sigma2_range),
       ylim = c(0, max(c(y_connu, y_inconnu))))
  lines(sigma2_range, y_inconnu, col = "red", lwd = 2)
  
  # Coloration des intervalles de confiance
  rect(ic_inf_var, 0, ic_sup_var, max(c(y_connu, y_inconnu)) * 0.8, 
       col = rgb(0, 0, 1, 0.2), border = NA)
  rect(ic_inf_var_n1, 0, ic_sup_var_n1, max(c(y_connu, y_inconnu)) * 0.8, 
       col = rgb(1, 0, 0, 0.2), border = NA)
  
  # Ligne verticale pour la vraie valeur
  abline(v = sigma_true^2, col = "black", lty = 2)
  
  # Légende
  legend("topright", 
         legend = c("Densité (moyenne connue)", 
                    "Densité (moyenne inconnue)",
                    "IC (moyenne connue)",
                    "IC (moyenne inconnue)",
                    "Variance réelle"),
         col = c("blue", "red", rgb(0, 0, 1, 0.7), rgb(1, 0, 0, 0.7), "black"),
         lty = c(1, 1, 1, 1, 2), lwd = c(2, 2, 10, 10, 1),
         cex = 0.8)
}

# Pour exécuter et visualiser le graphique, décommentez la ligne suivante:
comparer_IC_variances()

# Étudions l'impact de la taille de l'échantillon sur la largeur des intervalles de confiance
etudier_impact_taille <- function() {
  # Différentes tailles d'échantillon
  tailles <- c(10, 20, 30, 50, 100, 200, 500, 1000)
  largeurs_mu_connu <- numeric(length(tailles))
  largeurs_mu_inconnu <- numeric(length(tailles))
  couverture_mu_connu <- numeric(length(tailles))
  couverture_mu_inconnu <- numeric(length(tailles))
  
  n_simulations <- 1000  # Nombre de simulations pour chaque taille
  
  for (i in 1:length(tailles)) {
    n <- tailles[i]
    compteur_mu_connu <- 0
    compteur_mu_inconnu <- 0
    largeur_mu_connu_total <- 0
    largeur_mu_inconnu_total <- 0
    
    for (j in 1:n_simulations) {
      # Génération d'un échantillon
      echantillon <- rnorm(n, mean = mu_true, sd = sigma_true)
      
      # Variance biaisée (mu connu)
      variance_biaisee <- sum((echantillon - mu_true)^2) / n
      
      # Variance non biaisée (mu inconnu)
      variance_non_biaisee <- var(echantillon)
      
      # Calcul des quantiles
      chi2_inf <- qchisq(alpha/2, df = n)
      chi2_sup <- qchisq(1 - alpha/2, df = n)
      chi2_inf_n1 <- qchisq(alpha/2, df = n-1)
      chi2_sup_n1 <- qchisq(1 - alpha/2, df = n-1)
      
      # Intervalle de confiance (moyenne connue)
      ic_inf_var <- n * variance_biaisee / chi2_sup
      ic_sup_var <- n * variance_biaisee / chi2_inf
      
      # Intervalle de confiance (moyenne inconnue)
      ic_inf_var_n1 <- (n - 1) * variance_non_biaisee / chi2_sup_n1
      ic_sup_var_n1 <- (n - 1) * variance_non_biaisee / chi2_inf_n1
      
      # Vérification de la couverture
      if (ic_inf_var <= sigma_true^2 && sigma_true^2 <= ic_sup_var) {
        compteur_mu_connu <- compteur_mu_connu + 1
      }
      if (ic_inf_var_n1 <= sigma_true^2 && sigma_true^2 <= ic_sup_var_n1) {
        compteur_mu_inconnu <- compteur_mu_inconnu + 1
      }
      
      # Calcul des largeurs
      largeur_mu_connu_total <- largeur_mu_connu_total + (ic_sup_var - ic_inf_var)
      largeur_mu_inconnu_total <- largeur_mu_inconnu_total + (ic_sup_var_n1 - ic_inf_var_n1)
    }
    
    # Moyenne des résultats
    couverture_mu_connu[i] <- compteur_mu_connu / n_simulations
    couverture_mu_inconnu[i] <- compteur_mu_inconnu / n_simulations
    largeurs_mu_connu[i] <- largeur_mu_connu_total / n_simulations
    largeurs_mu_inconnu[i] <- largeur_mu_inconnu_total / n_simulations
  }
  
  # Configuration pour l'affichage en 1x2
  par(mfrow = c(1, 2))
  
  # Taux de couverture
  plot(tailles, couverture_mu_connu, type = "b", col = "blue", pch = 16, lwd = 2, 
       xlab = "Taille de l'échantillon", ylab = "Taux de couverture",
       main = "Taux de couverture des IC pour la variance",
       log = "x", ylim = c(0.7, 1))
  points(tailles, couverture_mu_inconnu, type = "b", col = "red", pch = 16, lwd = 2)
  abline(h = 1-alpha, col = "black", lty = 2)
  legend("bottomright", 
         legend = c("Moyenne connue", "Moyenne inconnue", "Niveau nominal (95%)"),
         col = c("blue", "red", "black"),
         lty = c(1, 1, 2), pch = c(16, 16, NA), lwd = c(2, 2, 1),
         cex = 0.8)
  
  # Largeur des intervalles
  plot(tailles, largeurs_mu_connu, type = "b", col = "blue", pch = 16, lwd = 2, 
       xlab = "Taille de l'échantillon", ylab = "Largeur moyenne",
       main = "Largeur moyenne des IC pour la variance",
       log = "x")
  points(tailles, largeurs_mu_inconnu, type = "b", col = "red", pch = 16, lwd = 2)
  legend("topright", 
         legend = c("Moyenne connue", "Moyenne inconnue"),
         col = c("blue", "red"),
         lty = c(1, 1), pch = c(16, 16), lwd = c(2, 2),
         cex = 0.8)
}
etudier_impact_taille()
