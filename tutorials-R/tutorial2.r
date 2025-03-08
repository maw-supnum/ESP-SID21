# TD - Tutorial - Théorie des Tests Statistiques
#
# Nous allons explorer :
# 1. Tests de comparaison de moyennes

# Importation des bibliothèques nécessaires
library(ggplot2)

# Configuration pour les graphiques
theme_set(theme_bw())

# -------------------------------------------------------------------------
# 1. Tests de comparaison de moyennes
# -------------------------------------------------------------------------
#
# Les tests de comparaison de moyennes permettent de déterminer s'il existe 
# une différence significative entre les moyennes de deux populations.
#
# 1.1 Théorie
#
# Nous avons deux approches principales pour comparer les moyennes de deux populations :
#
# Cas 1 : Lorsque les variances σ²ₓ et σ²ᵧ sont connues
#
# La fonction pivotale (statistique de test) est :
#
# U_0 = (X̄_n_X - Ȳ_n_Y) / sqrt((σ²_X/n_X) + (σ²_Y/n_Y)) ~ N(0, 1)
#
# Les différentes hypothèses et régions critiques sont :
#
# | (H₀)                | (H₁)      | Région Critique         |
# |---------------------|-----------|-------------------------|
# | μₓ = μᵧ              | μₓ ≠ μᵧ    | RCα = {|U₀| > z_{α/2}}  |
# | μₓ = μᵧ ou μₓ ≥ μᵧ    | μₓ < μᵧ    | RCα = {U₀ < -z_{α}}     |
# | μₓ = μᵧ ou μₓ ≤ μᵧ    | μₓ > μᵧ    | RCα = {U₀ > z_{α}}      |
#
# où z_β est le quantile d'ordre 1-β de la loi normale centrée réduite.
#
# Cas 2 : Lorsque les variances σ²ₓ et σ²ᵧ sont inconnues
#
# La fonction pivotale (statistique de test) devient :
#
# T_0 = (X̄_n_X - Ȳ_n_Y) / (ν * sqrt(1/n_X + 1/n_Y)) ~ T_m
#
# Avec ν² = ((n_X-1)S²_X + (n_Y-1)S²_Y) / (n_X+n_Y-2) et m = nₓ + nᵧ - 2 degrés de liberté.
#
# Les différentes hypothèses et régions critiques sont :
#
# | (H₀)                | (H₁)      | Région Critique            |
# |---------------------|-----------|----------------------------|
# | μₓ = μᵧ              | μₓ ≠ μᵧ    | RCα = {|T₀| > t_{m,α/2}}   |
# | μₓ = μᵧ ou μₓ ≥ μᵧ    | μₓ < μᵧ    | RCα = {T₀ < -t_{m,α}}      |
# | μₓ = μᵧ ou μₓ ≤ μᵧ    | μₓ > μᵧ    | RCα = {T₀ > t_{m,α}}       |
#
# où t_{m,β} est le quantile d'ordre 1-β de la loi de Student à m degrés de liberté.

# 1.2 Implémentation des tests de comparaison de moyennes
#
# Nous allons maintenant implémenter ces tests en R et les appliquer sur des exemples.

test_moyennes_variances_connues <- function(X, Y, sigma2_X, sigma2_Y, alpha=0.05, alternative='two-sided') {
  # Test de comparaison de moyennes lorsque les variances sont connues
  #
  # Paramètres :
  # X : vecteur numérique, échantillon de la population X
  # Y : vecteur numérique, échantillon de la population Y
  # sigma2_X : numeric, variance connue de la population X
  # sigma2_Y : numeric, variance connue de la population Y
  # alpha : numeric, niveau de signification (défaut: 0.05)
  # alternative : character, type d'hypothèse alternative ('two-sided', 'less', 'greater')
  #
  # Retourne :
  # Une liste contenant :
  # U0 : numeric, valeur de la statistique de test
  # p_value : numeric, p-valeur du test
  # critical_value : numeric, valeur critique
  # reject_H0 : logical, décision de rejeter H0 ou non
  
  # Calcul des tailles d'échantillon
  n_X <- length(X)
  n_Y <- length(Y)
  
  # Calcul des moyennes d'échantillon
  X_bar <- mean(X)
  Y_bar <- mean(Y)
  
  # Calcul de la statistique de test U0
  denominateur <- sqrt(sigma2_X/n_X + sigma2_Y/n_Y)
  U0 <- (X_bar - Y_bar) / denominateur
  
  # Calcul de la valeur critique et de la p-valeur en fonction du type d'alternative
  if (alternative == 'two-sided') {
    critical_value <- qnorm(1 - alpha/2)
    p_value <- 2 * (1 - pnorm(abs(U0)))
    reject_H0 <- abs(U0) > critical_value
  } else if (alternative == 'less') {
    critical_value <- qnorm(alpha)
    p_value <- pnorm(U0)
    reject_H0 <- U0 < critical_value
  } else if (alternative == 'greater') {
    critical_value <- qnorm(1 - alpha)
    p_value <- 1 - pnorm(U0)
    reject_H0 <- U0 > critical_value
  } else {
    stop("'alternative' doit être 'two-sided', 'less' ou 'greater'")
  }
  
  return(list(U0 = U0, p_value = p_value, critical_value = critical_value, reject_H0 = reject_H0))
}

test_moyennes_variances_inconnues <- function(X, Y, alpha=0.05, alternative='two-sided', equal_var=TRUE) {
  # Test de comparaison de moyennes lorsque les variances sont inconnues
  #
  # Paramètres :
  # X : vecteur numérique, échantillon de la population X
  # Y : vecteur numérique, échantillon de la population Y
  # alpha : numeric, niveau de signification (défaut: 0.05)
  # alternative : character, type d'hypothèse alternative ('two-sided', 'less', 'greater')
  # equal_var : logical, si TRUE, suppose que les variances sont égales (défaut: TRUE)
  #
  # Retourne :
  # Une liste contenant :
  # T0 : numeric, valeur de la statistique de test
  # p_value : numeric, p-valeur du test
  # critical_value : numeric, valeur critique
  # reject_H0 : logical, décision de rejeter H0 ou non
  # nu : numeric, estimation de l'écart-type poolé (ou NULL si equal_var=FALSE)
  
  # Calcul des tailles d'échantillon
  n_X <- length(X)
  n_Y <- length(Y)
  
  # Calcul des moyennes d'échantillon
  X_bar <- mean(X)
  Y_bar <- mean(Y)
  
  # Calcul des variances d'échantillon
  S2_X <- var(X)  # R utilise par défaut n-1 pour le calcul de la variance
  S2_Y <- var(Y)
  
  # Calcul des degrés de liberté
  df <- n_X + n_Y - 2
  
  if (equal_var) {
    # Calcul de nu^2 (variance poolée) sous l'hypothèse de variances égales
    nu2 <- ((n_X - 1) * S2_X + (n_Y - 1) * S2_Y) / df
    nu <- sqrt(nu2)
    
    # Calcul de T0
    T0 <- (X_bar - Y_bar) / (nu * sqrt(1/n_X + 1/n_Y))
  } else {
    # Test de Welch pour variances inégales
    var_X <- S2_X / n_X
    var_Y <- S2_Y / n_Y
    T0 <- (X_bar - Y_bar) / sqrt(var_X + var_Y)
    
    # Calcul des degrés de liberté de Welch
    df <- ((var_X + var_Y)^2) / ((var_X^2 / (n_X - 1)) + (var_Y^2 / (n_Y - 1)))
    nu <- NULL  # Non applicable dans ce cas
  }
  
  # Calcul de la valeur critique et de la p-valeur en fonction du type d'alternative
  if (alternative == 'two-sided') {
    critical_value <- qt(1 - alpha/2, df)
    p_value <- 2 * (1 - pt(abs(T0), df))
    reject_H0 <- abs(T0) > critical_value
  } else if (alternative == 'less') {
    critical_value <- qt(alpha, df)
    p_value <- pt(T0, df)
    reject_H0 <- T0 < critical_value
  } else if (alternative == 'greater') {
    critical_value <- qt(1 - alpha, df)
    p_value <- 1 - pt(T0, df)
    reject_H0 <- T0 > critical_value
  } else {
    stop("'alternative' doit être 'two-sided', 'less' ou 'greater'")
  }
  
  return(list(T0 = T0, p_value = p_value, critical_value = critical_value, reject_H0 = reject_H0, nu = nu))
}

# 1.3 Exemples d'application des tests de comparaison de moyennes
#
# Maintenant, appliquons ces tests sur des données simulées.

# Pour la reproductibilité
set.seed(0)

# Exemple 1: Deux populations avec des moyennes différentes et variances connues
mu_X <- 10
mu_Y <- 12
sigma_X <- 3
sigma_Y <- 2.5
n_X <- 30
n_Y <- 25

# Génération des échantillons
X <- rnorm(n_X, mu_X, sigma_X)
Y <- rnorm(n_Y, mu_Y, sigma_Y)

# Affichage des statistiques descriptives
cat("Échantillon X:\n")
cat(sprintf("  Moyenne: %.4f\n", mean(X)))
cat(sprintf("  Écart-type: %.4f\n", sd(X)))
cat(sprintf("  Taille: %d\n", length(X)))
cat("\nÉchantillon Y:\n")
cat(sprintf("  Moyenne: %.4f\n", mean(Y)))
cat(sprintf("  Écart-type: %.4f\n", sd(Y)))
cat(sprintf("  Taille: %d\n", length(Y)))

# Visualisation des distributions
df_combined <- data.frame(
  value = c(X, Y),
  group = c(rep("X", length(X)), rep("Y", length(Y)))
)

ggplot(df_combined, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = mean(X), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = mean(Y), linetype = "dashed", color = "red") +
  annotate("text", x = mean(X) + 0.5, y = 0.05, label = sprintf("Moyenne X: %.2f", mean(X)), color = "blue") +
  annotate("text", x = mean(Y) + 0.5, y = 0.05, label = sprintf("Moyenne Y: %.2f", mean(Y)), color = "red") +
  labs(title = "Distribution des échantillons X et Y", 
       x = "Valeur", 
       y = "Densité") +
  scale_fill_manual(values = c("X" = "blue", "Y" = "red")) +
  theme_minimal()

# Test de comparaison des moyennes avec variances connues
# Utilisons les vraies variances des populations
sigma2_X <- sigma_X^2
sigma2_Y <- sigma_Y^2

# Test bilatéral: H0: μₓ = μᵧ contre H1: μₓ ≠ μᵧ
test_result <- test_moyennes_variances_connues(
  X, Y, sigma2_X, sigma2_Y, alpha=0.05, alternative='two-sided')

cat("Test de comparaison des moyennes avec variances connues (test bilatéral)\n")
cat("H0: μₓ = μᵧ contre H1: μₓ ≠ μᵧ\n")
cat(sprintf("Statistique U0: %.4f\n", test_result$U0))
cat(sprintf("Valeur critique à α=0.05: ±%.4f\n", test_result$critical_value))
cat(sprintf("p-valeur: %.4f\n", test_result$p_value))
cat(sprintf("Décision: %s\n", ifelse(test_result$reject_H0, "Rejeter H0", "Ne pas rejeter H0")))
cat(sprintf("Conclusion: %s\n", ifelse(test_result$reject_H0, 
                                     "Il existe une différence significative entre les moyennes", 
                                     "Il n'y a pas de différence significative entre les moyennes")))

# Visualisation de la distribution de la statistique de test sous H0
x <- seq(-4, 4, length.out = 1000)
y <- dnorm(x, 0, 1)
df_norm <- data.frame(x = x, y = y)

ggplot(df_norm, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  geom_area(data = subset(df_norm, x >= test_result$critical_value | x <= -test_result$critical_value), 
            aes(x = x, y = y), fill = "red", alpha = 0.3) +
  geom_vline(xintercept = test_result$U0, color = "red", linetype = "dashed") +
  annotate("text", x = test_result$U0 + 0.3, y = 0.2, 
           label = sprintf("U0 = %.4f", test_result$U0), color = "red") +
  labs(title = 'Distribution de la statistique de test U0 sous H0',
       x = 'U0',
       y = 'Densité de probabilité') +
  theme_minimal()

# Test unilatéral à gauche: H0: μₓ = μᵧ ou μₓ ≥ μᵧ contre H1: μₓ < μᵧ
test_result <- test_moyennes_variances_connues(
  X, Y, sigma2_X, sigma2_Y, alpha=0.05, alternative='less')

cat("Test de comparaison des moyennes avec variances connues (test unilatéral à gauche)\n")
cat("H0: μₓ = μᵧ ou μₓ ≥ μᵧ contre H1: μₓ < μᵧ\n")
cat(sprintf("Statistique U0: %.4f\n", test_result$U0))
cat(sprintf("Valeur critique à α=0.05: %.4f\n", test_result$critical_value))
cat(sprintf("p-valeur: %.4f\n", test_result$p_value))
cat(sprintf("Décision: %s\n", ifelse(test_result$reject_H0, "Rejeter H0", "Ne pas rejeter H0")))
cat(sprintf("Conclusion: %s\n", ifelse(test_result$reject_H0, 
                                     "La moyenne de X est significativement inférieure à celle de Y", 
                                     "On ne peut pas conclure que la moyenne de X est inférieure à celle de Y")))

# Test unilatéral à droite: H0: μₓ = μᵧ ou μₓ ≤ μᵧ contre H1: μₓ > μᵧ
test_result <- test_moyennes_variances_connues(
  X, Y, sigma2_X, sigma2_Y, alpha=0.05, alternative='greater')

cat("\nTest de comparaison des moyennes avec variances connues (test unilatéral à droite)\n")
cat("H0: μₓ = μᵧ ou μₓ ≤ μᵧ contre H1: μₓ > μᵧ\n")
cat(sprintf("Statistique U0: %.4f\n", test_result$U0))
cat(sprintf("Valeur critique à α=0.05: %.4f\n", test_result$critical_value))
cat(sprintf("p-valeur: %.4f\n", test_result$p_value))
cat(sprintf("Décision: %s\n", ifelse(test_result$reject_H0, "Rejeter H0", "Ne pas rejeter H0")))
cat(sprintf("Conclusion: %s\n", ifelse(test_result$reject_H0, 
                                     "La moyenne de X est significativement supérieure à celle de Y", 
                                     "On ne peut pas conclure que la moyenne de X est supérieure à celle de Y")))

# 1.4 Test avec variances inconnues
#
# Maintenant, supposons que nous ne connaissons pas les variances des populations et utilisons notre deuxième fonction.

# Exemple 2: Deux populations avec des moyennes différentes et variances inconnues
# On utilise les mêmes données que précédemment mais on prétend ne pas connaître les variances

# Test bilatéral: H0: μₓ = μᵧ contre H1: μₓ ≠ μᵧ
test_result <- test_moyennes_variances_inconnues(
  X, Y, alpha=0.05, alternative='two-sided', equal_var=TRUE)

# Degrés de liberté
df <- length(X) + length(Y) - 2

cat("Test de Student pour la comparaison des moyennes avec variances inconnues (test bilatéral)\n")
cat("H0: μₓ = μᵧ contre H1: μₓ ≠ μᵧ\n")
cat(sprintf("Statistique T0: %.4f\n", test_result$T0))
cat(sprintf("Degrés de liberté: %d\n", df))
cat(sprintf("Écart-type poolé (ν): %.4f\n", test_result$nu))
cat(sprintf("Valeur critique à α=0.05: ±%.4f\n", test_result$critical_value))
cat(sprintf("p-valeur: %.4f\n", test_result$p_value))
cat(sprintf("Décision: %s\n", ifelse(test_result$reject_H0, "Rejeter H0", "Ne pas rejeter H0")))
cat(sprintf("Conclusion: %s\n", ifelse(test_result$reject_H0, 
                                     "Il existe une différence significative entre les moyennes", 
                                     "Il n'y a pas de différence significative entre les moyennes")))

# Visualisation de la distribution de la statistique de test sous H0
x <- seq(-4, 4, length.out = 1000)
y <- dt(x, df)
df_t <- data.frame(x = x, y = y)

ggplot(df_t, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  geom_area(data = subset(df_t, x >= test_result$critical_value | x <= -test_result$critical_value), 
            aes(x = x, y = y), fill = "red", alpha = 0.3) +
  geom_vline(xintercept = test_result$T0, color = "red", linetype = "dashed") +
  annotate("text", x = test_result$T0 + 0.3, y = 0.2, 
           label = sprintf("T0 = %.4f", test_result$T0), color = "red") +
  labs(title = sprintf('Distribution de la statistique de test T0 sous H0 (df=%d)', df),
       x = 'T0',
       y = 'Densité de probabilité') +
  theme_minimal()

# Test unilatéral à gauche: H0: μₓ = μᵧ ou μₓ ≥ μᵧ contre H1: μₓ < μᵧ
test_result <- test_moyennes_variances_inconnues(
  X, Y, alpha=0.05, alternative='less', equal_var=TRUE)

cat("Test de Student pour la comparaison des moyennes avec variances inconnues (test unilatéral à gauche)\n")
cat("H0: μₓ = μᵧ ou μₓ ≥ μᵧ contre H1: μₓ < μᵧ\n")
cat(sprintf("Statistique T0: %.4f\n", test_result$T0))
cat(sprintf("Valeur critique à α=0.05: %.4f\n", test_result$critical_value))
cat(sprintf("p-valeur: %.4f\n", test_result$p_value))
cat(sprintf("Décision: %s\n", ifelse(test_result$reject_H0, "Rejeter H0", "Ne pas rejeter H0")))
cat(sprintf("Conclusion: %s\n", ifelse(test_result$reject_H0, 
                                     "La moyenne de X est significativement inférieure à celle de Y", 
                                     "On ne peut pas conclure que la moyenne de X est inférieure à celle de Y")))

# Test unilatéral à droite: H0: μₓ = μᵧ ou μₓ ≤ μᵧ contre H1: μₓ > μᵧ
test_result <- test_moyennes_variances_inconnues(
  X, Y, alpha=0.05, alternative='greater', equal_var=TRUE)

cat("\nTest de Student pour la comparaison des moyennes avec variances inconnues (test unilatéral à droite)\n")
cat("H0: μₓ = μᵧ ou μₓ ≤ μᵧ contre H1: μₓ > μᵧ\n")
cat(sprintf("Statistique T0: %.4f\n", test_result$T0))
cat(sprintf("Valeur critique à α=0.05: %.4f\n", test_result$critical_value))
cat(sprintf("p-valeur: %.4f\n", test_result$p_value))
cat(sprintf("Décision: %s\n", ifelse(test_result$reject_H0, "Rejeter H0", "Ne pas rejeter H0")))
cat(sprintf("Conclusion: %s\n", ifelse(test_result$reject_H0, 
                                     "La moyenne de X est significativement supérieure à celle de Y", 
                                     "On ne peut pas conclure que la moyenne de X est supérieure à celle de Y")))

# Test unilatéral à droite: H0: μₓ = μᵧ ou μₓ ≤ μᵧ contre H1: μₓ > μᵧ
test_result <- test_moyennes_variances_inconnues(
  X, Y, alpha=0.05, alternative='greater', equal_var=FALSE)

cat("\nTest de Student pour la comparaison des moyennes avec variances inconnues (test unilatéral à droite)\n")
cat("H0: μₓ = μᵧ ou μₓ ≤ μᵧ contre H1: μₓ > μᵧ\n")
cat(sprintf("Statistique T0: %.4f\n", test_result$T0))
cat(sprintf("Valeur critique à α=0.05: %.4f\n", test_result$critical_value))
cat(sprintf("p-valeur: %.4f\n", test_result$p_value))
cat(sprintf("Décision: %s\n", ifelse(test_result$reject_H0, "Rejeter H0", "Ne pas rejeter H0")))
cat(sprintf("Conclusion: %s\n", ifelse(test_result$reject_H0, 
                                     "La moyenne de X est significativement supérieure à celle de Y", 
                                     "On ne peut pas conclure que la moyenne de X est supérieure à celle de Y")))

