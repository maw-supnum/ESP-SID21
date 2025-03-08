# TP sur la Théorie des Tests Statistiques
#
# Dans ce TP, nous allons explorer les tests paraméétriques sur des lois normales, en particulier :
# 1. Tests sur la moyenne d'une loi normale (avec variance connue ou inconnue)
# 2. Tests sur la variance d'une loi normale (avec moyenne connue)
#
# Ces tests sont fondamentaux en statistique inférentielle et constituent la base de nombreuses méthodes d'analyse statistique.

# Importation des bibliothèques nécessaires
library(ggplot2)

# Configuration pour de meilleurs graphiques
theme_set(theme_bw())

## 1. Tests sur la moyenne d'une loi normale
#
# Considérons l'hypothèse suivante :
#
# On suppose que l'on dispose de (X_1, X_2, ..., X_n) un échantillon issu de X qui suit une loi N(μ, σ²).
#
# Nous allons examiner les tests pour H_0 : μ = μ_0 contre différentes alternatives.

### 1.1 Génération d'un échantillon pour la démonstration

# Paramètres pour notre simulation
set.seed(0)  # Pour la reproductibilité
n <- 400  # Taille de l'échantillon
mu_true <- 5  # Vraie moyenne de la population
sigma_true <- 2  # Vraie écart-type de la population

# Génération d'un échantillon
sample_data <- rnorm(n, mu_true, sigma_true)

# Statistiques descriptives
x_bar <- mean(sample_data)
s <- sd(sample_data)

cat(sprintf("Moyenne de l'échantillon: %.4f\n", x_bar))
cat(sprintf("Écart-type de l'échantillon: %.4f\n", s))

# Visualisation de l'échantillon
p <- ggplot() +
  geom_histogram(aes(x = sample_data, y = after_stat(density)), bins = 15, alpha = 0.7, fill = "blue") +
  geom_density(aes(x = sample_data), color = "blue", linewidth = 1) +
  stat_function(fun = dnorm, args = list(mean = mu_true, sd = sigma_true), 
                color = "red", linewidth = 1) +
  stat_function(fun = dnorm, args = list(mean = x_bar, sd = s), 
                color = "green", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = mu_true, color = "red", alpha = 0.5) +
  geom_vline(xintercept = x_bar, color = "green", alpha = 0.5, linetype = "dashed") +
  labs(
    title = "Distribution de l'échantillon vs. Distribution théorique",
    x = "Valeur",
    y = "Densité"
  ) +
  theme_minimal()

print(p)

### 1.2 Test sur la moyenne avec variance connue
#
# Si σ² est connu, on utilise la fonction pivotale Z_0 = (X̄ - μ_0)/(σ/√n) ~ N(0, 1)
#
# Cas bilatéral: H_0 : μ = μ_0 contre H_1 : μ ≠ μ_0
# - Région critique: RC_α = {|Z_0| > z_{α/2}}
#
# Cas unilatéral à gauche: H_0 : μ = μ_0 (ou μ ≥ μ_0) contre H_1 : μ < μ_0
# - Région critique: RC_α = {Z_0 < -z_{α}}
#
# Cas unilatéral à droite: H_0 : μ = μ_0 (ou μ ≤ μ_0) contre H_1 : μ > μ_0
# - Région critique: RC_α = {Z_0 > z_{α}}

test_mean_known_variance <- function(sample, mu0, sigma_known, alpha = 0.05, alternative = "two.sided") {
  # Test sur la moyenne d'une loi normale avec variance connue.
  #
  # Parameters:
  # -----------
  # sample : vector
  #     Échantillon à tester
  # mu0 : numeric
  #     Valeur de la moyenne sous H0
  # sigma_known : numeric
  #     Écart-type connu de la population
  # alpha : numeric, optional
  #     Niveau de signification (default = 0.05)
  # alternative : character, optional
  #     Type d'alternative ("two.sided", "less", "greater", default = "two.sided")
  #
  # Returns:
  # --------
  # list
  #     Résultats du test (statistique, p-value, décision)
  
  n <- length(sample)
  x_bar <- mean(sample)
  z0 <- (x_bar - mu0) / (sigma_known / sqrt(n))
  
  # Calcul de la p-value selon l'alternative
  if (alternative == "two.sided") {
    p_value <- 2 * pnorm(-abs(z0))
    critical_value <- qnorm(1 - alpha/2)
    reject <- abs(z0) > critical_value
  } else if (alternative == "less") {
    p_value <- pnorm(z0)
    critical_value <- qnorm(alpha)
    reject <- z0 < critical_value
  } else if (alternative == "greater") {
    p_value <- pnorm(z0, lower.tail = FALSE)
    critical_value <- qnorm(1 - alpha)
    reject <- z0 > critical_value
  } else {
    stop("alternative must be 'two.sided', 'less' or 'greater'")
  }
  
  # Décision
  decision <- ifelse(reject, "Rejeter H0", "Ne pas rejeter H0")
  
  return(list(
    statistic = z0,
    p_value = p_value,
    critical_value = critical_value,
    decision = decision,
    reject = reject
  ))
}

# Test avec variance connue
mu0 <- 4.5  # Valeur hypothétique de la moyenne
sigma_known <- 2  # Variance connue (dans un cas réel, cette valeur serait fournie)
alpha <- 0.05  # Niveau de signification

# Test bilatéral
result_bilateral <- test_mean_known_variance(sample_data, mu0, sigma_known, alpha, "two.sided")
cat("\nTest bilatéral H0: μ = μ0 vs H1: μ ≠ μ0\n")
cat(sprintf("Statistique Z0 = %.4f\n", result_bilateral$statistic))
cat(sprintf("Valeur critique z%s = %.4f\n", alpha/2, result_bilateral$critical_value))
cat(sprintf("p-value = %.4f\n", result_bilateral$p_value))
cat(sprintf("Décision (alpha=%s): %s\n", alpha, result_bilateral$decision))

# Test unilatéral à gauche
result_left <- test_mean_known_variance(sample_data, mu0, sigma_known, alpha, "less")
cat("\nTest unilatéral à gauche H0: μ = μ0 (ou μ ≥ μ0) vs H1: μ < μ0\n")
cat(sprintf("Statistique Z0 = %.4f\n", result_left$statistic))
cat(sprintf("Valeur critique -z%s = %.4f\n", alpha, result_left$critical_value))
cat(sprintf("p-value = %.4f\n", result_left$p_value))
cat(sprintf("Décision (alpha=%s): %s\n", alpha, result_left$decision))

# Test unilatéral à droite
result_right <- test_mean_known_variance(sample_data, mu0, sigma_known, alpha, "greater")
cat("\nTest unilatéral à droite H0: μ = μ0 (ou μ ≤ μ0) vs H1: μ > μ0\n")
cat(sprintf("Statistique Z0 = %.4f\n", result_right$statistic))
cat(sprintf("Valeur critique z%s = %.4f\n", alpha, result_right$critical_value))
cat(sprintf("p-value = %.4f\n", result_right$p_value))
cat(sprintf("Décision (alpha=%s): %s\n", alpha, result_right$decision))

### Visualisation du test (variance connue)

plot_test_normal <- function(result, test_type = "two.sided", alpha = 0.05) {
  # Visualise la région critique et la statistique de test pour une loi normale.
  #
  # Parameters:
  # -----------
  # result : list
  #     Résultats du test
  # test_type : character, optional
  #     Type de test (default = "two.sided")
  # alpha : numeric, optional
  #     Niveau de signification (default = 0.05)
  
  z0 <- result$statistic
  z <- seq(-4, 4, length.out = 1000)
  pdf_z <- dnorm(z)
  
  p <- ggplot() +
    geom_line(aes(x = z, y = pdf_z), color = "blue", linewidth = 1) +
    labs(x = "z", y = "Densité")
  
  if (test_type == "two.sided") {
    # Région critique bilatérale
    z_crit <- qnorm(1 - alpha/2)
    idx_left <- z <= -z_crit
    idx_right <- z >= z_crit
    
    p <- p +
      geom_area(aes(x = z[idx_left], y = pdf_z[idx_left]), fill = "red", alpha = 0.5) +
      geom_area(aes(x = z[idx_right], y = pdf_z[idx_right]), fill = "red", alpha = 0.5) +
      geom_vline(xintercept = -z_crit, color = "red", linetype = "dashed", alpha = 0.7) +
      geom_vline(xintercept = z_crit, color = "red", linetype = "dashed", alpha = 0.7) +
      labs(title = "Test bilatéral sur la moyenne (H0: μ = μ0 vs H1: μ ≠ μ0)")
    
  } else if (test_type == "less") {
    # Région critique unilatérale à gauche
    z_crit <- qnorm(alpha)
    idx <- z <= z_crit
    
    p <- p +
      geom_area(aes(x = z[idx], y = pdf_z[idx]), fill = "red", alpha = 0.5) +
      geom_vline(xintercept = z_crit, color = "red", linetype = "dashed", alpha = 0.7) +
      labs(title = "Test unilatéral à gauche (H0: μ = μ0 ou μ ≥ μ0 vs H1: μ < μ0)")
    
  } else if (test_type == "greater") {
    # Région critique unilatérale à droite
    z_crit <- qnorm(1 - alpha)
    idx <- z >= z_crit
    
    p <- p +
      geom_area(aes(x = z[idx], y = pdf_z[idx]), fill = "red", alpha = 0.5) +
      geom_vline(xintercept = z_crit, color = "red", linetype = "dashed", alpha = 0.7) +
      labs(title = "Test unilatéral à droite (H0: μ = μ0 ou μ ≤ μ0 vs H1: μ > μ0)")
  }
  
  # Statistique calculée
  p <- p + geom_vline(xintercept = z0, color = "green", size = 1)
  
  # Décision
  decision_text <- ifelse(result$reject, "Décision: Rejeter H0", "Décision: Ne pas rejeter H0")
  decision_color <- ifelse(result$reject, "red", "green")
  
  p <- p + annotate("text", x = -3, y = 0.35, label = decision_text, 
                     hjust = 0, fontface = 2, color = decision_color)
  
  # Légende
  p <- p + annotate("text", x = -3, y = 0.32, label = sprintf("Z0 = %.4f", z0), 
                     hjust = 0, color = "green")
  
  if (test_type == "two.sided") {
    p <- p + annotate("text", x = -3, y = 0.29, label = sprintf("±z.crit = ±%.4f", z_crit), 
                       hjust = 0, color = "red")
  } else if (test_type == "less") {
    p <- p + annotate("text", x = -3, y = 0.29, label = sprintf("z.crit = %.4f", z_crit), 
                       hjust = 0, color = "red")
  } else if (test_type == "greater") {
    p <- p + annotate("text", x = -3, y = 0.29, label = sprintf("z.crit = %.4f", z_crit), 
                       hjust = 0, color = "red")
  }
  
  p <- p + theme_minimal()
  
  print(p)
}

# Visualisation des tests pour variance connue
plot_test_normal(result_bilateral, "two.sided", alpha)
plot_test_normal(result_left, "less", alpha)
plot_test_normal(result_right, "greater", alpha)

### 1.3 Test sur la moyenne avec variance inconnue
#
# Si σ² est inconnu, on utilise la fonction pivotale T_0 = (X̄ - μ_0)/(S/√n) ~ T(n-1)
#
# Cas bilatéral: H_0 : μ = μ_0 contre H_1 : μ ≠ μ_0
# - Région critique: RC_α = {|T_0| > t_{n-1,α/2}}
#
# Cas unilatéral à gauche: H_0 : μ = μ_0 (ou μ ≥ μ_0) contre H_1 : μ < μ_0
# - Région critique: RC_α = {T_0 < -t_{n-1,α}}
#
# Cas unilatéral à droite: H_0 : μ = μ_0 (ou μ ≤ μ_0) contre H_1 : μ > μ_0
# - Région critique: RC_α = {T_0 > t_{n-1,α}}

test_mean_unknown_variance <- function(sample, mu0, alpha = 0.05, alternative = "two.sided") {
  # Test sur la moyenne d'une loi normale avec variance inconnue.
  #
  # Parameters:
  # -----------
  # sample : vector
  #     Échantillon à tester
  # mu0 : numeric
  #     Valeur de la moyenne sous H0
  # alpha : numeric, optional
  #     Niveau de signification (default = 0.05)
  # alternative : character, optional
  #     Type d'alternative ("two.sided", "less", "greater", default = "two.sided")
  #
  # Returns:
  # --------
  # list
  #     Résultats du test (statistique, p-value, décision)
  
  n <- length(sample)
  x_bar <- mean(sample)
  s <- sd(sample)
  t0 <- (x_bar - mu0) / (s / sqrt(n))
  
  # Degrés de liberté
  df <- n - 1
  
  # Calcul de la p-value selon l'alternative
  if (alternative == "two.sided") {
    p_value <- 2 * pt(-abs(t0), df)
    critical_value <- qt(1 - alpha/2, df)
    reject <- abs(t0) > critical_value
  } else if (alternative == "less") {
    p_value <- pt(t0, df)
    critical_value <- qt(alpha, df)
    reject <- t0 < critical_value
  } else if (alternative == "greater") {
    p_value <- pt(t0, df, lower.tail = FALSE)
    critical_value <- qt(1 - alpha, df)
    reject <- t0 > critical_value
  } else {
    stop("alternative must be 'two.sided', 'less' or 'greater'")
  }
  
  # Décision
  decision <- ifelse(reject, "Rejeter H0", "Ne pas rejeter H0")
  
  return(list(
    statistic = t0,
    p_value = p_value,
    critical_value = critical_value,
    decision = decision,
    reject = reject,
    df = df
  ))
}

# Test avec variance inconnue
mu0 <- 4.5  # Valeur hypothétique de la moyenne
alpha <- 0.05  # Niveau de signification

# Test bilatéral
result_bilateral_t <- test_mean_unknown_variance(sample_data, mu0, alpha, "two.sided")
cat("\nTest bilatéral H0: μ = μ0 vs H1: μ ≠ μ0\n")
cat(sprintf("Statistique T0 = %.4f\n", result_bilateral_t$statistic))
cat(sprintf("Degrés de liberté = %d\n", result_bilateral_t$df))
cat(sprintf("Valeur critique t_%d,%s = %.4f\n", result_bilateral_t$df, alpha/2, result_bilateral_t$critical_value))
cat(sprintf("p-value = %.4f\n", result_bilateral_t$p_value))
cat(sprintf("Décision (alpha=%s): %s\n", alpha, result_bilateral_t$decision))

# Test unilatéral à gauche
result_left_t <- test_mean_unknown_variance(sample_data, mu0, alpha, "less")
cat("\nTest unilatéral à gauche H0: μ = μ0 (ou μ ≥ μ0) vs H1: μ < μ0\n")
cat(sprintf("Statistique T0 = %.4f\n", result_left_t$statistic))
cat(sprintf("Valeur critique -t_%d,%s = %.4f\n", result_left_t$df, alpha, result_left_t$critical_value))
cat(sprintf("p-value = %.4f\n", result_left_t$p_value))
cat(sprintf("Décision (alpha=%s): %s\n", alpha, result_left_t$decision))

# Test unilatéral à droite
result_right_t <- test_mean_unknown_variance(sample_data, mu0, alpha, "greater")
cat("\nTest unilatéral à droite H0: μ = μ0 (ou μ ≤ μ0) vs H1: μ > μ0\n")
cat(sprintf("Statistique T0 = %.4f\n", result_right_t$statistic))
cat(sprintf("Valeur critique t_%d,%s = %.4f\n", result_right_t$df, alpha, result_right_t$critical_value))
cat(sprintf("p-value = %.4f\n", result_right_t$p_value))
cat(sprintf("Décision (alpha=%s): %s\n", alpha, result_right_t$decision))

### Visualisation du test (variance inconnue)

plot_test_student <- function(result, test_type = "two.sided", alpha = 0.05) {
  # Visualise la région critique et la statistique de test pour une loi de Student.
  #
  # Parameters:
  # -----------
  # result : list
  #     Résultats du test
  # test_type : character, optional
  #     Type de test (default = "two.sided")
  # alpha : numeric, optional
  #     Niveau de signification (default = 0.05)
  
  t0 <- result$statistic
  df <- result$df
  x <- seq(-4, 4, length.out = 1000)
  pdf_t <- dt(x, df)
  
  p <- ggplot() +
    geom_line(aes(x = x, y = pdf_t), color = "blue", linewidth = 1) +
    labs(x = "t", y = "Densité")
  
  if (test_type == "two.sided") {
    # Région critique bilatérale
    t_crit <- qt(1 - alpha/2, df)
    idx_left <- x <= -t_crit
    idx_right <- x >= t_crit
    
    p <- p +
      geom_area(aes(x = x[idx_left], y = pdf_t[idx_left]), fill = "red", alpha = 0.5) +
      geom_area(aes(x = x[idx_right], y = pdf_t[idx_right]), fill = "red", alpha = 0.5) +
      geom_vline(xintercept = -t_crit, color = "red", linetype = "dashed", alpha = 0.7) +
      geom_vline(xintercept = t_crit, color = "red", linetype = "dashed", alpha = 0.7) +
      labs(title = "Test bilatéral sur la moyenne avec variance inconnue (H0: μ = μ0 vs H1: μ ≠ μ0)")
    
  } else if (test_type == "less") {
    # Région critique unilatérale à gauche
    t_crit <- qt(alpha, df)
    idx <- x <= t_crit
    
    p <- p +
      geom_area(aes(x = x[idx], y = pdf_t[idx]), fill = "red", alpha = 0.5) +
      geom_vline(xintercept = t_crit, color = "red", linetype = "dashed", alpha = 0.7) +
      labs(title = "Test unilatéral à gauche avec variance inconnue (H0: μ = μ0 ou μ ≥ μ0 vs H1: μ < μ0)")
    
  } else if (test_type == "greater") {
    # Région critique unilatérale à droite
    t_crit <- qt(1 - alpha, df)
    idx <- x >= t_crit
    
    p <- p +
      geom_area(aes(x = x[idx], y = pdf_t[idx]), fill = "red", alpha = 0.5) +
      geom_vline(xintercept = t_crit, color = "red", linetype = "dashed", alpha = 0.7) +
      labs(title = "Test unilatéral à droite avec variance inconnue (H0: μ = μ0 ou μ ≤ μ0 vs H1: μ > μ0)")
  }
  
  # Statistique calculée
  p <- p + geom_vline(xintercept = t0, color = "green", size = 1)
  
  # Décision
  decision_text <- ifelse(result$reject, "Décision: Rejeter H0", "Décision: Ne pas rejeter H0")
  decision_color <- ifelse(result$reject, "red", "green")
  
  p <- p + annotate("text", x = -3, y = 0.35, label = decision_text, 
                     hjust = 0, fontface = 2, color = decision_color)
  
  # Légende
  p <- p + annotate("text", x = -3, y = 0.32, label = sprintf("T0 = %.4f", t0), 
                     hjust = 0, color = "green")
  
  if (test_type == "two.sided") {
    p <- p + annotate("text", x = -3, y = 0.29, label = sprintf("±t.crit = ±%.4f", t_crit), 
                       hjust = 0, color = "red")
  } else if (test_type == "less") {
    p <- p + annotate("text", x = -3, y = 0.29, label = sprintf("t.crit = %.4f", t_crit), 
                       hjust = 0, color = "red")
  } else if (test_type == "greater") {
    p <- p + annotate("text", x = -3, y = 0.29, label = sprintf("t.crit = %.4f", t_crit), 
                       hjust = 0, color = "red")
  }
  
  p <- p + theme_minimal()
  
  print(p)
}

# Visualisation des tests pour variance inconnue
plot_test_student(result_bilateral_t, "two.sided", alpha)
plot_test_student(result_left_t, "less", alpha)
plot_test_student(result_right_t, "greater", alpha)

# Rendu de TP
#
# Maintenant on effectue un test sur la variance d'une loi normale avec moyenne connue.
# L'objectif sera de montrer qu'on ne peut pas rejeter l'hypothèse nulle pour un échantillon 
# généré à partir d'une loi normale de paramètres μ = 0, σ = 1 et α = 0.05. On fixe σ_0 = 2.
#
# Pour cela, vous devrez :
# 1. Générer un échantillon de taille n = 400 à partir d'une loi normale de paramètres μ = 0 et σ = 1
# 2. Réaliser les tests sur la variance avec σ_0 = 2 et α = 0.05
# 3. Visualiser les tests réalisés

