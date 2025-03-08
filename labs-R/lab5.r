# TP - Théorie des Tests Statistiques
#
# Nous allons explorer deux types de tests statistiques :
# 1. Tests sur une proportion
# 2. Tests d'indépendance

# Importation des bibliothèques nécessaires
library(ggplot2)
library(dplyr)
library(reshape2)

# Configuration pour de meilleurs graphiques
theme_set(theme_light())

# ----------------------------------------------------------------------
# 1. Tests sur une proportion
# ----------------------------------------------------------------------
#
# Ce type de test est similaire à celui d'une moyenne lorsque la variance est connue, 
# étant donné que nous avons une convergence asymptotique vers la loi normale centrée réduite.
#
# La fonction pivotale pour ce test est :
#
# U_0 = (p_hat_n - p_0) / sqrt(p_0 * (1 - p_0) / n) ~ N(0, 1)
#
# où :
# - p_hat_n est la proportion observée dans l'échantillon
# - p_0 est la proportion théorique sous l'hypothèse nulle
# - n est la taille de l'échantillon
#
# Selon les hypothèses, nous avons différentes régions critiques :
#
# | (H_0)                  | (H_1)         | Région Critique                       |
# | ---------------------- | ------------- | ------------------------------------- |
# | p = p_0                | p ≠ p_0       | RC_α = {|U_0| > z_α/2}                |
# | p = p_0 ou p ≥ p_0     | p < p_0       | RC_α = {U_0 < -z_α}                  |
# | p = p_0 ou p ≤ p_0     | p > p_0       | RC_α = {U_0 > z_α}                   |
#
# où z_β est le quantile d'ordre 1 - β de la loi normale centrée réduite.

test_proportion <- function(p_obs, p0, n, alternative = "two.sided", alpha = 0.05) {
  # Réalise un test sur une proportion.
  #
  # Paramètres :
  # - p_obs : proportion observée dans l'échantillon
  # - p0 : proportion théorique sous H0
  # - n : taille de l'échantillon
  # - alternative : "two.sided" (p≠p0), "less" (p<p0), "greater" (p>p0)
  # - alpha : niveau de signification
  #
  # Retourne :
  # - Une liste avec les résultats du test
  
  # Calcul de la statistique de test
  U0 <- (p_obs - p0) / sqrt(p0 * (1 - p0) / n)
  
  # Calcul de la p-valeur selon l'alternative
  if (alternative == "two.sided") {
    p_value <- 2 * (1 - pnorm(abs(U0)))
    critical_value <- qnorm(1 - alpha/2)
    rejet <- abs(U0) > critical_value
  } else if (alternative == "less") {
    p_value <- pnorm(U0)
    critical_value <- qnorm(alpha)
    rejet <- U0 < critical_value
  } else if (alternative == "greater") {
    p_value <- 1 - pnorm(U0)
    critical_value <- qnorm(1 - alpha)
    rejet <- U0 > critical_value
  }
  
  # Conclusion
  conclusion <- ifelse(rejet, "Rejet de H0", "Non-rejet de H0")
  
  return(list(
    U0 = U0,
    p_value = p_value,
    critical_value = critical_value,
    conclusion = conclusion,
    rejet = rejet
  ))
}

# ----------------------------------------------------------------------
# Exemple : Étudier à l'étranger
# ----------------------------------------------------------------------
#
# D'après l'exemple du support de cours :
#
# Le Ministre de l'enseignement supérieur en Mauritanie affirme qu'un bachelier 
# sur trois désire faire ses études à l'étranger. À la suite d'un sondage auprès 
# de 1000 nouveaux bacheliers, 280 désirent vouloir poursuivre leurs études à l'étranger. 
# Le ministre s'est-il trompé au risque de 10% ? 5% ? 1% ?
#
# Nous allons tester si la proportion réelle diffère de 1/3 (ce qui correspond à H1 : p ≠ p0).

# Données du problème
p0 <- 1/3  # Proportion théorique selon le ministre
n <- 1000  # Taille de l'échantillon
succes <- 280  # Nombre de bacheliers souhaitant étudier à l'étranger
p_obs <- succes / n  # Proportion observée

cat(sprintf("Proportion observée : %.4f\n", p_obs))
cat(sprintf("Proportion théorique : %.4f\n", p0))
cat(sprintf("Écart : %.4f\n", abs(p_obs - p0)))

# Test pour α = 10%
resultat_10 <- test_proportion(p_obs, p0, n, alternative = "two.sided", alpha = 0.10)
cat("Résultat du test pour α = 10% :\n")
for (nom in names(resultat_10)) {
  if (is.numeric(resultat_10[[nom]])) {
    cat(sprintf("%s : %.4f\n", nom, resultat_10[[nom]]))
  } else {
    cat(sprintf("%s : %s\n", nom, resultat_10[[nom]]))
  }
}

# Test pour α = 5%
resultat_5 <- test_proportion(p_obs, p0, n, alternative = "two.sided", alpha = 0.05)
cat("\nRésultat du test pour α = 5% :\n")
for (nom in names(resultat_5)) {
  if (is.numeric(resultat_5[[nom]])) {
    cat(sprintf("%s : %.4f\n", nom, resultat_5[[nom]]))
  } else {
    cat(sprintf("%s : %s\n", nom, resultat_5[[nom]]))
  }
}

# Test pour α = 1%
resultat_1 <- test_proportion(p_obs, p0, n, alternative = "two.sided", alpha = 0.01)
cat("\nRésultat du test pour α = 1% :\n")
for (nom in names(resultat_1)) {
  if (is.numeric(resultat_1[[nom]])) {
    cat(sprintf("%s : %.4f\n", nom, resultat_1[[nom]]))
  } else {
    cat(sprintf("%s : %s\n", nom, resultat_1[[nom]]))
  }
}

# Visualisation de la distribution normale et des zones critiques
visualiser_test_proportion <- function(U0, alpha_values = c(0.1, 0.05, 0.01)) {
  x <- seq(-4, 4, length.out = 1000)
  y <- dnorm(x)
  
  # Créer un data.frame pour ggplot
  df <- data.frame(x = x, y = y)
  
  # Créer le graphique de base
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(size = 1, color = "blue") +
    geom_vline(xintercept = U0, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = "Test de proportion: distribution de la statistique U0 sous H0",
      x = "U0",
      y = "Densité de probabilité"
    )
  
  # Ajouter les zones critiques pour chaque alpha
  colors <- c("lightblue", "lightgreen", "salmon")
  alphas <- c(0.5, 0.3, 0.2)
  
  for (i in seq_along(alpha_values)) {
    alpha <- alpha_values[i]
    critical <- qnorm(1 - alpha/2)
    
    # Ajouter zone critique à gauche
    p <- p + 
      geom_area(
        data = subset(df, x <= -critical),
        aes(x = x, y = y),
        fill = colors[i],
        alpha = alphas[i]
      )
    
    # Ajouter zone critique à droite
    p <- p + 
      geom_area(
        data = subset(df, x >= critical),
        aes(x = x, y = y),
        fill = colors[i],
        alpha = alphas[i]
      )
    
    # Ajouter lignes critiques
    p <- p + 
      geom_vline(xintercept = -critical, color = c("blue", "green", "red")[i], 
                linetype = "dashed", alpha = 0.7) +
      geom_vline(xintercept = critical, color = c("blue", "green", "red")[i], 
                linetype = "dashed", alpha = 0.7)
  }
  
  # Ajouter légende
  p <- p + 
    annotate("text", x = 3, y = 0.3, 
            label = sprintf("U0 = %.4f", U0), color = "red") +
    annotate("text", x = 3, y = 0.27, 
            label = sprintf("Seuil 10%% = %.4f", qnorm(1 - 0.1/2)), color = "blue") +
    annotate("text", x = 3, y = 0.24, 
            label = sprintf("Seuil 5%% = %.4f", qnorm(1 - 0.05/2)), color = "green") +
    annotate("text", x = 3, y = 0.21, 
            label = sprintf("Seuil 1%% = %.4f", qnorm(1 - 0.01/2)), color = "red")
  
  return(p)
}

# Pour visualiser le test
# print(visualiser_test_proportion(resultat_10$U0))

# ----------------------------------------------------------------------
# Interprétation
# ----------------------------------------------------------------------
#
# D'après les résultats obtenus, nous pouvons conclure :
#
# 1. Pour un risque α = 10%, on rejette H0. La p-valeur est inférieure à 0.10, 
#    ce qui signifie que l'écart observé est statistiquement significatif à ce niveau. 
#    Le ministre s'est donc trompé avec un risque de 10%.
#
# 2. Pour un risque α = 5%, on rejette H0. La p-valeur est inférieure à 0.05, 
#    ce qui signifie que l'écart observé est statistiquement significatif à ce niveau. 
#    Le ministre s'est donc trompé avec un risque de 5%.
#
# 3. Pour un risque α = 1%, on rejette H0. La p-valeur est inférieure à 0.01, 
#    ce qui signifie que l'écart observé est statistiquement significatif à ce niveau. 
#    Le ministre s'est donc trompé avec un risque de 1%.
#
# La statistique de test nous indique que la proportion observée (28%) est significativement 
# différente de la proportion affirmée par le ministre (33.33%).

# ----------------------------------------------------------------------
# 2. Tests d'indépendance
# ----------------------------------------------------------------------
#
# Ce type de test permet de déterminer s'il existe une relation entre deux 
# variables qualitatives ou quantitatives regroupées en classes.
#
# Soient X une variable à k modalités et Y une variable à p modalités. Nous voulons tester :
# - (H0) : X indépendant de Y 
# - contre (H1) : X dépendant de Y.
#
# Sous H0, nous avons : 
# - fX|Y = fX et fY|X = fY 
# - P(X = i, Y = j) = P(X = i)P(Y = j) pour tout i, j
#
# La statistique de test est définie par :
#
# δ² = ∑∑(nij - nij*)²/nij*
#
# où :
# - nij est l'effectif observé de la cellule (i,j)
# - nij* est l'effectif théorique sous H0, calculé comme (ni. × n.j) / n..
#
# Sous H0, cette statistique suit asymptotiquement une loi du χ² à (k-1)×(p-1) degrés de liberté.
#
# On rejette H0 lorsque δ² > z(k-1)×(p-1),α, où z(k-1)×(p-1),α est le quantile 
# d'ordre 1-α de la loi du χ² à (k-1)×(p-1) degrés de liberté.
#
# Note : Ce test est valide si moins de 20% des effectifs théoriques sont inférieurs à 5.

test_independance <- function(tableau_contingence, alpha = 0.05) {
  # Réalise un test d'indépendance du χ² à partir d'un tableau de contingence.
  #
  # Paramètres :
  # - tableau_contingence : matrice ou data.frame
  # - alpha : niveau de signification
  #
  # Retourne :
  # - Une liste avec les résultats du test
  
  # Conversion en matrice pour les calculs
  if (!is.matrix(tableau_contingence)) {
    tableau_contingence <- as.matrix(tableau_contingence)
  }
  
  # Calcul des marges
  n <- sum(tableau_contingence)
  marges_lignes <- rowSums(tableau_contingence)
  marges_colonnes <- colSums(tableau_contingence)
  
  # Calcul des effectifs théoriques
  tableau_theorique <- matrix(0, nrow = nrow(tableau_contingence), ncol = ncol(tableau_contingence))
  for (i in 1:nrow(tableau_contingence)) {
    for (j in 1:ncol(tableau_contingence)) {
      tableau_theorique[i, j] <- marges_lignes[i] * marges_colonnes[j] / n
    }
  }
  
  # Calcul de la statistique du χ²
  chi2_stat <- 0
  for (i in 1:nrow(tableau_contingence)) {
    for (j in 1:ncol(tableau_contingence)) {
      nij <- tableau_contingence[i, j]
      nij_star <- tableau_theorique[i, j]
      chi2_stat <- chi2_stat + (nij - nij_star)^2 / nij_star
    }
  }
  
  # Calcul des degrés de liberté
  k <- nrow(tableau_contingence)  # nb de lignes
  p <- ncol(tableau_contingence)  # nb de colonnes
  ddl <- (k - 1) * (p - 1)
  
  # Calcul de la p-valeur
  p_value <- 1 - pchisq(chi2_stat, ddl)
  
  # Valeur critique
  critical_value <- qchisq(1 - alpha, ddl)
  
  # Conclusion
  rejet <- chi2_stat > critical_value
  conclusion <- ifelse(rejet, "Rejet de H0 (dépendance)", "Non-rejet de H0 (indépendance)")
  
  # Vérification de la validité du test
  nb_effectifs_faibles <- sum(tableau_theorique < 5)
  pct_effectifs_faibles <- nb_effectifs_faibles / (k * p) * 100
  validite <- pct_effectifs_faibles < 20
  
  # Conversion en data.frame pour le retour
  tableau_theorique_df <- as.data.frame(tableau_theorique)
  colnames(tableau_theorique_df) <- colnames(tableau_contingence)
  rownames(tableau_theorique_df) <- rownames(tableau_contingence)
  
  return(list(
    chi2_stat = chi2_stat,
    p_value = p_value,
    ddl = ddl,
    critical_value = critical_value,
    conclusion = conclusion,
    tableau_theorique = tableau_theorique_df,
    validite = validite,
    pct_effectifs_faibles = pct_effectifs_faibles,
    rejet = rejet
  ))
}

# ----------------------------------------------------------------------
# Exemple : Taux de présence et notes
# ----------------------------------------------------------------------
#
# D'après l'exemple du support de cours, nous voulons tester si le taux de présence en cours influence les notes.
#
# Voici le tableau de contingence :
#
# | Notes        | Moins de 50% de présence | Plus de 50% de présence | Total |
# | ------------ | ------------------------ | ----------------------- | ----- |
# | 0 à 5 sur 20 | 0                        | 8                       | 8     |
# | 5 à 10 sur 20| 20                       | 12                      | 32    |
# | 10 à 15 sur 20| 12                      | 44                      | 56    |
# | 15 à 20 sur 20| 0                       | 24                      | 24    |
# | Total        | 32                       | 88                      | 120   |

# Création du tableau de contingence
donnees <- matrix(
  c(0, 8, 20, 12, 12, 44, 0, 24),
  nrow = 4,
  byrow = TRUE
)
rownames(donnees) <- c("0 à 5 sur 20", "5 à 10 sur 20", "10 à 15 sur 20", "15 à 20 sur 20")
colnames(donnees) <- c("Moins de 50% de présence", "Plus de 50% de présence")
tableau <- as.data.frame(donnees)

# Affichage du tableau
cat("Tableau de contingence :\n")
print(tableau)

# Ajout des totaux
tableau$Total <- rowSums(tableau)
tableau["Total", ] <- colSums(tableau)
cat("\nTableau avec totaux :\n")
print(tableau)

# Suppression des totaux pour l'analyse
tableau_analyse <- tableau[1:4, 1:2]

# Réalisation du test d'indépendance
resultat <- test_independance(tableau_analyse, alpha = 0.05)

# Affichage des résultats
cat("Résultats du test d'indépendance :\n")
for (nom in names(resultat)) {
  if (nom != "tableau_theorique") {
    if (is.numeric(resultat[[nom]])) {
      cat(sprintf("%s : %.4f\n", nom, resultat[[nom]]))
    } else {
      cat(sprintf("%s : %s\n", nom, resultat[[nom]]))
    }
  }
}

cat("\nTableau des effectifs théoriques :\n")
print(round(resultat$tableau_theorique, 2))

# Visualisation du tableau de contingence
visualiser_tableau_contingence <- function(tableau) {
  # Convertir le tableau en format long pour ggplot
  donnees_long <- melt(as.matrix(tableau))
  colnames(donnees_long) <- c("Notes", "Présence", "Effectif")
  
  # Créer le heatmap
  ggplot(donnees_long, aes(x = Présence, y = Notes, fill = Effectif)) +
    geom_tile() +
    geom_text(aes(label = Effectif), color = "black") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(
      title = "Tableau de contingence: Taux de présence vs Notes",
      x = "Taux de présence",
      y = "Notes",
      fill = "Effectifs observés"
    ) +
    theme_minimal()
}

# Pour visualiser le tableau de contingence
# print(visualiser_tableau_contingence(tableau_analyse))

# Visualisation de la distribution du Chi-2 et de la zone critique
visualiser_test_chi2 <- function(chi2_stat, ddl, critical_value) {
  x <- seq(0, max(20, chi2_stat * 1.5), length.out = 1000)
  y <- dchisq(x, ddl)
  
  # Créer un data.frame pour ggplot
  df <- data.frame(x = x, y = y)
  
  # Créer le graphique
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(size = 1, color = "blue") +
    geom_vline(xintercept = chi2_stat, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = "Test d'indépendance: distribution du χ² sous H0",
      x = "χ²",
      y = "Densité de probabilité"
    )
  
  # Ajouter la zone critique
  p <- p + 
    geom_area(
      data = subset(df, x >= critical_value),
      aes(x = x, y = y),
      fill = "lightcoral",
      alpha = 0.5
    )
  
  # Ajouter la ligne critique
  p <- p + 
    geom_vline(xintercept = critical_value, color = "blue", linetype = "dashed", size = 1)
  
  # Ajouter légende
  p <- p + 
    annotate("text", x = chi2_stat * 1.2, y = max(y) * 0.8, 
            label = sprintf("χ² calculé = %.4f", chi2_stat), color = "red") +
    annotate("text", x = critical_value * 1.2, y = max(y) * 0.7, 
            label = sprintf("Valeur critique = %.4f", critical_value), color = "blue")
  
  return(p)
}

# Pour visualiser le test du chi-2
# print(visualiser_test_chi2(resultat$chi2_stat, resultat$ddl, resultat$critical_value))

# ----------------------------------------------------------------------
# Rendu de TP
# ----------------------------------------------------------------------
#
# L'objectif du rendu de TP est de réaliser un test où l'on ne peut pas rejeter l'hypothèse nulle 
# dans le cas des tests sur une proportion.
# Pour cela, vous devez réaliser un test sur une proportion avec les paramètres suivants :
#
# - p_0 = 0.5 (proportion sous l'hypothèse nulle)
# - n = 100 (taille de l'échantillon)
# - p = 0.5 (proportion réelle)
# - α = 0.05 (seuil de signification)
#
#
# 1. Générer un échantillon de taille 100 avec une proportion de 0.5
# 2. Calculer la statistique de test
# 3. Calculer la p-valeur
# 4. Faire des visualisations pour montrer que l'hypothèse nulle n'est pas rejetée
# 5. Conclure

