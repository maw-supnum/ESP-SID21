# TP Estimation de Proportions avec Dataset Digits
#
# Dans ce TP, nous allons explorer l'estimation par intervalle de confiance pour une proportion

# Importation des bibliothèques nécessaires
library(ggplot2)
library(dplyr)
library(mlbench)  # Pour charger des datasets
library(caret)    # Pour le split

# Configuration pour la reproductibilité
set.seed(0)

# ------------------------------------------------------------------------------
# 1. Chargement et exploration des données
# ------------------------------------------------------------------------------
# Nous allons utiliser le dataset 'LetterRecognition' de mlbench et le transformer
# pour simuler un dataset digits (en gardant seulement les 10 premières classes)

# Charger le dataset
data(LetterRecognition)
# Nous allons transformer ce dataset pour simuler un dataset digits
# Gardons seulement les 10 premières lettres (A-J) et les transformer en chiffres (0-9)
digits_data <- LetterRecognition[LetterRecognition$lettr %in% LETTERS[1:10],]
# Conversion de A-J à 0-9
digits_data$target <- as.integer(as.factor(digits_data$lettr)) - 1

# Extraction des cibles
y <- digits_data$target
X <- as.matrix(digits_data[,2:17])  # Caractéristiques

# Affichage de quelques informations sur le dataset
cat(sprintf("Forme des données: %d x %d\n", nrow(X), ncol(X)))
cat(sprintf("Nombre d'échantillons: %d\n", length(y)))
cat("Distribution des classes:\n")
for (i in 0:9) {
  cat(sprintf("  Chiffre %d: %d échantillons\n", i, sum(y == i)))
}

# Visualisons quelques exemples (simulés puisque nous n'avons pas d'images)
# Dans une application réelle, vous pourriez visualiser des images réelles
# Plot simple de quelques exemples pour chaque chiffre
par(mfrow=c(2,5), mar=c(2,2,2,2))
for (i in 0:9) {
  # Sélection aléatoire d'un exemple pour chaque chiffre
  idx <- sample(which(y == i), 1)
  # Création d'une visualisation simple des caractéristiques
  barplot(X[idx,], main=paste("Chiffre:", y[idx]),
          col="gray", border=NA, ylim=c(0, max(X)))
}

# ------------------------------------------------------------------------------
# 2. Définition des événements binaires pour l'estimation de proportions
# ------------------------------------------------------------------------------
# Pour estimer des proportions, nous devons définir des événements binaires (succès/échec)
# Nous allons considérer plusieurs cas:
# 1. La proportion de chiffres pairs
# 2. La proportion de chiffres supérieurs ou égaux à 5
# 3. La proportion d'un chiffre spécifique (ex: le chiffre 7)

# Définition des événements binaires
est_pair <- as.integer(y %% 2 == 0)  # 1 si le chiffre est pair, 0 sinon
est_grand <- as.integer(y >= 5)      # 1 si le chiffre est >= 5, 0 sinon
est_sept <- as.integer(y == 7)       # 1 si le chiffre est 7, 0 sinon

# Calcul des proportions réelles dans l'ensemble complet
p_pair_reel <- mean(est_pair)
p_grand_reel <- mean(est_grand)
p_sept_reel <- mean(est_sept)

cat(sprintf("Proportion réelle de chiffres pairs: %.4f\n", p_pair_reel))
cat(sprintf("Proportion réelle de chiffres >= 5: %.4f\n", p_grand_reel))
cat(sprintf("Proportion réelle de chiffres = 7: %.4f\n", p_sept_reel))

# ------------------------------------------------------------------------------
# 3. Estimation par intervalle de confiance pour une proportion
# ------------------------------------------------------------------------------
# Nous allons maintenant diviser nos données en ensemble d'apprentissage et de test,
# puis estimer les proportions sur l'ensemble d'apprentissage et construire des
# intervalles de confiance.

# Division en ensembles d'apprentissage et de test
set.seed(0)
indices <- createDataPartition(y, p = 0.7, list = FALSE)
y_train <- y[indices]
y_test <- y[-indices]

# Calcul des événements binaires sur l'ensemble d'apprentissage
est_pair_train <- as.integer(y_train %% 2 == 0)
est_grand_train <- as.integer(y_train >= 5)
est_sept_train <- as.integer(y_train == 7)

# Taille de l'échantillon
taille_echantillon <- length(y_train)
cat(sprintf("Taille de l'échantillon d'apprentissage: %d\n", taille_echantillon))

# ------------------------------------------------------------------------------
# 3.1 Intervalle de confiance avec l'approximation normale
# ------------------------------------------------------------------------------
# Pour n assez grand (généralement n*p ≥ 5 et n*(1-p) ≥ 5), on peut utiliser
# l'approximation normale. L'erreur standard de p-chapeau est donnée par
# sqrt(p(1-p)/n).

# Fonction pour calculer l'intervalle de confiance avec l'approximation normale
intervalle_confiance_normal <- function(echantillon, alpha = 0.05) {
  p_hat <- mean(echantillon)
  n <- length(echantillon)
  z_alpha_2 <- qnorm(1 - alpha/2)
  erreur_standard <- sqrt(p_hat * (1 - p_hat) / n)
  
  ic_inf <- max(0, p_hat - z_alpha_2 * erreur_standard)
  ic_sup <- min(1, p_hat + z_alpha_2 * erreur_standard)
  
  return(list(p_hat = p_hat, ic_inf = ic_inf, ic_sup = ic_sup))
}

# Calcul des estimations et intervalles de confiance
result_pair <- intervalle_confiance_normal(est_pair_train)
result_grand <- intervalle_confiance_normal(est_grand_train)
result_sept <- intervalle_confiance_normal(est_sept_train)

# Affichage des résultats
cat("Méthode de l'approximation normale:\n")
cat(sprintf("Proportion de chiffres pairs:\n"))
cat(sprintf("  - Proportion réelle: %.4f\n", p_pair_reel))
cat(sprintf("  - Proportion estimée: %.4f\n", result_pair$p_hat))
cat(sprintf("  - Intervalle de confiance à 95%%: [%.4f, %.4f]\n", 
          result_pair$ic_inf, result_pair$ic_sup))
cat(sprintf("  - La proportion réelle est dans l'intervalle: %s\n", 
          result_pair$ic_inf <= p_pair_reel && p_pair_reel <= result_pair$ic_sup))

cat(sprintf("\nProportion de chiffres >= 5:\n"))
cat(sprintf("  - Proportion réelle: %.4f\n", p_grand_reel))
cat(sprintf("  - Proportion estimée: %.4f\n", result_grand$p_hat))
cat(sprintf("  - Intervalle de confiance à 95%%: [%.4f, %.4f]\n", 
          result_grand$ic_inf, result_grand$ic_sup))
cat(sprintf("  - La proportion réelle est dans l'intervalle: %s\n", 
          result_grand$ic_inf <= p_grand_reel && p_grand_reel <= result_grand$ic_sup))

cat(sprintf("\nProportion de chiffres = 7:\n"))
cat(sprintf("  - Proportion réelle: %.4f\n", p_sept_reel))
cat(sprintf("  - Proportion estimée: %.4f\n", result_sept$p_hat))
cat(sprintf("  - Intervalle de confiance à 95%%: [%.4f, %.4f]\n", 
          result_sept$ic_inf, result_sept$ic_sup))
cat(sprintf("  - La proportion réelle est dans l'intervalle: %s\n", 
          result_sept$ic_inf <= p_sept_reel && p_sept_reel <= result_sept$ic_sup))

# Visualisons les intervalles de confiance pour les trois proportions
# Fonction pour visualiser l'intervalle de confiance
visualiser_ic <- function(p_hat, ic_inf, ic_sup, p_reel, titre) {
  erreur_standard <- (ic_sup - ic_inf) / (2 * qnorm(0.975))
  
  # Création d'une séquence pour tracer la distribution
  x <- seq(max(0, p_hat - 4*erreur_standard), min(1, p_hat + 4*erreur_standard), length.out = 1000)
  y <- dnorm(x, p_hat, erreur_standard)
  
  # Préparation des données pour ggplot2
  df <- data.frame(x = x, y = y)
  df_ic <- data.frame(x = seq(max(0, ic_inf), min(1, ic_sup), length.out = 100),
                     y = dnorm(seq(max(0, ic_inf), min(1, ic_sup), length.out = 100), p_hat, erreur_standard))
  
  # Création du graphique
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(size = 1, color = "blue") +
    geom_area(data = df_ic, aes(x = x, y = y), fill = "skyblue", alpha = 0.4) +
    geom_vline(xintercept = p_reel, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = p_hat, color = "green", size = 1) +
    geom_vline(xintercept = max(0, ic_inf), color = "blue", linetype = "dotted", size = 1) +
    geom_vline(xintercept = min(1, ic_sup), color = "blue", linetype = "dotted", size = 1) +
    labs(title = titre,
         x = "Valeur de la proportion",
         y = "Densité") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Ajout de la légende
  p <- p + annotate("text", x = p_reel, y = max(y)/2, label = "Proportion réelle", 
                  color = "red", hjust = -0.1)
  p <- p + annotate("text", x = p_hat, y = max(y)/3, label = "Proportion estimée", 
                  color = "darkgreen", hjust = -0.1)
  
  print(p)
}

# Visualisation des intervalles de confiance
visualiser_ic(result_pair$p_hat, result_pair$ic_inf, result_pair$ic_sup, p_pair_reel, 
            'Intervalle de confiance pour la proportion de chiffres pairs')
visualiser_ic(result_grand$p_hat, result_grand$ic_inf, result_grand$ic_sup, p_grand_reel, 
            'Intervalle de confiance pour la proportion de chiffres >= 5')
visualiser_ic(result_sept$p_hat, result_sept$ic_inf, result_sept$ic_sup, p_sept_reel, 
            'Intervalle de confiance pour la proportion de chiffres = 7')

# ------------------------------------------------------------------------------
# 3.2 Intervalle de Wilson
# ------------------------------------------------------------------------------
# L'intervalle de Wilson est une méthode plus précise pour estimer l'intervalle
# de confiance d'une proportion, notamment pour les petits échantillons ou les
# proportions proches de 0 ou 1.

# Fonction pour calculer l'intervalle de Wilson
intervalle_wilson <- function(echantillon, alpha = 0.05) {
  p_hat <- mean(echantillon)
  n <- length(echantillon)
  X <- sum(echantillon)  # Nombre de succès
  z_alpha_2 <- qnorm(1 - alpha/2)
  
  # Terme central
  centre <- (X + z_alpha_2^2/2) / (n + z_alpha_2^2)
  
  # Terme de marge
  marge <- z_alpha_2 * sqrt((X * (n - X) / n + z_alpha_2^2/4) / (n + z_alpha_2^2))
  
  # Limites de l'intervalle
  ic_inf <- max(0, centre - marge)
  ic_sup <- min(1, centre + marge)
  
  return(list(p_hat = p_hat, ic_inf = ic_inf, ic_sup = ic_sup))
}

# Calcul des estimations et intervalles de Wilson
result_pair_w <- intervalle_wilson(est_pair_train)
result_grand_w <- intervalle_wilson(est_grand_train)
result_sept_w <- intervalle_wilson(est_sept_train)

# Affichage des résultats
cat("Méthode de Wilson:\n")
cat(sprintf("Proportion de chiffres pairs:\n"))
cat(sprintf("  - Proportion réelle: %.4f\n", p_pair_reel))
cat(sprintf("  - Proportion estimée: %.4f\n", result_pair_w$p_hat))
cat(sprintf("  - Intervalle de Wilson à 95%%: [%.4f, %.4f]\n", 
          result_pair_w$ic_inf, result_pair_w$ic_sup))
cat(sprintf("  - La proportion réelle est dans l'intervalle: %s\n", 
          result_pair_w$ic_inf <= p_pair_reel && p_pair_reel <= result_pair_w$ic_sup))

cat(sprintf("\nProportion de chiffres >= 5:\n"))
cat(sprintf("  - Proportion réelle: %.4f\n", p_grand_reel))
cat(sprintf("  - Proportion estimée: %.4f\n", result_grand_w$p_hat))
cat(sprintf("  - Intervalle de Wilson à 95%%: [%.4f, %.4f]\n", 
          result_grand_w$ic_inf, result_grand_w$ic_sup))
cat(sprintf("  - La proportion réelle est dans l'intervalle: %s\n", 
          result_grand_w$ic_inf <= p_grand_reel && p_grand_reel <= result_grand_w$ic_sup))

cat(sprintf("\nProportion de chiffres = 7:\n"))
cat(sprintf("  - Proportion réelle: %.4f\n", p_sept_reel))
cat(sprintf("  - Proportion estimée: %.4f\n", result_sept_w$p_hat))
cat(sprintf("  - Intervalle de Wilson à 95%%: [%.4f, %.4f]\n", 
          result_sept_w$ic_inf, result_sept_w$ic_sup))
cat(sprintf("  - La proportion réelle est dans l'intervalle: %s\n", 
          result_sept_w$ic_inf <= p_sept_reel && p_sept_reel <= result_sept_w$ic_sup))

# ------------------------------------------------------------------------------
# 4. Comparaison des méthodes pour différentes tailles d'échantillon
# ------------------------------------------------------------------------------
# Nous allons maintenant comparer les deux méthodes (normale et Wilson) pour
# différentes tailles d'échantillon et évaluer leur performance en termes de
# taux de couverture et de largeur d'intervalle.

# Différentes tailles d'échantillon à tester
tailles <- c(20, 50, 100, 200, 500, 1000)
alpha <- 0.05  # Niveau de confiance de 95%
n_simulations <- 500  # Nombre de simulations pour chaque taille

# Nous allons nous concentrer sur la proportion de chiffres = 7 (cas d'une proportion faible)
p_reel <- p_sept_reel
est_event <- est_sept

# Préparation des structures pour stocker les résultats
couverture_norm <- numeric(length(tailles))
couverture_wilson <- numeric(length(tailles))
largeurs_norm <- numeric(length(tailles))
largeurs_wilson <- numeric(length(tailles))

for (i in 1:length(tailles)) {
  n <- tailles[i]
  compteur_norm <- 0
  compteur_wilson <- 0
  largeur_norm_total <- 0
  largeur_wilson_total <- 0
  
  for (j in 1:n_simulations) {
    # Échantillonnage avec remplacement
    indices <- sample(length(est_event), size = n, replace = TRUE)
    echantillon <- est_event[indices]
    
    # Méthode normale
    result_n <- intervalle_confiance_normal(echantillon, alpha)
    
    # Méthode de Wilson
    result_w <- intervalle_wilson(echantillon, alpha)
    
    # Vérification de la couverture
    if (result_n$ic_inf <= p_reel && p_reel <= result_n$ic_sup) {
      compteur_norm <- compteur_norm + 1
    }
    if (result_w$ic_inf <= p_reel && p_reel <= result_w$ic_sup) {
      compteur_wilson <- compteur_wilson + 1
    }
    
    # Calcul des largeurs
    largeur_norm_total <- largeur_norm_total + (result_n$ic_sup - result_n$ic_inf)
    largeur_wilson_total <- largeur_wilson_total + (result_w$ic_sup - result_w$ic_inf)
  }
  
  # Moyenne des résultats
  couverture_norm[i] <- compteur_norm / n_simulations
  couverture_wilson[i] <- compteur_wilson / n_simulations
  largeurs_norm[i] <- largeur_norm_total / n_simulations
  largeurs_wilson[i] <- largeur_wilson_total / n_simulations
  
  cat(sprintf("Taille %d terminée.\n", n))
}

# Visualisation des résultats
# Création d'un graphique à deux panneaux
par(mfrow=c(1,2), mar=c(5,4,4,2))

# Taux de couverture
plot(tailles, couverture_norm, type="b", pch=16, col="blue", log="x",
     main="Taux de couverture des intervalles (Chiffre 7)",
     xlab="Taille de l'échantillon", ylab="Taux de couverture")
lines(tailles, couverture_wilson, type="b", pch=16, col="red")
abline(h=1-alpha, lty=2, col="black")
legend("bottomright", legend=c("Méthode normale", "Méthode de Wilson", "Niveau nominal (95%)"),
       col=c("blue", "red", "black"), lty=c(1, 1, 2), pch=c(16, 16, NA))

# Largeur des intervalles
plot(tailles, largeurs_norm, type="b", pch=16, col="blue", log="x",
     main="Largeur moyenne des intervalles (Chiffre 7)",
     xlab="Taille de l'échantillon", ylab="Largeur moyenne")
lines(tailles, largeurs_wilson, type="b", pch=16, col="red")
legend("topright", legend=c("Méthode normale", "Méthode de Wilson"),
       col=c("blue", "red"), lty=1, pch=16)

# ------------------------------------------------------------------------------
# Rendu de TP
# ------------------------------------------------------------------------------

# Faire la même analyse sur le dataset digits en considérant:
# - La proportion du chiffre 1
# - La proportion des chiffres entre 3 et 7 inclus

