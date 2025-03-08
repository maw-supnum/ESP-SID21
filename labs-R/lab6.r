# Régression Linéaire
# 
# Dans un modèle linéaire, nous supposons que la valeur à prédire y est une combinaison linéaire
# des variables d'entrée X ∈ R^p plus un bruit ε. Mathématiquement, cela peut être exprimé comme :
# 
# y = β₀ + β₁X + ε
# 
# où :
# - y est la variable dépendante (la valeur à prédire),
# - X est le vecteur de valeurs représentant les variables d'entrée,
# - β₀ est l'ordonnée à l'origine,
# - β₁ est le vecteur des coefficients pour les variables d'entrée,
# - ε est le terme d'erreur (bruit).
# 
# Nous pouvons supposer que β₀ = 0 et E(ε) = 0, car nous pouvons simplement étendre X à R^(p+1),
# en ajoutant une colonne de uns qui pourrait correspondre à la somme de l'ordonnée à l'origine
# et de l'espérance du terme d'erreur.

## Estimateur des Moindres Carrés

# importations
library(ggplot2)

# définir la graine aléatoire
set.seed(0)

# nombre d'échantillons
n <- 100

# nombre de variables d'entrée
p <- 1

# données d'entrée supposées suivre une distribution normale
# la première colonne est composée de uns (ordonnée à l'origine + espérance de l'erreur)
# nous définissons les données d'entrée avec une moyenne de 0 et un écart-type de 4
sigma <- 4
X <- cbind(rep(1, n), sigma * rnorm(n, mean = 0, sd = 1))

# générer le terme d'erreur
# nous définissons le terme d'erreur avec une moyenne de 0 et un écart-type de 3
epsilon <- 3 * rnorm(n)

# définir beta
beta <- c(3.2, 2.1)

# générer la variable de réponse
y <- X %*% beta + epsilon

# tracer les données
plot(X[, 2], y, xlab = "X", ylab = "y", main = "Données de régression")

# Nous pouvons déduire les estimations des moindres carrés des coefficients β₀ et β₁ en utilisant
# les formules suivantes : β* = (X^TX)^(-1)X^Ty

# déduire la pente et l'ordonnée à l'origine
beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y

# afficher la pente et l'ordonnée à l'origine estimées
cat("Pente estimée: ", beta_hat[2], "\n")
cat("Ordonnée à l'origine estimée: ", beta_hat[1], "\n")

# afficher la pente et l'ordonnée à l'origine réelles
cat("Pente réelle: ", beta[2], "\n")
cat("Ordonnée à l'origine réelle: ", beta[1], "\n")

# tracer les données et la droite ajustée
plot(X[, 2], y, xlab = "X", ylab = "y", main = "Régression linéaire")
abline(a = beta_hat[1], b = beta_hat[2], col = "red")

## Test de Fisher
#
# Nous voulons tester l'hypothèse nulle selon laquelle le coefficient β₁ est égal à zéro.
# Nous pouvons utiliser le test de Fisher pour cela. La statistique de test est donnée par :
#
# F = (R² × (n-p-1))/((1-R²)p)
#
# avec R² = Σᵢ(β*.Xᵢ - ȳ)² / Σᵢ(β*.Xᵢ - yᵢ)²

# calculer le R²
y_mean <- mean(y)
R_2 <- sum((y_mean - X %*% beta_hat)^2) / sum((y - y_mean)^2)

# calculer la statistique F
F_stat <- R_2 * (n - p - 1) / ((1 - R_2) * p)

# afficher les résultats
cat("Le seuil de rejet est : ", qf(0.95, p, n - p - 1), "\n")
cat("La statistique F est : ", F_stat, "\n")
cat("La p-valeur est : ", 1 - pf(F_stat, p, n - p - 1), "\n")

# Nous construisons un cadre où l'hypothèse nulle est vraie. Nous générons un ensemble de données
# avec une relation linéaire entre la variable d'entrée et la variable de sortie, et nous testons
# si le test de Fisher peut rejeter l'hypothèse nulle.

# y ne dépend pas de la deuxième colonne de X
y_test_2 <- epsilon

# déduire la pente et l'ordonnée à l'origine
beta_hat_2 <- solve(t(X) %*% X) %*% t(X) %*% y_test_2

# calculer le R²
y_test_2_mean <- mean(y_test_2)
R_2_2 <- sum((y_test_2_mean - X %*% beta_hat_2)^2) / sum((y_test_2 - y_test_2_mean)^2)

# calculer la statistique F
F_2 <- R_2_2 * (n - p - 1) / ((1 - R_2_2) * p)

# afficher les résultats
cat("Le seuil de rejet est : ", qf(0.95, p, n - p - 1), "\n")
cat("La statistique F est : ", F_2, "\n")
cat("La p-valeur est : ", 1 - pf(F_2, p, n - p - 1), "\n")

## Test de Student
#
# Nous pouvons également tester l'hypothèse nulle selon laquelle le coefficient β₁ est égal à zéro
# en utilisant le test de Student. La statistique de test est donnée par :
#
# t = β₁/σ₁

## Rendu de TP
#
# Générez trois cas :
# - Un cas où l'hypothèse nulle est vraie pour toutes les composantes, p=2 et n=100.
# - Un cas où l'hypothèse nulle est vraie pour une composante, p=2 et n=100.
# - Un cas où l'hypothèse nulle est fausse pour toutes les composantes, p=2 et n=100.
#
# Fixez la graine à 0 à des fins de reproductibilité.
