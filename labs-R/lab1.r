# ===============================================================================
# Distribution conjointe de deux variables aléatoires
# ===============================================================================

# Dans cette partie, nous allons construire la distribution conjointe de deux variables aléatoires.
# La première correspond au prix du riz et la seconde correspond au prix de la viande.
# Nous prévoyons d'aller dans une épicerie aléatoire de la ville pour faire nos achats quotidiens.

# Nous supposons qu'il y a 9 grands grossistes de riz dans la ville et 6 grands grossistes de viande.
# Chaque épicerie de la ville a un contrat avec l'un des grossistes de riz et l'un des grossistes de viande.
# Nous supposons que l'épicerie vendrait au même prix que celui du grossiste.
# Enfin, nous supposons qu'aucune épicerie n'a les mêmes grossistes de viande et de riz en même temps.

# Voici les prix fixés par les grossistes de riz (axe des x) : 101, 102, 103, 104, 105, 106, 107, 108, 109
# Voici les prix fixés par les grossistes de viande (axe des y) : 101, 102, 103, 104, 105, 106

# Voici comment les épiceries sont distribuées dans la ville :
# - 1 épicerie a un contrat avec le grossiste de riz 101 et le grossiste de viande 105 -> 101, 105
# - 1 épicerie a un contrat avec le grossiste de riz 102 et le grossiste de viande 104 -> 102, 104
# - 1 épicerie a un contrat avec le grossiste de riz 102 et le grossiste de viande 105 -> 102, 105
# - 1 épicerie a un contrat avec le grossiste de riz 103 et le grossiste de viande 103 -> 103, 103
# - 1 épicerie a un contrat avec le grossiste de riz 103 et le grossiste de viande 104 -> 103, 104
# - 1 épicerie a un contrat avec le grossiste de riz 103 et le grossiste de viande 105 -> 103, 105
# - 1 épicerie a un contrat avec le grossiste de riz 104 et le grossiste de viande 102 -> 104, 102
# - 1 épicerie a un contrat avec le grossiste de riz 104 et le grossiste de viande 103 -> 104, 103
# - 1 épicerie a un contrat avec le grossiste de riz 104 et le grossiste de viande 104 -> 104, 104
# - 1 épicerie a un contrat avec le grossiste de riz 104 et le grossiste de viande 105 -> 104, 105
# - 1 épicerie a un contrat avec le grossiste de riz 104 et le grossiste de viande 106 -> 104, 106
# - 1 épicerie a un contrat avec le grossiste de riz 105 et le grossiste de viande 101 -> 105, 101  
# - 1 épicerie a un contrat avec le grossiste de riz 105 et le grossiste de viande 102 -> 105, 102
# - 1 épicerie a un contrat avec le grossiste de riz 105 et le grossiste de viande 103 -> 105, 103
# - 1 épicerie a un contrat avec le grossiste de riz 105 et le grossiste de viande 104 -> 105, 104
# - 1 épicerie a un contrat avec le grossiste de riz 105 et le grossiste de viande 105 -> 105, 105
# - 1 épicerie a un contrat avec le grossiste de riz 105 et le grossiste de viande 106 -> 105, 106
# - 1 épicerie a un contrat avec le grossiste de riz 106 et le grossiste de viande 102 -> 106, 102
# - 1 épicerie a un contrat avec le grossiste de riz 106 et le grossiste de viande 103 -> 106, 103
# - 1 épicerie a un contrat avec le grossiste de riz 106 et le grossiste de viande 104 -> 106, 104
# - 1 épicerie a un contrat avec le grossiste de riz 106 et le grossiste de viande 105 -> 106, 105
# - 1 épicerie a un contrat avec le grossiste de riz 106 et le grossiste de viande 106 -> 106, 106  
# - 1 épicerie a un contrat avec le grossiste de riz 107 et le grossiste de viande 103 -> 107, 103
# - 1 épicerie a un contrat avec le grossiste de riz 107 et le grossiste de viande 104 -> 107, 104
# - 1 épicerie a un contrat avec le grossiste de riz 107 et le grossiste de viande 105 -> 107, 105
# - 1 épicerie a un contrat avec le grossiste de riz 108 et le grossiste de viande 104 -> 108, 104
# - 1 épicerie a un contrat avec le grossiste de riz 108 et le grossiste de viande 105 -> 108, 105
# - 1 épicerie a un contrat avec le grossiste de riz 109 et le grossiste de viande 105 -> 109, 105

# Importation des bibliothèques nécessaires
library(ggplot2)

# Création des combinaisons d'épiceries
combinations <- matrix(c(
  101, 105,
  102, 104,
  102, 105,
  103, 103,
  103, 104,
  103, 105,
  104, 102,
  104, 103,
  104, 104,
  104, 105,
  104, 106,
  105, 101,
  105, 102,
  105, 103,
  105, 104,
  105, 105,
  105, 106,
  106, 102,
  106, 103,
  106, 104,
  106, 105,
  106, 106,
  107, 103,
  107, 104,
  107, 105,
  108, 104,
  108, 105,
  109, 105
), ncol = 2, byrow = TRUE)

combinations_i <- combinations[, 1]
combinations_j <- combinations[, 2]

# Tracer les données dispersées
plot_data <- data.frame(x = combinations_i, y = combinations_j)
p <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_point() +
  labs(
    x = "prix du riz en MRU",
    y = "prix de la viande en MRU",
    title = "Nuage de points de x et y"
  ) +
  xlim(100, 110) +
  ylim(100, 107)
print(p)

# ===============================================================================
# Nous déduisons que nous établissons une distribution uniforme à partir des points ci-dessus
# - calculer la probabilité de chaque point ?
# - quel serait le montant moyenne à payer si nous allons dans une épicerie aléatoire de la ville ?
# ===============================================================================

# Nous avons la possibilité d'aller dans une épicerie choisie au hasard. Et nous avons que chaque point correspond à une et une seule épicerie.
# Ainsi, nous avons que la probabilité de chaque point est 1/28.
# Nous pouvons calculer la valeur moyenne en additionnant le produit du prix du riz et de la viande par la probabilité de chaque point.

# Maintenant, nous allons calculer le montant moyen à payer si nous achetons 1kg de riz et 1kg de viande
# D'abord, calculons la somme x + y pour chaque point (épicerie)
totals <- rowSums(combinations)

for (i in 1:length(totals)) {
  cat(sprintf("point : [%d, %d], total à payer : %d\n", combinations[i, 1], combinations[i, 2], totals[i]))
}

# En moyenne, nous devrons payer : S_i (x_i + y_i) / n
average_total <- mean(totals)
cat(sprintf("Montant total moyen à payer : %.2f\n", average_total))

# ===============================================================================
# La distribution marginale de x
# ===============================================================================

# Maintenant, nous ne sommes intéressés que par l'achat de riz (noté x).
# - Dessiner un histogramme de la distribution marginale de x.
# - Calculer et dessiner la distribution cumulée de x.
# - Quel est le prix moyen du riz ?

# Collecter les données
x <- combinations_i

# Tracer l'histogramme
hist(x, breaks = seq(100, 110, by = 1), 
     xlab = "prix du riz en MRU", 
     ylab = "nombre d'épiceries avec ce prix", 
     main = "Histogramme de x",
     xlim = c(99, 110))

# Tracer la probabilité
hist(x, breaks = seq(100, 110, by = 1), 
     xlab = "prix du riz en MRU", 
     ylab = "probabilité", 
     main = "Probabilité marginale de x",
     xlim = c(99, 110),
     ylim = c(0, 0.25),
     freq = FALSE)

# Calculer les probabilités pour chaque valeur unique de x
x_tab <- table(x)
x_probs <- x_tab / sum(x_tab)

for (i in 1:length(x_probs)) {
  cat(sprintf("Probabilité de x = %d est %d/%d = %.4f\n", 
              as.numeric(names(x_probs)[i]), 
              x_tab[i], 
              length(x), 
              x_probs[i]))
}

# Tracer la distribution cumulée
plot(ecdf(x), 
     xlab = "prix du riz en MRU", 
     ylab = "distribution cumulée", 
     main = "Distribution cumulée marginale de x",
     xlim = c(99, 110),
     ylim = c(0, 1.1))

# Le prix moyen à payer pour le riz si nous achetons 1kg de riz dans une épicerie aléatoire est S_i x_i / n
average_price <- mean(x)
cat(sprintf("Prix moyen à payer pour le riz : %.2f\n", average_price))

# ===============================================================================
# Distribution marginale de y (exercice) rendu de TP
# ===============================================================================

# Notre compagnon n'est intéressé que par l'achat de viande (notée y).
# - Dessiner un histogramme de la distribution marginale de y.
# - Calculer et dessiner la distribution cumulée de y.
# - Quel est le prix moyen que notre compagnon paierait pour la viande s'il allait dans une épicerie aléatoire de la ville ?

# ===============================================================================
# Distribution conditionnelle de x sachant y
# ===============================================================================

# Quand nous sommes finalement arrivés au centre-ville pour acheter du riz, 
# quelqu'un nous a dit que le prix de la viande dans la région est de 105.
# - Dessiner un histogramme de la distribution conditionnelle de x (le prix du riz) 
#   sachant que si nous devions acheter de la viande, elle coûterait 105 (y).
# - Calculer et dessiner la distribution cumulée de x sachant que y=105.
# - Quel est le prix moyen du riz sachant que le prix de la viande est de 105 ?

# Tracer l'histogramme pour la distribution conditionnelle de x sachant y=105
x_given_y_105 <- combinations_i[combinations_j == 105]

hist(x_given_y_105, breaks = seq(100, 110, by = 1), 
     xlab = "prix du riz en MRU", 
     ylab = "nombre d'épiceries avec ce prix", 
     main = "Histogramme de x sachant y=105",
     xlim = c(99, 110))

# Tracer la probabilité
hist(x_given_y_105, breaks = seq(100, 110, by = 1), 
     xlab = "prix du riz en MRU", 
     ylab = "probabilité", 
     main = "Probabilité conditionnelle de x sachant y=105",
     xlim = c(99, 110),
     freq = FALSE)

# Calculer les probabilités pour chaque valeur unique de x sachant y=105
x_given_y_105_tab <- table(x_given_y_105)
x_given_y_105_probs <- x_given_y_105_tab / sum(x_given_y_105_tab)

for (i in 1:length(x_given_y_105_probs)) {
  cat(sprintf("Probabilité de x = %d sachant y=105 est %d/%d = %.4f\n", 
              as.numeric(names(x_given_y_105_probs)[i]), 
              x_given_y_105_tab[i], 
              length(x_given_y_105), 
              x_given_y_105_probs[i]))
}

# Tracer la distribution cumulée de x sachant y=105
plot(ecdf(x_given_y_105), 
     xlab = "prix du riz en MRU", 
     ylab = "distribution cumulée", 
     main = "Distribution cumulée conditionnelle de x sachant y=105",
     xlim = c(99, 110),
     ylim = c(0, 1.1))

# Calculer le prix moyen à payer pour le riz sachant y=105
average_price_given_y_105 <- mean(x_given_y_105)
cat(sprintf("Sachant que si nous devions acheter 1kg de viande cela coûterait 105 MRU, le prix moyen à payer pour le riz est de %.2f MRU\n", average_price_given_y_105))

# ===============================================================================
# Distribution conditionnelle de y sachant x (exercice) rendu de TP
# ===============================================================================

# Notre compagnon est allé au centre-ville pour acheter de la viande 
# et quelqu'un lui a dit que le prix du riz dans la région est de 103.
# - Dessiner un histogramme de la distribution conditionnelle de y (le prix de la viande) 
#   sachant que s'il devait acheter du riz, cela lui coûterait 103 (x).
# - Calculer et dessiner la distribution cumulée de y sachant que x=103.
# - Quel est le prix moyen de la viande sachant que le prix du riz est de 103 ?

