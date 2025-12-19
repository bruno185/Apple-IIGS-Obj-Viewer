# 3D OBJ Viewer (Apple IIGS)

Un visualiseur 3D pour Apple IIGS écrit en ORCA/C utilisant une arithmétique Fixed32 (16.16).

## Description
Ce projet lit des fichiers OBJ simplifiés (sommets `v` et faces `f`), effectue des transformations 3D, réalise une projection perspective en 2D et dessine les polygones avec l'API QuickDraw de l'Apple IIGS.

## Fonctionnalités principales
- Lecture de fichiers OBJ (sommets et faces)
- Transformations 3D optimisées en Fixed32 (16.16)
- Tables trigonométriques pour vitesse
- Algorithme de peinture (painter's algorithm) avec tests Newell/Sancha
- Auto-scale non-destructif à l'import (optionnel) avec possibilité de revenir (touche `r`)
- Options interactives (angles, distance, palette de couleurs)

## Utilisation
1. Compiler avec votre chaîne ORCA/C ou outil compatible (exemples dans le dépôt). Exemple d'utilisation dans l'environnement de développement:

   iix compile GS3Dp.cc

2. Lancer l'exécutable généré. Le programme demandera le fichier OBJ à charger et proposera d'appliquer l'auto-scale.

3. Commandes clavier:
- Space : afficher les paramètres et redessiner
- A / Z : diminuer/augmenter la distance
- Flèches : ajuster les angles
- W / X : rotation écran
- C : changer palette
- R : annuler l'auto-scale (si appliqué)

## Remarques sur l'implémentation
- Les calculs critiques sont optimisés pour limiter les conversions flottantes et éviter les débordements (usage intensif de `Fixed32` et `Fixed64`).
- Un correctif d'orientation pour les exports OBJ Z-up a été ajouté (swap Y/Z à l'import) et peut être réverti manuellement.
- `computeDistanceToFit` a été optimisé pour utiliser des min/max en `Fixed32` puis une conversion unique.

## Crédits
- Auteur principal: Bruno
- Dédicace: *A tribute to Robert DONY* — Author of "Calcul des parties cachées" (Masson, 1986)

---

Si tu souhaites que j'ajoute un guide de compilation plus détaillé, des tests ou une licence, dis‑le et je l'ajoute.