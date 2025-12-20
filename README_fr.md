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
- Space : afficher les paramètres et redessiner (affiche si Auto-scale est ON et son facteur)
- A / Z : diminuer/augmenter la distance (pas de 10 %)
- +/- : appliquer l'auto-fit si aucun auto-scale, puis augmenter/diminuer la distance (pas de 10 %)
- Flèches : ajuster les angles
- W / X : rotation écran
- C : changer palette
- F : basculer l'affichage en contours de polygones (par défaut : rempli — polygones affichés remplis + contour)
- R : annuler l'auto-scale (si appliqué)
- K : éditer angles/distance sans recharger le modèle (ENTER peut déclencher l'auto-fit)

## Remarques sur l'implémentation
- Les calculs critiques sont optimisés pour limiter les conversions flottantes et éviter les débordements (usage intensif de `Fixed32` et `Fixed64`).
- Un correctif d'orientation pour les exports OBJ Z-up a été ajouté (swap Y/Z à l'import) et peut être réverti manuellement.
- L'auto‑fit utilise désormais une **sphère englobante** précomputée (centroïde + rayon) pour estimer la distance en **O(1)** via `computeDistanceFromBoundingSphere()` ; `computeDistanceToFit()` (scan par sommet) est conservée mais **dépréciée** en tant que solution de secours.

## Fonctions importantes & références fichier:ligne

<!-- FUNC_LIST_START -->
- 🔧 `void painter_newell_sancha(Model3D* model, int face_count)` — `GS3Dp.cc:654`
- 🔧 `Model3D* createModel3D(void)` — `GS3Dp.cc:927`
- 🔧 `void destroyModel3D(Model3D* model)` — `GS3Dp.cc:1245`
- 🔧 `int loadModel3D(Model3D* model, const char* filename)` — `GS3Dp.cc:1306`
- 🔧 `void computeModelBoundingSphere(Model3D* model)` — `GS3Dp.cc:563`
- 🔧 `Fixed32 computeDistanceFromBoundingSphere(Model3D* model, float margin)` — `GS3Dp.cc:564`
- 🔧 `void getObserverParams(ObserverParams* params, Model3D* model)` — `GS3Dp.cc:574`
- 🔧 `void processModelFast(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:1498`
- 🔧 `void processModelWireframe(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:1603`
- 🔧 `int readVertices(const char* filename, VertexArrays3D* vtx, int max_vertices)` — `GS3Dp.cc:499`
- 🔧 `int readFaces_model(const char* filename, Model3D* model)` — `GS3Dp.cc:1760`
- 🔧 `void projectTo2D(VertexArrays3D* vtx, int angle_w_deg)` — `GS3Dp.cc:560`
- 🔧 `void calculateFaceDepths(Model3D* model, Face3D* faces, int face_count)` — `GS3Dp.cc:628`
- 🔧 `Fixed32 computeDistanceToFit(VertexArrays3D* vtx, float margin)` — `GS3Dp.cc:573`
- 🔧 `void autoScaleModel(Model3D* model, float target_max_dim, float min_scale, float max_scale, int center_flag)` — `GS3Dp.cc:577`
- 🔧 `void revertAutoScaleModel(Model3D* model)` — `GS3Dp.cc:578`
- 🔧 `void backupModelCoords(Model3D* model)` — `GS3Dp.cc:588`
- 🔧 `void freeBackupModelCoords(Model3D* model)` — `GS3Dp.cc:589`
- 🔧 `void fitModelToView(Model3D* model, ObserverParams* params, float target_max_dim, float margin, float percentile, int center_flag)` — `GS3Dp.cc:581`
- 🔧 `void drawPolygons(Model3D* model, int* vertex_count, int face_count, int vertex_count_total)` — `GS3Dp.cc:627`
- 🔧 `int main()` — `GS3Dp.cc:2689`
<!-- FUNC_LIST_END -->

## Crédits
- Auteur principal: Bruno
- Dédicace: *A tribute to Robert DONY* — Author of "Calcul des parties cachées" (Masson, 1986)

---

Si tu souhaites que j'ajoute un guide de compilation plus détaillé, des tests ou une licence, dis‑le et je l'ajoute.