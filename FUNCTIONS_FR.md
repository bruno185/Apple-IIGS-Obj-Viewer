# Liste des fonctions principales (résumé en français)

Ce fichier répertorie les fonctions majeures trouvées dans `GS3Dp.cc` et donne un court résumé de leur rôle. Il est destiné à faciliter la restructuration du code.

---

## Fonctions utilitaires

- `sin_deg_int(int deg)` / `cos_deg_int(int deg)`
  - Retourne sinus/cosinus pour un angle entier (degrés) via table de lookup (Fixed32).

- `fixed32_to_ascii(Fixed32 val, char *buf, int bufsize, int decimals)`
  - Convertit un `Fixed32` (16.16) en chaîne ASCII sans utiliser `%f` (précision contrôlée).

- `is_ascii_numeric(const char* s)`
  - Vérifie qu'une chaîne ne contient que des caractères numériques valides (chiffres, point, e/E, signes).

- `normalize_deg(int deg)`
  - Ramène un angle en degrés dans l'intervalle [0, 359].

- Macros/fonctions `FIXED_*` (p.ex. `FLOAT_TO_FIXED`, `FIXED_TO_FLOAT`, `FIXED_MUL_64`)
  - Opérations et conversions pour l'arithmétique en virgule fixe 16.16.

---

## Gestion du modèle 3D

- `Model3D* createModel3D(void)`
  - Alloue et initialise les structures (sommets, faces, buffers) pour un modèle 3D.

- `void destroyModel3D(Model3D* model)`
  - Libère toute la mémoire allouée par `createModel3D`.

- `int loadModel3D(Model3D* model, const char* filename)`
  - Charge un fichier OBJ (sommets + faces), met à jour les compteurs et calcule la bounding sphere.

- `void backupModelCoords(Model3D* model)` / `void freeBackupModelCoords(Model3D* model)`
  - Sauvegarde / restaure ou libère une copie des coordonnées du modèle (utilisé pour transformations temporaires).

- `void computeModelBoundingSphere(Model3D* model)`
  - Calcule le centre et rayon approximatif du modèle (pour auto-fit, heuristiques).

- `Fixed32 computeDistanceFromBoundingSphere(Model3D* model, float margin)`
  - Estime une distance d'observation basée sur la bounding sphere et une marge.

---

## Lecture / parsing

- `int readVertices(const char* filename, VertexArrays3D* vtx, int max)`
  - Lit les sommets (v) depuis un OBJ, les convertit en `Fixed32` et les stocke.

- `int readFaces_model(const char* filename, Model3D* model)`
  - Lit les faces (f) depuis un OBJ, remplit le buffer d'indices, gère la mémoire en chunks.

---

## Projection et traitement géométrique

- `void projectTo2D(VertexArrays3D* vtx, int angle_w_deg)`
  - Transforme les coordonnées observateur en coordonnées écran (projection perspective).

- `void calculateFaceDepths(Model3D* model, Face3D* faces, int face_count)`
  - Calcule `z_min`, `z_max`, `z_mean`, les coefficients de plan (`a,b,c,d`), bounding box 2D et `display_flag` par face.

- `void computeDistanceToFit(Model3D* model, ObserverParams* params)`
  - (Helper) calcule la distance nécessaire pour faire tenir le modèle dans la vue.

- `void autoScaleModel(Model3D* model, float target_max_dim, float min_scale, float max_scale, int center_flag)`
  - Applique un redimensionnement automatique et éventuellement recentre le modèle.

- `void revertAutoScaleModel(Model3D* model)`
  - Restaure les coordonnées d'origine après un `autoScale`.

- `void fitModelToView(Model3D* model, ObserverParams* params, float target_max_dim, float margin, float percentile, int center_flag)`
  - Procédure haut niveau pour adapter modèle et paramètres d'observateur afin de fit l'écran.

---

## Algorithmes d'ordonnancement / painter

- `void painter_newell_sancha_fast(Model3D* model, int face_count)`
  - Variante rapide (optimisée) du painter (corrections limitées).

- `void painter_newell_sancha(Model3D* model, int face_count)`
  - Algorithme central d'ordonnancement: tri initial par `z_mean`, puis corrections itératives (tests 1..7), journalisation (`PAINTLOG.TXT`), détection d'oscillation et stratégies de secours (Pascal fallback, bail-out, stable sort).

- **(Non présentes)** `int simulateDevant_obs(...)` / `int simulateDerriere_obs(...)`
  - Référence dans la documentation : fonctions de simulation des tests Pascal. **Elles ne sont pas présentes** dans l'état courant du dépôt ; il faut les réimplémenter si l'on veut comparer/comprendre l'ordre de Pascal en espace observateur.

- **(Non présente)** `int compute_pascal_order(FaceArrays3D* faces, VertexArrays3D* vtx, int n, int* out_idx)`
  - Devrait calculer l'ordre de rendu théorique à la Pascal (full pairwise insertion) pour diagnostic / fallback. **Absente** actuellement — proposition : ajouter une implémentation qui simule TestDevant/TestDerriere pour toutes les paires et retourne un ordre stable pour comparaison.

- `int cmp_faces_by_zmean(const void* a, const void* b)`
  - Comparateur utilisé par `qsort` pour trier faces selon `z_mean` (tie-breaker: indice).

---

## Sortie et diagnostic

- `void dumpFaceEquationsCSV(Model3D* model, const char* csv_filename)`
  - Exporte les coefficients de chaque face (`a,b,c,d` + indices et z des 3 sommets) dans `equ.csv` (ou fichier passé en paramètre).

---

## Rendu & UI

- `void drawPolygons(Model3D* model, int* vertex_count, int face_count, int vertex_count_total)`
  - Dessine les polygones triés sur l'écran (en utilisant QuickDraw/l'API IIGS).

- `void DoColor(...)` / `void DoText(...)` (helpers)
  - Routines d'affichage de couleur et de texte pour l'interface.

- `void getObserverParams(ObserverParams* params, Model3D* model)`
  - Interface textuelle pour saisir/éditer les paramètres d'observation (angles, distance, etc.).

- `void processModelFast(...)` / `void processModelWireframe(...)`
  - Pipelines de traitement (préparation des données et rendu rapide ou wireframe).

---

## Entrée / sortie principale

- `int main()`
  - Point d'entrée, boucle principale UI, gestion des interactions clavier, charge les modèles et appelle les pipelines.

---

Notes:
- Cette liste est volontairement concise; certaines fonctions utilitaires ou variantes (suffixées `_tmp` ou versions de debug) ne sont pas détaillées ici.
- Dites-moi si vous voulez :
  - une version plus détaillée (paramètres et types pour chaque fonction),
  - ou que je sépare les fonctions par modules (ex: `io`, `math`, `painter`, `render`).
