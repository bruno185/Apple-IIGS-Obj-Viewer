# 3D OBJ Viewer (Apple IIGS)

A 3D viewer for the Apple IIGS written in ORCA/C using Fixed32 (16.16) arithmetic.

## Description
This project reads simplified OBJ files (vertices `v` and faces `f`), performs 3D transformations, projects to 2D with perspective, and draws polygons using the Apple IIGS QuickDraw API.

## Main features
- Read OBJ files (vertices and faces)
- Optimized 3D transforms using Fixed32 (16.16)
- Fast trig tables for speed
- Painter's algorithm with Newell/Sancha overlap tests
- Non-destructive auto-scale on import (optional) with revert (press `r`)
- Interactive options (angles, distance, color palette)

## Usage
1. Compile with your ORCA/C toolchain or a compatible tool (examples in the repo). Example:

   iix compile GS3Dp.cc

2. Run the generated executable. The program will ask for an OBJ filename and optionally apply auto-scale.

3. Keyboard controls:
- Space: show parameters and redraw (shows whether auto-scale is ON and its factor)
- A / Z: decrease/increase distance (10% steps)
- +/-: apply auto-fit if none present then increase/decrease distance (10% steps)
- Arrow keys: adjust angles
- W / X: rotate screen
- C: change palette
- F: toggle fast painter (default: ON — tests 1-3 only)
- P: toggle frame-only polygons (default: OFF)
- R: revert auto-scale (if applied)
- K: edit angles/distance interactively without reloading the model (ENTER may trigger auto-fit)

## Implementation notes
- Critical computations are optimized to reduce floating conversions and avoid overflow (heavy use of `Fixed32` and `Fixed64`).
- An orientation fix was added for OBJ Z-up exports (swap Y/Z at import) and can be reverted manually.
- Auto‑fit uses a precomputed **bounding sphere** (centroid + radius) for O(1) distance estimates via `computeDistanceFromBoundingSphere()`; `computeDistanceToFit()` (per‑vertex scan) is retained but **deprecated** as a fallback.

## Important functions & file:line references

<!-- FUNC_LIST_START -->
- 🔧 `void painter_newell_sancha_fast(Model3D* model, int face_count)` — `GS3Dp.cc:659`
- 🔧 `void painter_newell_sancha(Model3D* model, int face_count)` — `GS3Dp.cc:705`
- 🔧 `void dumpFaceEquationsCSV(Model3D* model, const char* csv_filename)` — `GS3Dp.cc:2174`
- 🔧 `Model3D* createModel3D(void)` — `GS3Dp.cc:1036`
- 🔧 `void destroyModel3D(Model3D* model)` — `GS3Dp.cc:1354`
- 🔧 `int loadModel3D(Model3D* model, const char* filename)` — `GS3Dp.cc:1415`
- 🔧 `void computeModelBoundingSphere(Model3D* model)` — `GS3Dp.cc:564`
- 🔧 `Fixed32 computeDistanceFromBoundingSphere(Model3D* model, float margin)` — `GS3Dp.cc:565`
- 🔧 `void getObserverParams(ObserverParams* params, Model3D* model)` — `GS3Dp.cc:575`
- 🔧 `void processModelFast(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:1606`
- 🔧 `void processModelWireframe(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:1687`
- 🔧 `int readVertices(const char* filename, VertexArrays3D* vtx, int max_vertices)` — `GS3Dp.cc:500`
- 🔧 `int readFaces_model(const char* filename, Model3D* model)` — `GS3Dp.cc:1821`
- 🔧 `void projectTo2D(VertexArrays3D* vtx, int angle_w_deg)` — `GS3Dp.cc:561`
- 🔧 `void calculateFaceDepths(Model3D* model, Face3D* faces, int face_count)` — `GS3Dp.cc:629`
- 🔧 `Fixed32 computeDistanceToFit(VertexArrays3D* vtx, float margin)` — `GS3Dp.cc:574`
- 🔧 `void autoScaleModel(Model3D* model, float target_max_dim, float min_scale, float max_scale, int center_flag)` — `GS3Dp.cc:578`
- 🔧 `void revertAutoScaleModel(Model3D* model)` — `GS3Dp.cc:579`
- 🔧 `void backupModelCoords(Model3D* model)` — `GS3Dp.cc:589`
- 🔧 `void freeBackupModelCoords(Model3D* model)` — `GS3Dp.cc:590`
- 🔧 `void fitModelToView(Model3D* model, ObserverParams* params, float target_max_dim, float margin, float percentile, int center_flag)` — `GS3Dp.cc:582`
- 🔧 `void drawPolygons(Model3D* model, int* vertex_count, int face_count, int vertex_count_total)` — `GS3Dp.cc:628`
- 🔧 `int main()` — `GS3Dp.cc:2860`
<!-- FUNC_LIST_END -->

## Credits
- Main author: Bruno
- Tribute: *A tribute to Robert DONY* — Author of "Calcul des parties cachées" (Masson, 1986)

---

If you'd like a more detailed build guide, tests, or a license added, say so and I'll add it.