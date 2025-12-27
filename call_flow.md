# Call Flow (start at `main`) ✅

> Visual style:
> - **Function** names are bolded and use `code` formatting
> - *Purpose* and short description follow on the same line
> - Subcalls are indented and prefixed with `→` for readability

---

## Main flow

### Startup / Load

- 🔧 **`main()`** — Entry point
  - → 🔧 **`createModel3D()`** — allocate arrays, init metadata (auto-scale / backup pointers)
  - → 🔧 **`loadModel3D(model, filename)`** — orchestrates model loading
    - → 🔧 **`readVertices()`** — read `v x y z` into Fixed arrays
    - → 🔧 **`readFaces_model()`** — parse faces into packed index buffers
    - → 🔧 **`computeModelBoundingSphere()`** — compute model bounding sphere (centroid + radius) **O(n)** (done once)

### Parameter parsing / Auto-fit

- 🔧 **`getObserverParams(&params, model)`** — interactive parameter parsing
  - reads angles H/V/W and screen rotation
  - If user presses ENTER for distance: **auto-fit path**
    - 🔧 **`fitModelToView(model, params, target_max_dim, margin, percentile, center_flag)`** — fits model using sampled radii
      - sample vertices (up to `max_samples`) to compute centroid + squared radii
      - quickselect percentile radius → compute scale and center
      - `backupModelCoords(model)` and *apply scale+center* to vertices (non-destructive)
      - update `model->auto_scale`, `auto_center_*`, `auto_scaled`
      - update model bounding sphere (`bs_cx/bs_cy/bs_cz/bs_r`) accordingly — now kept in sync **O(1)**
      - Distance estimation is disabled here; `params->distance` must be set explicitly by the user
      - **FALLBACK:** none — per-vertex auto-fit removed
  - else: set `params->distance` from user input

---

## Main render loop

- Enter main render loop (`bigloop`):
  - 🔧 **`processModelFast(model, &params, filename)`** — runs every frame; ultra-fast transformation + projection
    - precompute trig products (Fixed32)
    - For each vertex (tight Fixed32 loop): transform → compute `xo/yo/zo` → project to `x2d/y2d`
    - **Note:** runtime auto-fit is *not* applied inside `processModelFast` — any autoscale is applied ahead of time by 🔧 **`fitModelToView()`**, which modifies model coordinates and updates the bounding sphere; `processModelFast` operates on those (possibly scaled) model vertices.
    - 🔧 **`calculateFaceDepths()`** — compute per-face `z_min/z_max/z_mean`, display flags, planar coefficients (Newell)
    - 🔧 **`painter_newell_sancha()`** — sort faces by depth and correct ambiguous order (qsort + corrections)

    - 🔧 **`processModelWireframe(model, &params, filename)`** — lightweight wireframe processing: transforms and projects vertices and sets face visibility only (no per-face depth/sorting); used for wireframe/frame-only display to improve speed.

- 🔧 **`drawPolygons(model, faces, face_count, vert_count)`** — render loop
  - uses sorted faces; for each face builds QuickDraw polygon from `vtx->x2d/y2d` (these already reflect any model-space auto-scaling applied by `fitModelToView`)
  - Fill + Frame polygon (QuickDraw)

### UI / Input handling

- `startgraph()` / render / `endgraph()` / `DoText()` / optional `DoColor()`
- Keys: Space (info), N (new model), Arrows / A Z (angles/distance), `r` (revert autoscale), `K` (edit angles/distance without reloading model), `+`/`-` (adjust autoscale)
- `K` invokes `getObserverParams(&params, model)` interactively and applies new angles/distance without requiring a reload.
- `+`/`-` behavior: if the model is not yet auto-scaled, these keys first perform an **auto-fit** (`fitModelToView()`); they then increase or decrease the current `params->distance` and update `model->auto_scale`. Automatic recomputation via bounding-sphere is disabled; distance adjustments are manual. The action prints a short message (e.g., "Distance increased" / "Distance decreased").

---

## Supporting / fallback functions (short purpose)

- 🔧 **`projectTo2D()`** — standalone projection helper (used in places)
- 🔧 **`computeModelBoundingSphere(model)`** — O(n) at load, stores `bs_*` in `Model3D`
- 🔧 **`autoScaleModel()` / `revertAutoScaleModel()`** — non-destructive scaling helpers (backup + apply + revert)
- 🔧 **`backupModelCoords()` / `freeBackupModelCoords()`** — support for non-destructive transforms
- 🔧 **`processModelWireframe(model, &params, filename)`** — lightweight wireframe processing (transform & project only) used for fast wireframe/frame-only rendering
- 🔧 **`destroyModel3D(Model3D* model)`** — frees all memory allocated by `createModel3D()`; must be called to avoid leaks
- 🔧 **`readVertices()` / `readFaces_model()`** — file parsing helpers

**Notes:** `adjustDistanceFast()` (earlier fast-adjust prototype) was removed — automatic distance estimation has been disabled; distance adjustments are manual.

---

## Notes

- This Markdown preserves the visual grouping: **function name** + **purpose/action**, and uses `→` to show call flow.
- If you want more detail (e.g., exact file:line references or prototypes), I can append them in a table or add a per-function section with signatures.

---

## Function signatures & file:line references

- 🔧 `void painter_newell_sancha(Model3D* model, int face_count)` — `GS3Dp.cc:654`
- 🔧 `Model3D* createModel3D(void)` — `GS3Dp.cc:927`
- 🔧 `void destroyModel3D(Model3D* model)` — `GS3Dp.cc:1245`
- 🔧 `int loadModel3D(Model3D* model, const char* filename)` — `GS3Dp.cc:1306`
- 🔧 `void computeModelBoundingSphere(Model3D* model)` — `GS3Dp.cc:1339`
- 🔧 `void getObserverParams(ObserverParams* params, Model3D* model)` — `GS3Dp.cc:1402`
- 🔧 `void processModelFast(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:1498`
- 🔧 `void processModelWireframe(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:1603`
- 🔧 `int readVertices(const char* filename, VertexArrays3D* vtx, int max_vertices)` — `GS3Dp.cc:1706`
- 🔧 `int readFaces_model(const char* filename, Model3D* model)` — `GS3Dp.cc:1760`
- 🔧 `void projectTo2D(VertexArrays3D* vtx, int angle_w_deg)` — `GS3Dp.cc:1875`
- 🔧 `void calculateFaceDepths(Model3D* model, Face3D* faces, int face_count)` — `GS3Dp.cc:1940`
- 🔧 `void autoScaleModel(Model3D* model, float target_max_dim, float min_scale, float max_scale, int center_flag)` — `GS3Dp.cc:2112`
- 🔧 `void revertAutoScaleModel(Model3D* model)` — `GS3Dp.cc:2177`
- 🔧 `void backupModelCoords(Model3D* model)` — `GS3Dp.cc:2207`
- 🔧 `void freeBackupModelCoords(Model3D* model)` — `GS3Dp.cc:2234`
- 🔧 `void fitModelToView(Model3D* model, ObserverParams* params, float target_max_dim, float margin, float percentile, int center_flag)` — `GS3Dp.cc:2243`
- 🔧 `void drawPolygons(Model3D* model, int* vertex_count, int face_count, int vertex_count_total)` — `GS3Dp.cc:2516`
- 🔧 `int main()` — `GS3Dp.cc:2689`

---

*Generated from `call_flow.txt` and converted to Markdown for readability.*
