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
      - set `params->distance = computeDistanceFromBoundingSphere(model, margin)` (fast **O(1)**)
      - **FALLBACK:** 🔧 **`computeDistanceToFit()`** (slower per-vertex) only if sphere not valid
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

- 🔧 **`drawPolygons(model, faces, face_count, vert_count)`** — render loop
  - uses sorted faces; for each face builds QuickDraw polygon from `vtx->x2d/y2d` (these already reflect any model-space auto-scaling applied by `fitModelToView`)
  - Fill + Frame polygon (QuickDraw)

### UI / Input handling

- `startgraph()` / render / `endgraph()` / `DoText()` / optional `DoColor()`
- Keys: Space (info), N (new model), Arrows / A Z (angles/distance), `r` (revert autoscale), `K` (edit angles/distance without reloading model), `+`/`-` (adjust autoscale)
- `K` invokes `getObserverParams(&params, model)` interactively and applies new angles/distance without requiring a reload.
- `+`/`-` behavior: if the model is not yet auto-scaled, these keys first perform an **auto-fit** (`fitModelToView()`); they then increase or decrease the current `params->distance`, update `model->auto_scale` and **scale the bounding sphere proportionally** (O(1)), and recompute `params->distance` using `computeDistanceFromBoundingSphere()` (fast). The action prints a short message (e.g., "Distance increased" / "Distance decreased").

---

## Supporting / fallback functions (short purpose)

- 🔧 **`projectTo2D()`** — standalone projection helper (used in places)
- 🔧 **`computeDistanceToFit(vtx, margin)`** — slower per-vertex fit (kept as fallback / DEPRECATED)
- 🔧 **`computeModelBoundingSphere(model)`** — O(n) at load, stores `bs_*` in `Model3D`
- 🔧 **`computeDistanceFromBoundingSphere(model, margin)`** — O(1) camera distance estimate using the bounding sphere
- 🔧 **`autoScaleModel()` / `revertAutoScaleModel()`** — non-destructive scaling helpers (backup + apply + revert)
- 🔧 **`backupModelCoords()` / `freeBackupModelCoords()`** — support for non-destructive transforms
- 🔧 **`readVertices()` / `readFaces_model()`** — file parsing helpers

**Notes:** `adjustDistanceFast()` (earlier fast-adjust prototype) was removed — use `+`/`-` behavior and `computeDistanceFromBoundingSphere()` for distance adjustments.

---

## Notes

- This Markdown preserves the visual grouping: **function name** + **purpose/action**, and uses `→` to show call flow.
- If you want more detail (e.g., exact file:line references or prototypes), I can append them in a table or add a per-function section with signatures.

---

## Function signatures & file:line references

- 🔧 `int main()` — `GS3Dp.cc:2637`
- 🔧 `Model3D* createModel3D(void)` — `GS3Dp.cc:939`
- 🔧 `int loadModel3D(Model3D* model, const char* filename)` — `GS3Dp.cc:1318`
- 🔧 `void computeModelBoundingSphere(Model3D* model)` — `GS3Dp.cc:1351`
- 🔧 `Fixed32 computeDistanceFromBoundingSphere(Model3D* model, float margin)` — `GS3Dp.cc:1379`
- 🔧 `void getObserverParams(ObserverParams* params, Model3D* model)` — `GS3Dp.cc:1414`
- 🔧 `void fitModelToView(Model3D* model, ObserverParams* params, float target_max_dim, float margin, float percentile, int center_flag)` — `GS3Dp.cc:2198`
- 🔧 `void processModelFast(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:1511`
- 🔧 `void drawPolygons(Model3D* model, int* vertex_count, int face_count, int vertex_count_total)` — `GS3Dp.cc:2471`
- 🔧 `Fixed32 computeDistanceToFit(VertexArrays3D* vtx, float margin)` — `GS3Dp.cc:2017`

---

*Generated from `call_flow.txt` and converted to Markdown for readability.*
