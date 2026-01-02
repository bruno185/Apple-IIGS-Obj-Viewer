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
    - 🔧 **`calculateFaceDepths()`** — compute per-face `z_min/z_max/z_mean`, display flags, planar coefficients (Newell). Optionally performs observer-space back-face culling (plane D <= 0) when the `B` toggle is enabled.
    - 🔧 **`painter_newell_sancha()`** — sort faces by depth and correct ambiguous order (qsort + corrections)
      - Collects *inconclusive pairs* (pairs of faces where order is ambiguous) into an **in-memory buffer** for later inspection; note: the buffer `inconclusive_pairs` is now a global buffer (preallocated for performance) used by diagnostic and framing helpers.

    - 🔧 **`processModelWireframe(model, &params, filename)`** — lightweight wireframe processing: transforms and projects vertices and sets face visibility only (no per-face depth/sorting); used for wireframe/frame-only display to improve speed.

- 🔧 **`drawPolygons(model, faces, face_count, vert_count)`** — render loop
  - uses sorted faces; for each face builds QuickDraw polygon from `vtx->x2d/y2d` (these already reflect any model-space auto-scaling applied by `fitModelToView`)
  - **Allocates QuickDraw polygon handle sized using `MAX_FACE_VERTICES` to support polygons up to that many vertices (default 6).**
  - Fill + Frame polygon (QuickDraw)
  - There is also a helper `frameInconclusivePairs()` that can frame (in white) polygons listed in the `inconclusive_pairs` buffer for debugging/diagnostic display.

### UI / Input handling

- `startgraph()` / render / `endgraph()` / `DoText()` / optional `DoColor()`
- Keys: Space (info), N (new model), Arrows / A Z (angles/distance), `K` (edit angles/distance without reloading model), `+`/`-` (adjust autoscale), `B` (toggle back-face culling: observer-space d<=0 test)
- `K` invokes `getObserverParams(&params, model)` interactively and applies new angles/distance without requiring a reload.
- `+`/`-` behavior: if the model is not yet auto-scaled, these keys first perform an **auto-fit** (`fitModelToView()`); they then increase or decrease the current `params->distance` and update `model->auto_scale`. Automatic recomputation via bounding-sphere is disabled; distance adjustments are manual. The action prints a short message (e.g., "Distance increased" / "Distance decreased").

---

## Supporting / fallback functions (short purpose)


- Bounding-sphere metric removed; auto-fit uses bbox heuristic (O(n) at load remains but is simpler)

- 🔧 **Non-destructive backup support removed** — per-vertex backup API was removed to simplify flow
- 🔧 **`processModelWireframe(model, &params, filename)`** — lightweight wireframe processing (transform & project only) used for fast wireframe/frame-only rendering
- 🔧 **`destroyModel3D(Model3D* model)`** — frees all memory allocated by `createModel3D()`; must be called to avoid leaks
- 🔧 **`readVertices()` / `readFaces_model()`** — file parsing helpers
- 🔧 **`frameInconclusivePairs(Model3D* model)`** — utility: frames in white all polygons currently recorded in the global `inconclusive_pairs` buffer (diagnostic; no runtime side effects beyond rendering) 

**Notes:** `adjustDistanceFast()` (earlier fast-adjust prototype) was removed — automatic distance estimation has been disabled; distance adjustments are manual.

---

## Notes

- This Markdown preserves the visual grouping: **function name** + **purpose/action**, and uses `→` to show call flow.
- If you want more detail (e.g., exact file:line references or prototypes), I can append them in a table or add a per-function section with signatures.

---

## Function signatures & file:line references

### Core / Public APIs
- 🔧 `Model3D* createModel3D(void)` — `GS3Dp.cc:1514`
- 🔧 `void destroyModel3D(Model3D* model)` — `GS3Dp.cc:1835`
- 🔧 `int loadModel3D(Model3D* model, const char* filename)` — `GS3Dp.cc:1895`
- 🔧 `int readVertices(const char* filename, VertexArrays3D* vtx, int max_vertices, Model3D* owner)` — `GS3Dp.cc:2262`
- 🔧 `int readFaces_model(const char* filename, Model3D* model)` — `GS3Dp.cc:2392`

### Rendering / Painter
- 🔧 `void processModelFast(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:2037`
- 🔧 `void processModelWireframe(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:2152`
- 🔧 `void calculateFaceDepths(Model3D* model, Face3D* faces, int face_count)` — `GS3Dp.cc:2546`
- 🔧 `void painter_newell_sancha(Model3D* model, int face_count)` — `GS3Dp.cc:703`
- 🔧 `void painter_newell_sancha_fast(Model3D* model, int face_count)` — `GS3Dp.cc:681`
- 🔧 `void painter_newell_sancha_float(Model3D* model, int face_count)` — `GS3Dp.cc:1178`
- 🔧 `void drawPolygons(Model3D* model, int* vertex_count, int face_count, int vertex_count_total)` — `GS3Dp.cc:2956`
- 🔧 `void drawFace(Model3D* model, int face_id)` — `GS3Dp.cc:2848`
- 🔧 `void frameInconclusivePairs(Model3D* model)` — `GS3Dp.cc:3093`

### Utilities / Helpers
- 🔧 `void fitModelToView(Model3D* model, ObserverParams* params, float target_max_dim, float margin, float percentile, int center_flag)` — `GS3Dp.cc:2804`
- 🔧 `void dumpFaceEquationsCSV(Model3D* model, const char* csv_filename, int alt_format)` — `GS3Dp.cc:2681`
- 🔧 `void getObserverParams(ObserverParams* params, Model3D* model)` — `GS3Dp.cc:1947`
- 🔧 `void compute2DFromObserver(Model3D* model, int angle_w)` — `GS3Dp.cc:3192`
- 🔧 `void DoColor()` — `GS3Dp.cc:3219`
- 🔧 `void DoText()` — `GS3Dp.cc:3248`

### Internal / Small helpers (static/inline)
- 🔧 `static void ensure_vertex_capacity(int vcount)` — `GS3Dp.cc:97`
- 🔧 `static void ensure_face_capacity(int face_count)` — `GS3Dp.cc:109`
- 🔧 `static void ensure_order_capacity(int face_count)` — `GS3Dp.cc:131`
- 🔧 `static inline Fixed32 sin_deg_int(int deg)` — `GS3Dp.cc:250`
- 🔧 `static inline Fixed32 cos_deg_int(int deg)` — `GS3Dp.cc:255`
- 🔧 `static inline int FIXED_ROUND_TO_INT(Fixed32 x)` — `GS3Dp.cc:266`
- 🔧 `static inline int normalize_deg(int deg)` — `GS3Dp.cc:295`
- 🔧 `static int cmp_faces_by_zmean(const void* pa, const void* pb)` — `GS3Dp.cc:653`

---

*All functions listed are present in `GS3Dp.cc` as of this commit; helper functions are grouped separately.*

---

*Generated from `call_flow.txt` and converted to Markdown for readability.*
