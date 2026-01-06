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

### Parameter parsing / Auto-fit

- 🔧 **`getObserverParams(&params, model)`** — interactive parameter parsing
  - reads angles H/V/W and screen rotation
  - If user presses ENTER for distance: a precomputed auto-fit suggestion may be applied if available; on-the-fly auto-fit is disabled.
  - else: set `params->distance` from user input

---

## Main render loop

- Enter main render loop (`bigloop`):
  - 🔧 **`processModelFast(model, &params, filename)`** — runs every frame; ultra-fast transformation + projection
    - precompute trig products (Fixed32)
    - For each vertex (tight Fixed32 loop): transform → compute `xo/yo/zo` → project to `x2d/y2d`
    - **Note:** runtime auto-fit is *not* applied inside `processModelFast` — any autoscale must be applied ahead of time; the auto-fit helper has been removed and archived to `chutier.txt`. `processModelFast` operates on already-scaled model vertices.
    - 🔧 **`calculateFaceDepths()`** — compute per-face `z_min/z_max/z_mean`, display flags, planar coefficients (Newell). Optionally performs observer-space back-face culling (plane D <= 0) when the `B` toggle is enabled.
    - 🔧 **`painter_newell_sancha()`** — sort faces by depth and correct ambiguous order (qsort + corrections). When back-face culling is enabled, the painter builds and sorts a list limited to faces with `display_flag == 1` (visible faces), performs order corrections only on that sub-list for efficiency and correctness, and appends culled faces afterward to preserve `sorted_face_indices` stability.
      - Collects *inconclusive pairs* (pairs of faces where order is ambiguous) into an **in-memory buffer** for later inspection; note: the buffer `inconclusive_pairs` is now a global buffer (preallocated for performance) used by diagnostic and framing helpers.

    - ⚠️ **Wireframe mode** — handled via the `framePolyOnly` flag and `drawPolygons()` (no separate `processModelWireframe()` function in the current codebase). Use `framePolyOnly=1` to render wireframe previews.

- 🔧 **`drawPolygons(model, faces, face_count, vert_count)`** — render loop
  - uses sorted faces; for each face builds QuickDraw polygon from `vtx->x2d/y2d` (these already reflect any model-space auto-scaling applied ahead of time; auto-fit is disabled — archived implementation is in `chutier.txt`). When back-face culling is active, `faces->sorted_face_indices` contains visible faces first (sorted) followed by culled faces; `drawPolygons` still checks `display_flag` and skips any face with `display_flag == 0` at draw time.
  - Fill + Frame polygon (QuickDraw)
  - There is also a helper `frameInconclusivePairs()` that can frame (in white) polygons listed in the `inconclusive_pairs` buffer for debugging/diagnostic display.

### UI / Input handling

- `startgraph()` / render / `endgraph()` / `DoText()` / optional `DoColor()`
- Keys: Space (info), N (new model), Arrows / A Z (angles/distance), `K` (edit angles/distance without reloading model), `+`/`-` (adjust autoscale), `B` (toggle back-face culling: observer-space d<=0 test)
  - **New keys (recent additions):**
    - `D` / `d`: Inspect faces placed BEFORE a selected face in the painter order and report misplaced faces (previews in orange)
    - `S` / `s`: Inspect faces placed AFTER a selected face that should be BEFORE it (previews in pink)
    - `O` / `o`: Interactive **projected polygon overlap** inspector — prompts for two face IDs, reports YES/NO if their 2D projections overlap, and optionally previews the two faces (green/orange)
    - `L` / `l`: Label mode — show the model with each face's ID drawn at its polygon center
- `K` invokes `getObserverParams(&params, model)` interactively and applies new angles/distance without requiring a reload.
- `+`/`-` behavior: if the model is not yet auto-scaled, these keys previously performed an **auto-fit** via `fitModelToView()`; auto-fit has been disabled (see `chutier.txt` for the archived implementation). They now only increase or decrease the current `params->distance` and update `model->auto_scale`. Automatic recomputation via bounding-sphere is disabled; distance adjustments are manual. The action prints a short message (e.g., "Distance increased" / "Distance decreased").

---

## Supporting / fallback functions (short purpose)


- Bounding-sphere metric removed; auto-fit uses bbox heuristic (O(n) at load remains but is simpler)

- 🔧 **Non-destructive backup support removed** — per-vertex backup API was removed to simplify flow
- 🔧 **Wireframe mode** — implemented via the `framePolyOnly` flag and `drawPolygons()` (no separate `processModelWireframe()` function)
- 🔧 **`destroyModel3D(Model3D* model)`** — frees all memory allocated by `createModel3D()`; must be called to avoid leaks
- 🔧 **`readVertices()` / `readFaces_model()`** — file parsing helpers
- 🔧 **`frameInconclusivePairs(Model3D* model)`** — utility: frames in white all polygons currently recorded in the global `inconclusive_pairs` buffer (diagnostic; no runtime side effects beyond rendering) 
- 🔧 **`projected_polygons_overlap(Model3D* model, int f1, int f2)`** — screen-space test that returns 1 if two faces' projected 2D polygons *overlap* (proper edge intersection or containment), **0 if disjoint**. Important: *touching-only* cases (shared edge or single-vertex contact) are considered **NON-overlap** and return 0. The algorithm uses integer segment intersection (proper intersection only) then ray-casting containment; points on edges are treated as outside.
- 🔧 `void inspect_faces_before(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:2032` — interactive wrapper bound to `D`/`d`: prints a compact per-face diagnostic and offers a wireframe preview that highlights faces placed BEFORE a selected face (misplaced faces highlighted).
- 🔧 `void inspect_faces_after(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:2219` — interactive wrapper bound to `S`/`s`: prints faces placed AFTER a selected face that should be BEFORE it and offers a wireframe preview with highlights.
- 🔧 `void inspect_polygons_overlap(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:2371` — interactive wrapper bound to `O`/`o`: prompts for two face IDs, reports overlap status (YES/NO), and optionally previews the two faces (green/orange); default on ENTER shows the entire model in wireframe with the two faces highlighted.
- 🔧 `void display_model_face_ids(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:2449` — label mode bound to `L`/`l`: draws the model in wireframe and overlays each face's ID centered on that face (uses `drawFace(..., show_index=1)`), useful for debugging face ordering and references.

**Notes:** Automatic distance estimation has been disabled; distance adjustments are manual.

---

## Notes

- This Markdown preserves the visual grouping: **function name** + **purpose/action**, and uses `→` to show call flow.
- If you want more detail (e.g., exact file:line references or prototypes), I can append them in a table or add a per-function section with signatures.

---

## Function signatures & file:line references

### Core / Public APIs
- 🔧 `Model3D* createModel3D(void)` — `GS3Dp.cc:2540`
- 🔧 `void destroyModel3D(Model3D* model)` — `GS3Dp.cc:2860`
- 🔧 `int loadModel3D(Model3D* model, const char* filename)` — `GS3Dp.cc:2920`
- 🔧 `int readVertices(const char* filename, VertexArrays3D* vtx, int max_vertices, Model3D* owner)` — `GS3Dp.cc:3220`
- 🔧 `int readFaces_model(const char* filename, Model3D* model)` — `GS3Dp.cc:3351`

### Rendering / Painter
- 🔧 `void processModelFast(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:3077`
- 🔧 `void calculateFaceDepths(Model3D* model, Face3D* faces, int face_count)` — `GS3Dp.cc:3517`
- 🔧 `void painter_newell_sancha_fast(Model3D* model, int face_count)` — `GS3Dp.cc:732`
- 🔧 `void painter_newell_sancha(Model3D* model, int face_count)` — `GS3Dp.cc:772`
- 🔧 `void painter_newell_sancha_float(Model3D* model, int face_count)` — `GS3Dp.cc:1325`
- 🔧 `void drawFace(Model3D* model, int face_id, int fillPenPat, int show_index)` — `GS3Dp.cc:3794`
- 🔧 `void drawPolygons(Model3D* model, int* vertex_count, int face_count, int vertex_count_total)` — `GS3Dp.cc:3913`
- 🔧 `void frameInconclusivePairs(Model3D* model)` — `GS3Dp.cc:4055`

> Note: Wireframe mode is handled by the `framePolyOnly` flag and using `drawPolygons()`; there is no separate `processModelWireframe()` implementation in the current file.

### Utilities / Helpers
- 🔧 `void dumpFaceEquationsCSV(Model3D* model, const char* csv_filename, int alt_format)` — `GS3Dp.cc:3671`
- 🔧 `void debug_two_faces(Model3D* model, int f1, int f2)` — `GS3Dp.cc:762` — small debug helper to preview two faces side-by-side (used by painter diagnostics)
- 🔧 `void getObserverParams(ObserverParams* params, Model3D* model)` — `GS3Dp.cc:2971`
- 🔧 `void compute2DFromObserver(Model3D* model, int angle_w)` — `GS3Dp.cc:4155`
- 🔧 `void DoColor()` — `GS3Dp.cc:4183`
- 🔧 `void DoText()` — `GS3Dp.cc:4213`
- 🔧 `void inspect_faces_before(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:2032`
- 🔧 `void inspect_faces_after(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:2219`
- 🔧 `void inspect_polygons_overlap(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:2371`
- 🔧 `void display_model_face_ids(Model3D* model, ObserverParams* params, const char* filename)` — `GS3Dp.cc:2441`

### Internal / Small helpers (static/inline)
- 🔧 `static void ensure_vertex_capacity(int vcount)` — `GS3Dp.cc:98`
- 🔧 `static void ensure_face_capacity(int face_count)` — `GS3Dp.cc:110`
- 🔧 `static void ensure_order_capacity(int face_count)` — `GS3Dp.cc:132`
- 🔧 `static inline Fixed32 sin_deg_int(int deg)` — `GS3Dp.cc:261`
- 🔧 `static inline Fixed32 cos_deg_int(int deg)` — `GS3Dp.cc:266`
- 🔧 `static inline int FIXED_ROUND_TO_INT(Fixed32 x)` — `GS3Dp.cc:277`
- 🔧 `int prepare_inspector_sort(Model3D* m, int fc)` — `GS3Dp.cc:2040` — internal helper used by `inspect_faces_before/after` to prepare sorting lists.
- 🔧 `static inline int normalize_deg(int deg)` — `GS3Dp.cc:303`
- 🔧 `static int cmp_faces_by_zmean(const void* pa, const void* pb)` — `GS3Dp.cc:673`

---

*All functions listed are present in `GS3Dp.cc` as of this commit; helper functions are grouped separately.*

---

*Generated from `call_flow.txt` and converted to Markdown for readability.*
