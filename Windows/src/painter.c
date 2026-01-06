#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <string.h>
#include <windows.h> /* for polygon region helpers used by V4 */

// Minimal structures (mirrors those in obj.c)
typedef struct { float x,y,z; } Vertex;
// Extended Face struct to match GS3Dp: add plane coeffs, bbox and display flag
typedef struct { int *indices; int count; float z_min,z_mean,z_max; float plane_a, plane_b, plane_c, plane_d; long long plane_a_i, plane_b_i, plane_c_i, plane_d_i; int minx, maxx, miny, maxy; int display_flag; } Face;
typedef struct { Vertex* verts; int vert_count; Face* faces; int face_count; } Model;

// Observer params (match GS3Dp semantics: angle_h, angle_v, angle_w, distance)
static int s_angle_h = 30; static int s_angle_v = 20; static int s_angle_w = 0; static float s_distance = 300.0f;

void set_observer_params(float ah, float av, float aw, float dist) { s_angle_h = (int)roundf(ah); s_angle_v = (int)roundf(av); s_angle_w = (int)roundf(aw); s_distance = dist; }

// Projection parameters set by win32_main to ensure bbox uses same cx/cy/scale
static float s_proj_cx = 0.0f, s_proj_cy = 0.0f, s_proj_scale = 200.0f;
void set_projection_params(float cx, float cy, float scale) { s_proj_cx = cx; s_proj_cy = cy; s_proj_scale = scale; }

// Transform model into observer (GS3Dp float port of processModelFast per-vertex transform)
typedef struct { float xo, yo, zo; } ObsVertex;

/* Forward-declare V4 helpers so they can be used earlier for arrondi when computing planes */
static float V4_arrondi(float r);

// Exposed: compute observer-space vertices for rendering / comparison
void compute_obs_vertices(Model* m, ObsVertex* out) {
    // Port of GS3Dp.cc per-vertex transform using float trig (angles in degrees)
    float cos_h = cosf(s_angle_h * (M_PI/180.0f));
    float sin_h = sinf(s_angle_h * (M_PI/180.0f));
    float cos_v = cosf(s_angle_v * (M_PI/180.0f));
    float sin_v = sinf(s_angle_v * (M_PI/180.0f));
    float cos_h_cos_v = cos_h * cos_v;
    float sin_h_cos_v = sin_h * cos_v;
    float cos_h_sin_v = cos_h * sin_v;
    float sin_h_sin_v = sin_h * sin_v;
    /* Pascal multiplies distance by 4 in the fixed pipeline; replicate that here for parity */
    float distance = s_distance * 4.0f;
    /* Emulate Pascal Fixed32 (16.16) pipeline exactly to reproduce integer rounding/shift
       behavior used by TestComplet. Compute trig values in fixed format once and use
       64-bit intermediate multiplies with >>16 shifts (fixed_mul_64 semantics). */
    int cos_h_f = (int)lroundf(cos_h * 65536.0f);
    int sin_h_f = (int)lroundf(sin_h * 65536.0f);
    int cos_v_f = (int)lroundf(cos_v * 65536.0f);
    int sin_v_f = (int)lroundf(sin_v * 65536.0f);
    int cos_h_cos_v_f = (int)(((long long)cos_h_f * (long long)cos_v_f) >> 16);
    int sin_h_cos_v_f = (int)(((long long)sin_h_f * (long long)cos_v_f) >> 16);
    int cos_h_sin_v_f = (int)(((long long)cos_h_f * (long long)sin_v_f) >> 16);
    int sin_h_sin_v_f = (int)(((long long)sin_h_f * (long long)sin_v_f) >> 16);
    int distance_f = (int)lroundf(distance * 65536.0f);

    for (int i = 0; i < m->vert_count; ++i) {
        float xf = m->verts[i].x;
        float yf = m->verts[i].y;
        float zf = m->verts[i].z;
        /* To match Pascal we convert coordinates to 16.16 ints and perform fixed ops */
        int x_f = (int)lroundf(xf * 65536.0f);
        int y_f = (int)lroundf(yf * 65536.0f);
        int z_f = (int)lroundf(zf * 65536.0f);
        long long term1_f = ((long long)x_f * (long long)cos_h_cos_v_f) >> 16;
        long long term2_f = ((long long)y_f * (long long)sin_h_cos_v_f) >> 16;
        long long term3_f = ((long long)z_f * (long long)sin_v_f) >> 16;
        long long zo_f = -term1_f - term2_f - term3_f + (long long)distance_f;
        out[i].zo = (float)((double)zo_f / 65536.0);
        if (zo_f > 0) {
            long long xo_f = -(((long long)x_f * (long long)sin_h_f) >> 16) + (((long long)y_f * (long long)cos_h_f) >> 16);
            long long yo_f = -(((long long)x_f * (long long)cos_h_sin_v_f) >> 16) - (((long long)y_f * (long long)sin_h_sin_v_f) >> 16) + (((long long)z_f * (long long)cos_v_f) >> 16);
            out[i].xo = (float)((double)xo_f / 65536.0);
            out[i].yo = (float)((double)yo_f / 65536.0);
        } else {
            out[i].xo = 0.0f; out[i].yo = 0.0f;
        }
    }
}

// Project using GS3Dp convention: screen x = xo/zo, y = yo/zo (zo is depth)
static void project_vertex(ObsVertex* v, float* px, float* py) {
    if (v->zo == 0.0f) v->zo = 1e-6f;
    *px = v->xo / v->zo; *py = v->yo / v->zo;
}

// Provide a small screen-coordinate helper for the Win32 build (renderer.c is not part of WIN32 target)
void screen_coords_from_proj(float px, float py, int winw, int winh, float scale, float centerx, float centery, int* sx, int* sy) {
    *sx = (int)(winw*0.5f + (px - centerx) * scale);
    *sy = (int)(winh*0.5f - (py - centery) * scale);
}

// Forward declarations used by the inspector (avoid implicit-int conflicts)
int projected_polygons_overlap(Model* m, int f1, int f2);
void set_highlight_faces(int* faces, int count, int color);

// Simple interactive inspector wrapper (console) to match GS3Dp's inspect_polygons_overlap behaviour
void inspect_polygons_overlap(Model* m, void* params_unused, const char* filename) {
    if (!m) { printf("No model loaded\n"); return; }
    printf("Enter two face IDs (f1 f2) or empty to cancel: "); char buf[128]; if (fgets(buf, sizeof(buf), stdin)) {
        int f1=-1,f2=-1; if (sscanf(buf, "%d %d", &f1, &f2) == 2) {
            int ov = projected_polygons_overlap(m, f1, f2);
            printf("Projected overlap: %s\n", ov ? "YES" : "NO");
            if (ov) {
                printf("Press Y to preview highlighting faces, any other key to skip: "); int c = getchar(); if (c=='Y' || c=='y') { int arr[2] = {f1,f2}; set_highlight_faces(arr, 2, 10); }
            }
        } else { printf("Invalid input\n"); }
    }
}

/* per-instr trace file pointer (declared before helper so helper can reference it) */
static FILE* g_v4_pvf = NULL; /* optional per-instruction trace file (set by V4 when opened) */

/* Per-instruction tracer helper (no-op if g_v4_pvf is NULL). Lines are simple CSV for easy parsing.
   Format: instr,kind,pass,i,j,f1,f2,stage,idx,expr,subexpr, value
*/
static void V4_instr_log(const char* kind, int pass, int i, int j, int f1, int f2, const char* stage, int idx, const char* expr, const char* subexpr, float value) {
    if (!g_v4_pvf) return;
    fprintf(g_v4_pvf, "instr,%s,%d,%d,%d,%d,%d,%s,%d,%s,%s,%.9f\n", kind, pass, i, j, f1, f2, stage, idx, expr, subexpr, value);
}

// Pair comparator helpers and GS3Dp-like calculateFaceDepths + insertion-based Newell/Sancha
static Model* g_model = NULL; static ObsVertex* g_obsv = NULL; // for comparator
/* Runtime selection: 1 = original adjacent-swap (V1), 2 = pairwise compare (V2) */
int g_painter_order_version = 1;

// Painter mode selection to match GS3Dp (FAST/FIXED/FLOAT)
#ifndef PAINTER_MODE_FAST
#define PAINTER_MODE_FAST 0
#define PAINTER_MODE_FIXED 1
#define PAINTER_MODE_FLOAT 2
#endif
int g_painter_mode = 0; /* default = PAINTER_MODE_FAST */
void set_painter_mode(int mode) {
    if (mode < 0 || mode > 2) return;
    g_painter_mode = mode;
    // Map painter mode to internal ordering variant for parity with GS3Dp
    if (g_painter_mode == PAINTER_MODE_FAST) g_painter_order_version = 1;
    else if (g_painter_mode == PAINTER_MODE_FIXED) g_painter_order_version = 2;
    else if (g_painter_mode == PAINTER_MODE_FLOAT) g_painter_order_version = 4;

    char buf[128]; snprintf(buf, sizeof(buf), "Painter mode set to %d, order_version=%d\r\n", g_painter_mode, g_painter_order_version);
    // write small log for diagnostics
    FILE* lf = fopen("viewer_win32.log", "a"); if (lf) { fprintf(lf, "%s", buf); fclose(lf); }

    // If a model is loaded, recompute observer vertices and painter order immediately
    if (g_model && g_obsv) {
        compute_obs_vertices(g_model, g_obsv);
        // recompute order using the new ordering selection
        // allocate an order buffer on stack if g_model->face_count is reasonable
        int *ord = (int*)malloc(sizeof(int) * g_model->face_count);
        if (ord) {
            compute_painter_order(g_model, ord);
            // If there is a global order buffer in the UI, it will be re-computed by caller; here we just leave recomputed result
            free(ord);
        }
    }
}

/* Runtime back-face culling toggle - default enabled (match GS3Dp behavior)
   Use set_cull_back_faces/get_cull_back_faces to toggle at runtime from the UI */
static int s_cull_back_faces = 1;
void set_cull_back_faces(int v) { s_cull_back_faces = v ? 1 : 0; }
int get_cull_back_faces(void) { return s_cull_back_faces; }

/* Highlighting helpers (used by inspectors to preview faces) */
static int* g_highlight_faces = NULL; static int g_highlight_count = 0; static int g_highlight_color = 0;
void set_highlight_faces(int* faces, int count, int color) {
    if (g_highlight_faces) free(g_highlight_faces);
    if (count <= 0) { g_highlight_faces = NULL; g_highlight_count = 0; g_highlight_color = 0; return; }
    g_highlight_faces = (int*)malloc(sizeof(int) * count);
    memcpy(g_highlight_faces, faces, sizeof(int) * count);
    g_highlight_count = count; g_highlight_color = color;
}
void clear_highlight_faces(void) { if (g_highlight_faces) free(g_highlight_faces); g_highlight_faces = NULL; g_highlight_count = 0; g_highlight_color = 0; }
int face_is_highlighted(int f) { for (int i=0;i<g_highlight_count;i++) if (g_highlight_faces[i]==f) return g_highlight_color; return 0; }

/* Inconclusive pair buffer: collect pairs of faces where the painter tests are inconclusive
   so they can be highlighted later for diagnostics. Stored as pairs [f1,f2, f1,f2, ...] */
static int* g_inconclusive_pairs = NULL; static int g_inconclusive_count = 0; static int g_inconclusive_capacity = 0;

void clear_inconclusive_pairs(void) {
    if (g_inconclusive_pairs) free(g_inconclusive_pairs);
    g_inconclusive_pairs = NULL; g_inconclusive_count = 0; g_inconclusive_capacity = 0;
}

void add_inconclusive_pair(int f1, int f2) {
    if (f1 < 0 || f2 < 0) return;
    // Avoid duplicates (unordered pair)
    for (int i = 0; i < g_inconclusive_count; ++i) {
        int a = g_inconclusive_pairs[2*i], b = g_inconclusive_pairs[2*i+1];
        if ((a==f1 && b==f2) || (a==f2 && b==f1)) return; // already recorded
    }
    if (g_inconclusive_count + 1 > g_inconclusive_capacity) {
        int newcap = (g_inconclusive_capacity == 0) ? 8 : g_inconclusive_capacity * 2;
        g_inconclusive_pairs = (int*)realloc(g_inconclusive_pairs, sizeof(int) * 2 * newcap);
        g_inconclusive_capacity = newcap;
    }
    g_inconclusive_pairs[2 * g_inconclusive_count] = f1;
    g_inconclusive_pairs[2 * g_inconclusive_count + 1] = f2;
    g_inconclusive_count++;
}

int get_inconclusive_pair_count(void) { return g_inconclusive_count; }
int* get_inconclusive_pairs(void) { return g_inconclusive_pairs; }

/* point_in_triangle removed to match ORCA/GS3Dp (no local per-triangle depth tests there) */

// Quick comparator for initial qsort by z_mean descending (stable tie-breaker by index)
static int compar_face_qsort(const void* pa, const void* pb) {
    int a = *(const int*)pa; int b = *(const int*)pb;
    float za = g_model->faces[a].z_mean; float zb = g_model->faces[b].z_mean;
    if (za > zb) return -1; if (za < zb) return 1; if (a < b) return -1; if (a > b) return 1; return 0;
}

/* Helper: determine whether face f1 should be before face f2 according to
   the same Tests 4..7 used by the painter. Returns 1 if f1 should be before f2,
   -1 if f2 should be before f1, 0 if inconclusive. */
static int face_should_be_before(Model* m, int f1, int f2) {
    if (!m) return 0;
    // quick z-range checks
    if (m->faces[f2].z_max <= m->faces[f1].z_min) return 1;   // f1 in front (no swap needed)
    if (m->faces[f1].z_max <= m->faces[f2].z_min) return -1;  // f2 in front (swap)
    // bbox quick rejection
    if (m->faces[f1].maxx <= m->faces[f2].minx || m->faces[f2].maxx <= m->faces[f1].minx) return 0;
    if (m->faces[f1].maxy <= m->faces[f2].miny || m->faces[f2].maxy <= m->faces[f1].miny) return 0;
    // If FIXED mode is selected, use integer Fixed arithmetic identical to GS3Dp for Tests 4..7
    if (g_painter_mode == PAINTER_MODE_FIXED) {
        int n1 = m->faces[f1].count, n2 = m->faces[f2].count;
        long long a1_i = m->faces[f1].plane_a_i, b1_i = m->faces[f1].plane_b_i, c1_i = m->faces[f1].plane_c_i, d1_i = m->faces[f1].plane_d_i;
        long long a2_i = m->faces[f2].plane_a_i, b2_i = m->faces[f2].plane_b_i, c2_i = m->faces[f2].plane_c_i, d2_i = m->faces[f2].plane_d_i;
        // Fixed epsilon = FLOAT_TO_FIXED(0.01f)
        int fixed_eps = (int)lroundf(0.01f * 65536.0f); long long eps64 = (long long)fixed_eps;
        int k;
        int obs_side1 = 0, obs_side2 = 0, side = 0, all_same_side = 0, all_opposite_side = 0;
        // Test 4 (all same side -> f2 before f1)
        obs_side1 = 0; if (d1_i > eps64) obs_side1 = 1; else if (d1_i < -eps64) obs_side1 = -1; else goto skipT4_fixed;
        all_same_side = 1;
        for (k = 0; k < n2; ++k) {
            int v = m->faces[f2].indices[k]; int xo_i = (int)lroundf(g_obsv[v].xo * 65536.0f); int yo_i = (int)lroundf(g_obsv[v].yo * 65536.0f); int zo_i = (int)lroundf(g_obsv[v].zo * 65536.0f);
            long long acc = 0;
            acc  = (((long long)a1_i * (long long)xo_i) >> 16);
            acc += (((long long)b1_i * (long long)yo_i) >> 16);
            acc += (((long long)c1_i * (long long)zo_i) >> 16);
            acc += (long long)d1_i;
            if (acc > eps64) side = 1; else if (acc < -eps64) side = -1; else continue;
            if (obs_side1 != side) { all_same_side = 0; break; }
        }
        if (all_same_side) return -1;
    skipT4_fixed: ;
        // Test 5 (all opposite side -> f2 before f1)
        obs_side2 = 0; if (d2_i > eps64) obs_side2 = 1; else if (d2_i < -eps64) obs_side2 = -1; else goto skipT5_fixed;
        all_opposite_side = 1;
        for (k = 0; k < n1; ++k) {
            int v = m->faces[f1].indices[k]; int xo_i = (int)lroundf(g_obsv[v].xo * 65536.0f); int yo_i = (int)lroundf(g_obsv[v].yo * 65536.0f); int zo_i = (int)lroundf(g_obsv[v].zo * 65536.0f);
            long long acc = 0;
            acc  = (((long long)a2_i * (long long)xo_i) >> 16);
            acc += (((long long)b2_i * (long long)yo_i) >> 16);
            acc += (((long long)c2_i * (long long)zo_i) >> 16);
            acc += (long long)d2_i;
            if (acc > eps64) side = 1; else if (acc < -eps64) side = -1; else continue;
            if (obs_side2 == side) { all_opposite_side = 0; break; }
        }
        if (all_opposite_side) return -1;
    skipT5_fixed: ;
        // Test 6 (all opposite side -> swap)
        obs_side1 = 0; if (d1_i > eps64) obs_side1 = 1; else if (d1_i < -eps64) obs_side1 = -1; else goto skipT6_fixed;
        all_opposite_side = 1;
        for (k = 0; k < n2; ++k) {
            int v = m->faces[f2].indices[k]; int xo_i = (int)lroundf(g_obsv[v].xo * 65536.0f); int yo_i = (int)lroundf(g_obsv[v].yo * 65536.0f); int zo_i = (int)lroundf(g_obsv[v].zo * 65536.0f);
            long long acc = 0;
            acc  = (((long long)a1_i * (long long)xo_i) >> 16);
            acc += (((long long)b1_i * (long long)yo_i) >> 16);
            acc += (((long long)c1_i * (long long)zo_i) >> 16);
            acc += (long long)d1_i;
            if (acc > eps64) side = 1; else if (acc < -eps64) side = -1; else continue;
            if (obs_side1 == side) { all_opposite_side = 0; break; }
        }
        if (all_opposite_side) return 1;
    skipT6_fixed: ;
        // Test 7 (all same side -> swap)
        obs_side2 = 0; if (d2_i > eps64) obs_side2 = 1; else if (d2_i < -eps64) obs_side2 = -1; else goto skipT7_fixed;
        all_same_side = 1;
        for (k = 0; k < n1; ++k) {
            int v = m->faces[f1].indices[k]; int xo_i = (int)lroundf(g_obsv[v].xo * 65536.0f); int yo_i = (int)lroundf(g_obsv[v].yo * 65536.0f); int zo_i = (int)lroundf(g_obsv[v].zo * 65536.0f);
            long long acc = 0;
            acc  = (((long long)a2_i * (long long)xo_i) >> 16);
            acc += (((long long)b2_i * (long long)yo_i) >> 16);
            acc += (((long long)c2_i * (long long)zo_i) >> 16);
            acc += (long long)d2_i;
            if (acc > eps64) side = 1; else if (acc < -eps64) side = -1; else continue;
            if (obs_side2 != side) { all_same_side = 0; break; }
        }
        if (all_same_side) return 1;
    skipT7_fixed: ;
        return 0;
    }

    // Use plane coeffs and Devant/Derriere style tests
    float a1 = m->faces[f1].plane_a, b1 = m->faces[f1].plane_b, c1 = m->faces[f1].plane_c, d1 = m->faces[f1].plane_d;
    float a2 = m->faces[f2].plane_a, b2 = m->faces[f2].plane_b, c2 = m->faces[f2].plane_c, d2 = m->faces[f2].plane_d;
    int n1 = m->faces[f1].count, n2 = m->faces[f2].count;
    // Test 6: if f2 opposite side wrt plane f1 -> f2 behind f1 => f1 should be before f2 (return 1)
    float maxabs6 = fabsf(d1); for (int kk=0; kk<n2; ++kk) { int v = m->faces[f2].indices[kk]; float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1; if (fabsf(tv) > maxabs6) maxabs6 = fabsf(tv); }
    float eps_rel6 = maxabs6 * 1e-6f; if (eps_rel6 < 0.0001f) eps_rel6 = 0.0001f;
    int pos6=0, neg6=0; for (int kk=0; kk<n2; ++kk) { int v = m->faces[f2].indices[kk]; float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1; if (tv > eps_rel6) pos6++; else if (tv < -eps_rel6) neg6++; }
    int thr6 = (3 * n2 + 3) / 4; if ((d1 > eps_rel6 && neg6 >= thr6) || (d1 < -eps_rel6 && pos6 >= thr6)) return 1;
    // Test 7: if f1 same side wrt plane f2 -> swap
    float maxabs7 = fabsf(d2); for (int kk=0; kk<n1; ++kk) { int v = m->faces[f1].indices[kk]; float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2; if (fabsf(tv) > maxabs7) maxabs7 = fabsf(tv); }
    float eps_rel7 = maxabs7 * 1e-6f; if (eps_rel7 < 0.0001f) eps_rel7 = 0.0001f;
    int pos7=0, neg7=0; for (int kk=0; kk<n1; ++kk) { int v = m->faces[f1].indices[kk]; float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2; if (tv > eps_rel7) pos7++; else if (tv < -eps_rel7) neg7++; }
    int thr7 = (3 * n1 + 3) / 4; if ((d2 > eps_rel7 && pos7 >= thr7) || (d2 < -eps_rel7 && neg7 >= thr7)) return -1;
    return 0; // inconclusive
}

/* Simple interactive inspectors (console based): ask for face id and report mismatches
   They also set highlight faces for preview via set_highlight_faces(). */
void inspect_faces_before(Model* m) {
    if (!m) { printf("No model loaded\n"); return; }
    printf("Enter target face id (or empty to cancel): "); char buf[128]; if (!fgets(buf, sizeof(buf), stdin)) return; int target = -1; if (sscanf(buf, "%d", &target) != 1) { printf("Cancelled\n"); return; }
    if (target < 0 || target >= m->face_count) { printf("Invalid face id\n"); return; }
    // Build current order by mean depth
    int* order = (int*)malloc(sizeof(int)*m->face_count); for (int i=0;i<m->face_count;i++) order[i] = i; qsort(order, m->face_count, sizeof(int), compar_face_qsort);
    // Find faces that are placed BEFORE target but should be AFTER
    int misplaced_count = 0; int *misplaced = (int*)malloc(sizeof(int)*m->face_count);
    for (int i=0;i<m->face_count;i++) {
        int f = order[i]; if (f==target) break; // only faces before target
        int res = face_should_be_before(m, f, target);
        if (res == -1) { misplaced[misplaced_count++] = f; }
    }
    printf("Found %d faces BEFORE target that should be AFTER\n", misplaced_count);
    if (misplaced_count>0) {
        for (int i=0;i<misplaced_count;i++) printf("  %d\n", misplaced[i]);
        set_highlight_faces(misplaced, misplaced_count, 6); // orange
    }
    free(misplaced); free(order);
}

void inspect_faces_after(Model* m) {
    if (!m) { printf("No model loaded\n"); return; }
    printf("Enter target face id (or empty to cancel): "); char buf[128]; if (!fgets(buf, sizeof(buf), stdin)) return; int target = -1; if (sscanf(buf, "%d", &target) != 1) { printf("Cancelled\n"); return; }
    if (target < 0 || target >= m->face_count) { printf("Invalid face id\n"); return; }
    int* order = (int*)malloc(sizeof(int)*m->face_count); for (int i=0;i<m->face_count;i++) order[i] = i; qsort(order, m->face_count, sizeof(int), compar_face_qsort);
    int misplaced_count = 0; int *misplaced = (int*)malloc(sizeof(int)*m->face_count);
    int seen_target = 0;
    for (int i=0;i<m->face_count;i++) {
        int f = order[i]; if (f==target) { seen_target = 1; continue; }
        if (!seen_target) continue;
        int res = face_should_be_before(m, f, target);
        if (res == 1) { misplaced[misplaced_count++] = f; }
    }
    printf("Found %d faces AFTER target that should be BEFORE\n", misplaced_count);
    if (misplaced_count>0) {
        for (int i=0;i<misplaced_count;i++) printf("  %d\n", misplaced[i]);
        set_highlight_faces(misplaced, misplaced_count, 12); // pink-ish
    }
    free(misplaced); free(order);
}

/* face_depth_at_point_tri removed to match ORCA/GS3Dp (no local per-triangle depth tests there) */


static void calculateFaceDepths(Model* model) {
    if (!model) return;
    for (int i = 0; i < model->face_count; ++i) {
        Face* face = &model->faces[i]; float zmin = 1e30f, zmax = -1e30f, sum = 0.0f; int n = face->count; int minx = 999999, maxx = -999999, miny = 999999, maxy = -999999; int display_flag = 1;
        for (int k = 0; k < n; ++k) {
            int vid = face->indices[k]; if (vid < 0 || vid >= model->vert_count) continue;
            ObsVertex ov = g_obsv[vid];
            float zo = ov.zo;
            /* Back-face culling optional: when enabled, faces with any vertex zo < 0 are marked non-displayable */
            if (s_cull_back_faces && zo < 0.0f) display_flag = 0;
            if (zo < zmin) zmin = zo;
            if (zo > zmax) zmax = zo;
            sum += zo;
            float px = (ov.zo == 0.0f) ? ov.xo : (ov.xo / ov.zo);
            float py = (ov.zo == 0.0f) ? ov.yo : (ov.yo / ov.zo);
            int sx = (int)lroundf((s_proj_cx - px) * -s_proj_scale + (s_proj_scale*0.5f));
            int sy = (int)lroundf((s_proj_cy - py) * -s_proj_scale + (s_proj_scale*0.5f));
            if (sx < minx) minx = sx; if (sx > maxx) maxx = sx; if (sy < miny) miny = sy; if (sy > maxy) maxy = sy;
        }
        if (!display_flag || n < 3) {
            face->plane_a = face->plane_b = face->plane_c = face->plane_d = 0.0f;
        } else {
            int idx0 = face->indices[0], idx1 = face->indices[1], idx2 = face->indices[2]; if (idx0<0||idx1<0||idx2<0) { face->plane_a = face->plane_b = face->plane_c = face->plane_d = 0.0f; }
            else {
                ObsVertex A = g_obsv[idx0], B = g_obsv[idx1], C = g_obsv[idx2];
                /* Compute plane coefficients using Fixed32 16.16 pipeline (Pascal exact) to match TestComplet */
                /* Convert observer-space floats to fixed 16.16 */
                int x1f = (int)lroundf(A.xo * 65536.0f); int y1f = (int)lroundf(A.yo * 65536.0f); int z1f = (int)lroundf(A.zo * 65536.0f);
                int x2f = (int)lroundf(B.xo * 65536.0f); int y2f = (int)lroundf(B.yo * 65536.0f); int z2f = (int)lroundf(B.zo * 65536.0f);
                int x3f = (int)lroundf(C.xo * 65536.0f); int y3f = (int)lroundf(C.yo * 65536.0f); int z3f = (int)lroundf(C.zo * 65536.0f);
                if (g_v4_pvf) {
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "x1f", "", (float)x1f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "y1f", "", (float)y1f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "z1f", "", (float)z1f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "x2f", "", (float)x2f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "y2f", "", (float)y2f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "z2f", "", (float)z2f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "x3f", "", (float)x3f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "y3f", "", (float)y3f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "z3f", "", (float)z3f);
                }
                /* Helper macros inlined: fixed_mul_64 and fixed_sub/add semantics */
                long long a_f = 0;
                long long term1 = ((long long)y1f * (long long)(z2f - z3f)) >> 16;
                long long term2 = ((long long)y2f * (long long)(z3f - z1f)) >> 16;
                long long term3 = ((long long)y3f * (long long)(z1f - z2f)) >> 16;
                a_f = term1 + term2 + term3;
                long long b_f = 0;
                b_f = -((long long)x1f * (long long)(z2f - z3f) >> 16) + ((long long)x2f * (long long)(z1f - z3f) >> 16) - ((long long)x3f * (long long)(z1f - z2f) >> 16);
                long long c_f = 0;
                c_f = ((long long)x1f * (long long)(y2f - y3f) >> 16) - ((long long)x2f * (long long)(y1f - y3f) >> 16) + ((long long)x3f * (long long)(y1f - y2f) >> 16);
                long long t1f = ((long long)y2f * (long long)z3f >> 16) - ((long long)y3f * (long long)z2f >> 16);
                long long t2f = ((long long)y1f * (long long)z3f >> 16) - ((long long)y3f * (long long)z1f >> 16);
                long long t3f = ((long long)y1f * (long long)z2f >> 16) - ((long long)y2f * (long long)z1f >> 16);
                long long d_f = -((long long)x1f * t1f >> 16) + ((long long)x2f * t2f >> 16) - ((long long)x3f * t3f >> 16);
                if (g_v4_pvf) {
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "term1", "", (float)term1);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "term2", "", (float)term2);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "term3", "", (float)term3);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "a_f", "", (float)a_f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "b_f", "", (float)b_f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "c_f", "", (float)c_f);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "d_f", "", (float)d_f);
                    /* Also emit a compact integer-only line so tools can parse fixed ints reliably */
                    if (g_v4_pvf) {
                        fprintf(g_v4_pvf, "fixed_face,%d,%d,%d,%d,%lld,%lld,%lld,%lld\n", i, x1f, y1f, z1f, a_f, b_f, c_f, d_f);
                    }
                }
                /* Save integer fixed-plane coeffs for FIXED-mode comparisons */
                face->plane_a_i = a_f; face->plane_b_i = b_f; face->plane_c_i = c_f; face->plane_d_i = d_f;
                /* Convert fixed to float (16.16) and apply V4 arrondi */
                float a = (float)((double)a_f / 65536.0);
                float b = (float)((double)b_f / 65536.0);
                float c = (float)((double)c_f / 65536.0);
                float d = (float)((double)d_f / 65536.0);
                if (g_v4_pvf) {
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "a_f_float", "", a);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "b_f_float", "", b);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "c_f_float", "", c);
                    V4_instr_log("plane_fixed", 0, i, 0, i, 0, "plane", -1, "d_f_float", "", d);
                }
                face->plane_a = V4_arrondi(a); face->plane_b = V4_arrondi(b); face->plane_c = V4_arrondi(c); face->plane_d = V4_arrondi(d);
            }
        }
        face->z_min = (n>0)?zmin:0.0f; face->z_max = (n>0)?zmax:0.0f; face->z_mean = (n>0)?(sum/n):0.0f; face->minx = (n>0)?minx:0; face->maxx = (n>0)?maxx:0; face->miny = (n>0)?miny:0; face->maxy = (n>0)?maxy:0; face->display_flag = display_flag;
    }
}

int compute_painter_order_V1(Model* m, int* order_out) {
    if (!m || !order_out) return 0;
    g_model = m;
    ObsVertex* obs = malloc(sizeof(ObsVertex)*m->vert_count); g_obsv = obs; compute_obs_vertices(m, obs);
    calculateFaceDepths(m);
    const char* dbg_env = getenv("DEBUG_DUMP_FACE_EQ");
    if (dbg_env && dbg_env[0]) {
        const char* tmp = getenv("TEMP"); char fn[1024]; if (tmp) snprintf(fn, sizeof(fn), "%s\\equ_windows.csv", tmp); else snprintf(fn, sizeof(fn), "equ_windows.csv");
        FILE* out = fopen(fn, "w"); if (out) {
            fprintf(out, "# META: angle_h=%d,angle_v=%d,angle_w=%d,distance=%.6f\n", s_angle_h, s_angle_v, s_angle_w, s_distance);
            fprintf(out, "face,a,b,c,d,z_min,z_mean,z_max,minx,maxx,miny,maxy,display_flag,vertex_indices\n");
            for (int i=0;i<m->face_count;i++) {
                Face* f = &m->faces[i]; fprintf(out, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,\"",
                    i, f->plane_a, f->plane_b, f->plane_c, f->plane_d, f->z_min, f->z_mean, f->z_max, f->minx, f->maxx, f->miny, f->maxy, f->display_flag);
                for (int k=0;k<f->count;k++) { fprintf(out, "%d ", f->indices[k]); }
                fprintf(out, "\"\n");
            }
            fclose(out);
        }
    }
    int face_count = m->face_count;
    // Build initial order: visible faces first when display_flag indicates non-visibility
    int visible_count = face_count;
    int idx = 0;
    for (int i=0;i<face_count;i++) {
        if (m->faces[i].display_flag) order_out[idx++] = i;
    }
    visible_count = idx;
    int tail = visible_count;
    for (int i=0;i<face_count;i++) {
        if (!m->faces[i].display_flag) order_out[tail++] = i;
    }
    // Sort only visible faces by z_mean
    qsort(order_out, visible_count, sizeof(int), compar_face_qsort);
    int swapped = 0;
    int ordered_pairs_capacity = face_count * 4;
    typedef struct { int face1; int face2; } OrderedPair;
    OrderedPair* ordered_pairs = NULL; int ordered_pairs_count = 0;
    if (ordered_pairs_capacity > 0) ordered_pairs = (OrderedPair*)malloc(sizeof(OrderedPair) * ordered_pairs_capacity);
    do {
        swapped = 0;
        for (int i=0;i<visible_count-1;i++) {
            int f1 = order_out[i]; int f2 = order_out[i+1];
            int already_ordered = 0;
            for (int p=0;p<ordered_pairs_count;p++) { if (ordered_pairs[p].face1==f1 && ordered_pairs[p].face2==f2) { already_ordered = 1; break; } }
            if (already_ordered) continue;
            if (m->faces[f2].z_max <= m->faces[f1].z_min) continue;
            if (m->faces[f1].z_max <= m->faces[f2].z_min) {
                order_out[i] = f2; order_out[i+1] = f1; swapped = 1;
                if (ordered_pairs != NULL && ordered_pairs_count < ordered_pairs_capacity) { ordered_pairs[ordered_pairs_count].face1 = order_out[i]; ordered_pairs[ordered_pairs_count].face2 = order_out[i+1]; ordered_pairs_count++; }
                continue;
            }
            if (m->faces[f1].maxx <= m->faces[f2].minx || m->faces[f2].maxx <= m->faces[f1].minx) continue;
            if (m->faces[f1].maxy <= m->faces[f2].miny || m->faces[f2].maxy <= m->faces[f1].miny) continue;
            // Tests ported from GS3Dp.cc (Tests 4..7). If any test decides a swap is needed, perform adjacent swap and record ordered pair.
            {
                // Prepare plane coeffs and counts
                float a1 = m->faces[f1].plane_a;
                float b1 = m->faces[f1].plane_b;
                float c1 = m->faces[f1].plane_c;
                float d1 = m->faces[f1].plane_d;

                float a2 = m->faces[f2].plane_a;
                float b2 = m->faces[f2].plane_b;
                float c2 = m->faces[f2].plane_c;
                float d2 = m->faces[f2].plane_d;

                int n1 = m->faces[f1].count;
                int n2 = m->faces[f2].count;

                // Test 4: f2 same side as observer wrt plane f1 -> ordered (no swap)
                float maxabs = fabsf(d1);
                for (int kk = 0; kk < n2; kk++) {
                    int v = m->faces[f2].indices[kk];
                    float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1;
                    if (fabsf(tv) > maxabs) {
                        maxabs = fabsf(tv);
                    }
                }

                float eps_rel_f = maxabs * 1e-6f;
                if (eps_rel_f < 0.0001f) {
                    eps_rel_f = 0.0001f;
                }

                int obs_side1 = 0;
                if (d1 > eps_rel_f) {
                    obs_side1 = 1;
                } else if (d1 < -eps_rel_f) {
                    obs_side1 = -1;
                } else {
                    goto skipT4;
                }

                int pos = 0;
                int neg = 0;
                int zero = 0;
                for (int kk = 0; kk < n2; kk++) {
                    int v = m->faces[f2].indices[kk];
                    float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1;
                    if (fabsf(tv) <= eps_rel_f) {
                        zero++;
                    } else if (tv > 0) {
                        pos++;
                    } else {
                        neg++;
                    }
                }

                int thr = (3 * n2 + 3) / 4;
                if ((obs_side1 == 1 && pos >= thr) || (obs_side1 == -1 && neg >= thr)) {
                    /* ordered correctly */
                    continue;
                }
            skipT4: ;

                // Test 5: f1 opposite side wrt plane f2 -> ordered (no swap)
                float maxabs2 = fabsf(d2);
                for (int kk = 0; kk < n1; kk++) {
                    int v = m->faces[f1].indices[kk];
                    float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2;
                    if (fabsf(tv) > maxabs2) {
                        maxabs2 = fabsf(tv);
                    }
                }

                float eps_rel2 = maxabs2 * 1e-6f;
                if (eps_rel2 < 0.0001f) {
                    eps_rel2 = 0.0001f;
                }

                int obs_side2 = 0;
                if (d2 > eps_rel2) {
                    obs_side2 = 1;
                } else if (d2 < -eps_rel2) {
                    obs_side2 = -1;
                } else {
                    goto skipT5;
                }

                int pos2 = 0;
                int neg2 = 0;
                int zero2 = 0;
                for (int kk = 0; kk < n1; kk++) {
                    int v = m->faces[f1].indices[kk];
                    float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2;
                    if (fabsf(tv) <= eps_rel2) {
                        zero2++;
                    } else if (tv > 0) {
                        pos2++;
                    } else {
                        neg2++;
                    }
                }

                int thr2 = (3 * n1 + 3) / 4;
                if ((obs_side2 == 1 && neg2 >= thr2) || (obs_side2 == -1 && pos2 >= thr2)) {
                    /* ordered correctly */
                    continue;
                }
            skipT5: ;

                // Test 6: f2 opposite side wrt plane f1 -> swap
                float maxabs6 = fabsf(d1);
                for (int kk = 0; kk < n2; kk++) {
                    int v = m->faces[f2].indices[kk];
                    float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1;
                    if (fabsf(tv) > maxabs6) {
                        maxabs6 = fabsf(tv);
                    }
                }

                float eps_rel6 = maxabs6 * 1e-6f;
                if (eps_rel6 < 0.0001f) {
                    eps_rel6 = 0.0001f;
                }

                obs_side1 = 0;
                if (d1 > eps_rel6) {
                    obs_side1 = 1;
                } else if (d1 < -eps_rel6) {
                    obs_side1 = -1;
                } else {
                    goto skipT6;
                }

                int pos6 = 0;
                int neg6 = 0;
                int zero6 = 0;
                for (int kk = 0; kk < n2; kk++) {
                    int v = m->faces[f2].indices[kk];
                    float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1;
                    if (fabsf(tv) <= eps_rel6) {
                        zero6++;
                    } else if (tv > 0) {
                        pos6++;
                    } else {
                        neg6++;
                    }
                }

                int thr6 = (3 * n2 + 3) / 4;
                if ((obs_side1 == 1 && neg6 >= thr6) || (obs_side1 == -1 && pos6 >= thr6)) {
                    goto do_swap;
                }
            skipT6: ;

                // Test 7: f1 same side wrt plane f2 -> swap
                float maxabs7 = fabsf(d2);
                for (int kk = 0; kk < n1; kk++) {
                    int v = m->faces[f1].indices[kk];
                    float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2;
                    if (fabsf(tv) > maxabs7) {
                        maxabs7 = fabsf(tv);
                    }
                }

                float eps_rel7 = maxabs7 * 1e-6f;
                if (eps_rel7 < 0.0001f) {
                    eps_rel7 = 0.0001f;
                }

                obs_side2 = 0;
                if (d2 > eps_rel7) {
                    obs_side2 = 1;
                } else if (d2 < -eps_rel7) {
                    obs_side2 = -1;
                } else {
                    goto skipT7;
                }

                int pos7 = 0;
                int neg7 = 0;
                int zero7 = 0;
                for (int kk = 0; kk < n1; kk++) {
                    int v = m->faces[f1].indices[kk];
                    float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2;
                    if (fabsf(tv) <= eps_rel7) {
                        zero7++;
                    } else if (tv > 0) {
                        pos7++;
                    } else {
                        neg7++;
                    }
                }

                int thr7 = (3 * n1 + 3) / 4;
                if ((obs_side2 == 1 && pos7 >= thr7) || (obs_side2 == -1 && neg7 >= thr7)) {
                    goto do_swap;
                } else {
                    /* inconclusive */
                    add_inconclusive_pair(f1, f2);
                    goto skipT7;
                }

            do_swap: {
                int a = order_out[i];
                int b = order_out[i+1];
                order_out[i] = b;
                order_out[i+1] = a;
                swapped = 1;
                if (ordered_pairs != NULL && ordered_pairs_count < ordered_pairs_capacity) {
                    ordered_pairs[ordered_pairs_count].face1 = order_out[i];
                    ordered_pairs[ordered_pairs_count].face2 = order_out[i+1];
                    ordered_pairs_count++;
                }
            }
            skipT7: ;
            }

        }
    } while (swapped);
    if (ordered_pairs) free(ordered_pairs);
    free(obs); g_obsv = NULL; g_model = NULL; return 1;
}

// New pairwise comparison variant: compare every face with every other face
// and perform swaps when tests indicate the opposite ordering; ordered pair
// recording and the tests are identical to the original function.
int compute_painter_order_V2(Model* m, int* order_out) {
    if (!m || !order_out) return 0;
    g_model = m;
    ObsVertex* obs = malloc(sizeof(ObsVertex)*m->vert_count); g_obsv = obs; compute_obs_vertices(m, obs);
    calculateFaceDepths(m);

    const char* dbg_env = getenv("DEBUG_DUMP_FACE_EQ");
    if (dbg_env && dbg_env[0]) {
        const char* tmp = getenv("TEMP"); char fn[1024]; if (tmp) snprintf(fn, sizeof(fn), "%s\\equ_windows_v2.csv", tmp); else snprintf(fn, sizeof(fn), "equ_windows_v2.csv");
        FILE* out = fopen(fn, "w"); if (out) {
            fprintf(out, "# META: angle_h=%d,angle_v=%d,angle_w=%d,distance=%.6f\n", s_angle_h, s_angle_v, s_angle_w, s_distance);
            fprintf(out, "face,a,b,c,d,z_min,z_mean,z_max,minx,maxx,miny,maxy,display_flag,vertex_indices\n");
            for (int i=0;i<m->face_count;i++) {
                Face* f = &m->faces[i]; fprintf(out, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,\"",
                    i, f->plane_a, f->plane_b, f->plane_c, f->plane_d, f->z_min, f->z_mean, f->z_max, f->minx, f->maxx, f->miny, f->maxy, f->display_flag);
                for (int k=0;k<f->count;k++) { fprintf(out, "%d ", f->indices[k]); }
                fprintf(out, "\"\n");
            }
            fclose(out);
        }
    }

    for (int i=0;i<m->face_count;i++) order_out[i] = i;
    qsort(order_out, m->face_count, sizeof(int), compar_face_qsort);

    int face_count = m->face_count;
    int ordered_pairs_capacity = face_count * 4;
    typedef struct { int face1; int face2; } OrderedPair;
    OrderedPair* ordered_pairs = NULL; int ordered_pairs_count = 0;
    if (ordered_pairs_capacity > 0) ordered_pairs = (OrderedPair*)malloc(sizeof(OrderedPair) * ordered_pairs_capacity);

    // Compare every pair (i,j) with i<j
    for (int i = 0; i < face_count; ++i) {
        for (int j = i+1; j < face_count; ++j) {
            int f1 = order_out[i]; int f2 = order_out[j];
            int already_ordered = 0;
            for (int p = 0; p < ordered_pairs_count; ++p) { if (ordered_pairs[p].face1==f1 && ordered_pairs[p].face2==f2) { already_ordered = 1; break; } }
            if (already_ordered) continue;

            if (m->faces[f2].z_max <= m->faces[f1].z_min) continue;
            if (m->faces[f1].z_max <= m->faces[f2].z_min) {
                int tmp = order_out[i]; order_out[i] = order_out[j]; order_out[j] = tmp;
                if (ordered_pairs != NULL && ordered_pairs_count < ordered_pairs_capacity) { ordered_pairs[ordered_pairs_count].face1 = order_out[i]; ordered_pairs[ordered_pairs_count].face2 = order_out[j]; ordered_pairs_count++; }
                continue;
            }

            if (m->faces[f1].maxx <= m->faces[f2].minx || m->faces[f2].maxx <= m->faces[f1].minx) continue;
            if (m->faces[f1].maxy <= m->faces[f2].miny || m->faces[f2].maxy <= m->faces[f1].miny) continue;

            // Same set of Tests 4..7 as original
            {
                float a1 = m->faces[f1].plane_a;
                float b1 = m->faces[f1].plane_b;
                float c1 = m->faces[f1].plane_c;
                float d1 = m->faces[f1].plane_d;

                float a2 = m->faces[f2].plane_a;
                float b2 = m->faces[f2].plane_b;
                float c2 = m->faces[f2].plane_c;
                float d2 = m->faces[f2].plane_d;

                int n1 = m->faces[f1].count;
                int n2 = m->faces[f2].count;

                // Test 4: f2 same side as observer wrt plane f1 -> ordered (no swap)
                float maxabs = fabsf(d1);
                for (int kk = 0; kk < n2; kk++) {
                    int v = m->faces[f2].indices[kk]; float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1; if (fabsf(tv) > maxabs) maxabs = fabsf(tv);
                }
                float eps_rel_f = maxabs * 1e-6f; if (eps_rel_f < 0.0001f) eps_rel_f = 0.0001f;
                int obs_side1 = 0; if (d1 > eps_rel_f) obs_side1 = 1; else if (d1 < -eps_rel_f) obs_side1 = -1; else goto skipT4_v2;
                int pos = 0, neg = 0, zero = 0;
                for (int kk = 0; kk < n2; kk++) { int v = m->faces[f2].indices[kk]; float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1; if (fabsf(tv) <= eps_rel_f) zero++; else if (tv > 0) pos++; else neg++; }
                int thr = (3 * n2 + 3) / 4;
                if ((obs_side1 == 1 && pos >= thr) || (obs_side1 == -1 && neg >= thr)) { continue; }
            skipT4_v2: ;

                // Test 5: f1 opposite side wrt plane f2 -> ordered (no swap)
                float maxabs2 = fabsf(d2);
                for (int kk = 0; kk < n1; kk++) { int v = m->faces[f1].indices[kk]; float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2; if (fabsf(tv) > maxabs2) maxabs2 = fabsf(tv); }
                float eps_rel2 = maxabs2 * 1e-6f; if (eps_rel2 < 0.0001f) eps_rel2 = 0.0001f;
                int obs_side2 = 0; if (d2 > eps_rel2) obs_side2 = 1; else if (d2 < -eps_rel2) obs_side2 = -1; else goto skipT5_v2;
                int pos2 = 0, neg2 = 0, zero2 = 0;
                for (int kk = 0; kk < n1; kk++) { int v = m->faces[f1].indices[kk]; float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2; if (fabsf(tv) <= eps_rel2) zero2++; else if (tv > 0) pos2++; else neg2++; }
                int thr2 = (3 * n1 + 3) / 4;
                if ((obs_side2 == 1 && neg2 >= thr2) || (obs_side2 == -1 && pos2 >= thr2)) { continue; }
            skipT5_v2: ;

                // Test 6: f2 opposite side wrt plane f1 -> swap
                float maxabs6 = fabsf(d1);
                for (int kk = 0; kk < n2; kk++) { int v = m->faces[f2].indices[kk]; float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1; if (fabsf(tv) > maxabs6) maxabs6 = fabsf(tv); }
                float eps_rel6 = maxabs6 * 1e-6f; if (eps_rel6 < 0.0001f) eps_rel6 = 0.0001f;
                obs_side1 = 0; if (d1 > eps_rel6) obs_side1 = 1; else if (d1 < -eps_rel6) obs_side1 = -1; else goto skipT6_v2;
                int pos6 = 0, neg6 = 0, zero6 = 0; for (int kk = 0; kk < n2; kk++) { int v = m->faces[f2].indices[kk]; float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1; if (fabsf(tv) <= eps_rel6) zero6++; else if (tv > 0) pos6++; else neg6++; }
                int thr6 = (3 * n2 + 3) / 4;
                if ((obs_side1 == 1 && neg6 >= thr6) || (obs_side1 == -1 && pos6 >= thr6)) { goto do_swap_v2; }
            skipT6_v2: ;

                // Test 7: f1 same side wrt plane f2 -> swap
                float maxabs7 = fabsf(d2);
                for (int kk = 0; kk < n1; kk++) { int v = m->faces[f1].indices[kk]; float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2; if (fabsf(tv) > maxabs7) maxabs7 = fabsf(tv); }
                float eps_rel7 = maxabs7 * 1e-6f; if (eps_rel7 < 0.0001f) eps_rel7 = 0.0001f;
                obs_side2 = 0; if (d2 > eps_rel7) obs_side2 = 1; else if (d2 < -eps_rel7) obs_side2 = -1; else goto skipT7_v2;
                int pos7 = 0, neg7 = 0, zero7 = 0; for (int kk = 0; kk < n1; kk++) { int v = m->faces[f1].indices[kk]; float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2; if (fabsf(tv) <= eps_rel7) zero7++; else if (tv > 0) pos7++; else neg7++; }
                int thr7 = (3 * n1 + 3) / 4;
                if ((obs_side2 == 1 && pos7 >= thr7) || (obs_side2 == -1 && neg7 >= thr7)) { goto do_swap_v2; } else { add_inconclusive_pair(f1,f2); goto skipT7_v2; }

            do_swap_v2: {
                int a = order_out[i]; int b = order_out[j]; order_out[i] = b; order_out[j] = a; if (ordered_pairs != NULL && ordered_pairs_count < ordered_pairs_capacity) { ordered_pairs[ordered_pairs_count].face1 = order_out[i]; ordered_pairs[ordered_pairs_count].face2 = order_out[j]; ordered_pairs_count++; }
            }
            skipT7_v2: ;
            }

        }
    }

    if (ordered_pairs) free(ordered_pairs);
    free(obs); g_obsv = NULL; g_model = NULL; return 1;
}

// V4 helpers (moved to file scope so code is C, not C++)
// Pas de 'arrondi' exact comme Delphi : Arrondi(r) => (fabs(r) < tolerance) ? 0.0 : r
static float V4_tolerance = 0.01f;
static float V4_arrondi(float r) { if (fabsf(r) < V4_tolerance) return 0.0f; return r; }
static int V4_signe(float r) { return (r >= 0.0f) ? 1 : -1; }

static void compute_plane_world(Model* m_local, int fidx, float* a, float* b, float* c, float* d) {
    int n = m_local->faces[fidx].count;
    if (n < 3) { *a = *b = *c = *d = 0.0f; return; }
    int i0 = m_local->faces[fidx].indices[0];
    int i1 = m_local->faces[fidx].indices[1];
    int i2 = m_local->faces[fidx].indices[2];
    // Use observer-space transformed coordinates (g_obsv) to match Pascal's Fixed32 pipeline
    float x1 = m_local->verts[i0].x, y1 = m_local->verts[i0].y, z1 = m_local->verts[i0].z;
    float x2 = m_local->verts[i1].x, y2 = m_local->verts[i1].y, z2 = m_local->verts[i1].z;
    float x3 = m_local->verts[i2].x, y3 = m_local->verts[i2].y, z3 = m_local->verts[i2].z;
    // Compute plane coefficients using observer-space coordinates (preferred) to match Pascal Fixed32 pipeline.
    // This avoids relying on possibly stale precomputed values and ensures arrondi is applied uniformly.
    if (g_obsv) {
        x1 = g_obsv[i0].xo; y1 = g_obsv[i0].yo; z1 = g_obsv[i0].zo;
        x2 = g_obsv[i1].xo; y2 = g_obsv[i1].yo; z2 = g_obsv[i1].zo;
        x3 = g_obsv[i2].xo; y3 = g_obsv[i2].yo; z3 = g_obsv[i2].zo;
        if (g_v4_pvf) {
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "i0", "", (float)i0);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "i1", "", (float)i1);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "i2", "", (float)i2);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "x1", "", x1);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "y1", "", y1);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "z1", "", z1);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "x2", "", x2);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "y2", "", y2);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "z2", "", z2);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "x3", "", x3);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "y3", "", y3);
            V4_instr_log("plane_verts", 0, 0, 0, fidx, 0, "plane", -1, "z3", "", z3);
        }
    }
    /* Use precomputed coefficients from calculateFaceDepths to ensure consistency
       between depth computation and the pairwise Devant/Derriere tests.  Recomputing
       from floats here caused subtle translation/rounding differences that flipped
       ordering for some pairs. */
    Face* f = &m_local->faces[fidx];
    *a = f->plane_a; *b = f->plane_b; *c = f->plane_c; *d = f->plane_d;
    if (g_v4_pvf) {
        /* log arrondi coefficients for this plane (already computed) */
        V4_instr_log("plane_arr", 0, 0, 0, fidx, 0, "plane", -1, "a_arr", "", *a);
        V4_instr_log("plane_arr", 0, 0, 0, fidx, 0, "plane", -1, "b_arr", "", *b);
        V4_instr_log("plane_arr", 0, 0, 0, fidx, 0, "plane", -1, "c_arr", "", *c);
        V4_instr_log("plane_arr", 0, 0, 0, fidx, 0, "plane", -1, "d_arr", "", *d);
    }
    // Apply Delphi Arrondi tolerance already done above
}

static int eval_side_world(Model* m_local, int fP, int fQ, float xobs, float yobs, float zobs) {
    float Ap,Bp,Cp,Dp; compute_plane_world(m_local, fP, &Ap,&Bp,&Cp,&Dp);
    int pos = 0, neg = 0;
    int nq = m_local->faces[fQ].count;
    for (int k = 0; k < nq; ++k) {
        int vid = m_local->faces[fQ].indices[k];
        /* Evaluate plane using observer-space coordinates (xo,yo,zo) to match Pascal semantics */
        float xo_v = g_obsv[vid].xo, yo_v = g_obsv[vid].yo, zo_v = g_obsv[vid].zo;
        float ax = Ap * xo_v;
        float by = Bp * yo_v;
        float cz = Cp * zo_v;
        float raw = ax + by + cz + Dp;
        float tempo = V4_arrondi(raw);
        if (g_v4_pvf) {
            /* log per-vertex component breakdown for observer-space evaluation */
            V4_instr_log("vertex_world", 0, 0, 0, fP, fQ, "vertex", vid, "ax", "", ax);
            V4_instr_log("vertex_world", 0, 0, 0, fP, fQ, "vertex", vid, "by", "", by);
            V4_instr_log("vertex_world", 0, 0, 0, fP, fQ, "vertex", vid, "cz", "", cz);
            V4_instr_log("vertex_world", 0, 0, 0, fP, fQ, "vertex", vid, "raw", "", raw);
            V4_instr_log("vertex_world", 0, 0, 0, fP, fQ, "vertex", vid, "tempo_arr", "", tempo);
        }
        if (tempo > 0.0f) pos++; else if (tempo < 0.0f) neg++;
    }
    if (pos > 0 && neg > 0) return 0; // indecisive
    int cotepoint = 0; if (pos>0) cotepoint = 1; else if (neg>0) cotepoint = -1;
    float ax = Ap * xobs;
    float by = Bp * yobs;
    float cz = Cp * zobs;
    float raw = ax + by + cz + Dp;
    float valobs = V4_arrondi(raw);
    if (g_v4_pvf) {
        V4_instr_log("valobs", 0, 0, 0, fP, fQ, "valobs", -1, "ax", "", ax);
        V4_instr_log("valobs", 0, 0, 0, fP, fQ, "valobs", -1, "by", "", by);
        V4_instr_log("valobs", 0, 0, 0, fP, fQ, "valobs", -1, "cz", "", cz);
        V4_instr_log("valobs", 0, 0, 0, fP, fQ, "valobs", -1, "raw", "", raw);
        V4_instr_log("valobs", 0, 0, 0, fP, fQ, "valobs", -1, "valobs_arr", "", valobs);
    }
    int position = V4_signe(valobs);
    if (position == cotepoint) return 1; else return -1;
}

static int limite_region_gdi(Model* m_local, int f1, int f2) {
    int n1 = m_local->faces[f1].count; int n2 = m_local->faces[f2].count;
    POINT *pt1 = (POINT*)malloc(sizeof(POINT)*n1);
    POINT *pt2 = (POINT*)malloc(sizeof(POINT)*n2);
    for (int k=0;k<n1;++k) {
        int vid = m_local->faces[f1].indices[k]; float px = g_obsv[vid].xo / (g_obsv[vid].zo==0.0f?1e-6f:g_obsv[vid].zo); float py = g_obsv[vid].yo / (g_obsv[vid].zo==0.0f?1e-6f:g_obsv[vid].zo);
        int sx = (int)lroundf((s_proj_cx - px) * -s_proj_scale + (s_proj_scale*0.5f));
        int sy = (int)lroundf((s_proj_cy - py) * -s_proj_scale + (s_proj_scale*0.5f));
        pt1[k].x = sx; pt1[k].y = sy;
    }
    for (int k=0;k<n2;++k) {
        int vid = m_local->faces[f2].indices[k]; float px = g_obsv[vid].xo / (g_obsv[vid].zo==0.0f?1e-6f:g_obsv[vid].zo); float py = g_obsv[vid].yo / (g_obsv[vid].zo==0.0f?1e-6f:g_obsv[vid].zo);
        int sx = (int)lroundf((s_proj_cx - px) * -s_proj_scale + (s_proj_scale*0.5f));
        int sy = (int)lroundf((s_proj_cy - py) * -s_proj_scale + (s_proj_scale*0.5f));
        pt2[k].x = sx; pt2[k].y = sy;
    }
    HRGN r1 = CreatePolygonRgn(pt1, n1, WINDING);
    HRGN r2 = CreatePolygonRgn(pt2, n2, WINDING);
    HRGN r = CreatePolygonRgn(pt2, n2, WINDING);
    int comb = CombineRgn(r, r1, r2, RGN_AND);
    DeleteObject(r1); DeleteObject(r2); DeleteObject(r);
    free(pt1); free(pt2);
    return (comb == NULLREGION) ? 1 : 0; // 1 = disjoint, 0 = overlap
}

/* Return 1 if two faces' projected polygons overlap (proper overlap), 0 otherwise */
int projected_polygons_overlap(Model* m, int f1, int f2) {
    if (!m) return 0;
    if (f1 < 0 || f1 >= m->face_count || f2 < 0 || f2 >= m->face_count) return 0;
    int disjoint = limite_region_gdi(m, f1, f2);
    return disjoint ? 0 : 1;
}

/* Display simple face id markers and print centers to console for debugging */
void display_model_face_ids(Model* m) {
    if (!m) return;
    ObsVertex* obs = (ObsVertex*)malloc(sizeof(ObsVertex) * m->vert_count);
    compute_obs_vertices(m, obs);
    printf("Face ID list (face:center_x,center_y in screen coords):\n");
    for (int i=0;i<m->face_count;i++) {
        Face* f = &m->faces[i]; if (f->count<=0) continue;
        float cx = 0.0f, cy = 0.0f; int n = f->count;
        for (int k=0;k<n;k++) { int vi = f->indices[k]; ObsVertex v = obs[vi]; float px = (v.zo==0.0f)?v.xo:(v.xo/v.zo); float py = (v.zo==0.0f)?v.yo:(v.yo/v.zo); cx += px; cy += py; }
        cx /= (float)n; cy /= (float)n;
        int sx, sy; screen_coords_from_proj(cx, cy, 1024, 768, s_proj_scale, s_proj_cx, s_proj_cy, &sx, &sy);
        printf("%d: %d,%d\n", i, sx, sy);
    }
    free(obs);
}



/* Version 4: implement Delphi TestComplet algorithm
   - Initialize ordering by mean depth (qsort)
   - Compare every pair i<j, applying in order:
       LimiteXY (fast bbox), LimiteRegion (polygon intersection),
       Devant(p,Q), Derriere(p,Q), Devant(Q,p), Derriere(Q,p)
     If Devant(Q,p) or Derriere(Q,p) indicate swap, perform DoPermut(0)
     else if all tests inconclusive, count as indeterminate
   - DoPermut avoids repeating the same swap; repeated blocked swaps are tracked
   - Repeat iterations until blocked-perms do not increase (matches Pascal loop)
*/
int compute_painter_order_V4(Model* m, int* order_out) {
    if (!m || !order_out) return 0;
    g_model = m;
    ObsVertex* obs = malloc(sizeof(ObsVertex)*m->vert_count); g_obsv = obs; compute_obs_vertices(m, obs);
    calculateFaceDepths(m);

    // compute observer position in world coords
    float cos_h = cosf(s_angle_h * (M_PI/180.0f)); float sin_h = sinf(s_angle_h * (M_PI/180.0f));
    float cos_v = cosf(s_angle_v * (M_PI/180.0f)); float sin_v = sinf(s_angle_v * (M_PI/180.0f));
    float xobs = s_distance * cos_h * cos_v;
    float yobs = s_distance * sin_h * cos_v;
    float zobs = s_distance * sin_v;    // apply Delphi Arrondi semantics to observer coordinates
    float xobs_a = V4_arrondi(xobs);
    float yobs_a = V4_arrondi(yobs);
    float zobs_a = V4_arrondi(zobs);
    // Initialize order via mean depth (descending)
    for (int i=0;i<m->face_count;++i) order_out[i] = i;
    qsort(order_out, m->face_count, sizeof(int), compar_face_qsort);

    int nbpermut = 0; int nbpermutb = 0; int savnbpermutb = 0; int nbiteration = 0;
    typedef struct { int a,b; } Pair; Pair *tabpermut = NULL; int tabpermut_count = 0; Pair *tabpermutb = NULL; int tabpermutb_count = 0;

    // Open a trace log to record each considered pair and action taken (help parity debugging)
    char exe_path[MAX_PATH]; GetModuleFileNameA(NULL, exe_path, MAX_PATH); char exe_dir[MAX_PATH]; strncpy(exe_dir, exe_path, MAX_PATH); char* lastbs = strrchr(exe_dir, '\\'); if (lastbs) *lastbs = '\0'; char tracefn[1024]; snprintf(tracefn, sizeof(tracefn), "%s\\pair_swaps_v4.csv", exe_dir); FILE* trf = fopen(tracefn, "w"); if (trf) { fprintf(trf, "pass,i,j,f1,f2,limiteXY,limiteRegion,dev_pq,dev_qp,action,order\n"); }
    // Open a detailed per-value trace to help find exact line-level divergences with Pascal
    char valfn[1024]; snprintf(valfn, sizeof(valfn), "%s\\pair_values_v4.csv", exe_dir); FILE* pvf = fopen(valfn, "w"); if (pvf) {
        /* extended header: include per-term components (ax,by,cz,raw) for both obs and world tempos */
        fprintf(pvf, "pass,i,j,f1,f2,stage,idx,xo,yo,zo,wx,wy,wz,ax_obs,by_obs,cz_obs,raw_obs,tempo_obs,ax_w,by_w,cz_w,raw_w,tempo_world,valobs,Ap,Bp,Cp,Dp,Aq,Bq,Cq,Dq\n");
        g_v4_pvf = pvf; /* enable per-instruction logging to same file */
        if (g_v4_pvf) {
            /* dump observer-space vertices once for parity checks */
            for (int vid=0; vid < m->vert_count; ++vid) {
                V4_instr_log("obs_dump", 0, 0, 0, -1, 0, "obs", vid, "xo", "", g_obsv[vid].xo);
                V4_instr_log("obs_dump", 0, 0, 0, -1, 0, "obs", vid, "yo", "", g_obsv[vid].yo);
                V4_instr_log("obs_dump", 0, 0, 0, -1, 0, "obs", vid, "zo", "", g_obsv[vid].zo);
            }
        }
    }

    do {
        nbiteration++;
        savnbpermutb = nbpermutb;
        nbpermutb = 0; free(tabpermutb); tabpermutb = NULL; tabpermutb_count = 0;
        for (int i = 0; i < m->face_count - 1; ++i) {
            int j = i + 1;
            int f1 = order_out[i], f2 = order_out[j];
            // LimiteXY quick bbox test
            if (m->faces[f1].maxx <= m->faces[f2].minx || m->faces[f2].maxx <= m->faces[f1].minx || m->faces[f1].maxy <= m->faces[f2].miny || m->faces[f2].maxy <= m->faces[f1].miny) continue;
            // LimiteRegion
            int limRegion = limite_region_gdi(m, f1, f2);
            if (limRegion) {
                if (trf) {
                    // write trace line (no further tests because regions disjoint)
                    fprintf(trf, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,\"", nbiteration, i, j, f1, f2, (m->faces[f1].maxx <= m->faces[f2].minx || m->faces[f2].maxx <= m->faces[f1].minx || m->faces[f1].maxy <= m->faces[f2].miny || m->faces[f2].maxy <= m->faces[f1].miny)?1:0, limRegion, 0, 0, "ok_no_overlap_region");
                    for (int kk=0; kk<m->face_count; ++kk) fprintf(trf, "%d ", order_out[kk]); fprintf(trf, "\"\n");
                }
                continue;
            }
            // Devant/Derriere tests using world planes (compute both for trace)
            // Compute plane coeffs here so we can log per-vertex values and compare to Pascal
            float Ap,Bp,Cp,Dp; compute_plane_world(m, f1, &Ap,&Bp,&Cp,&Dp);
            float Aq,Bq,Cq,Dq; compute_plane_world(m, f2, &Aq,&Bq,&Cq,&Dq);
            if (pvf) {
                // write plane coefficients
                fprintf(pvf, "%d,%d,%d,%d,%d,planes,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", nbiteration, i, j, f1, f2, Ap,Bp,Cp,Dp,Aq,Bq,Cq,Dq);
                // per-vertex tempos for f2 evaluated on plane f1
                for (int k=0; k < m->faces[f2].count; ++k) {
                    int vid = m->faces[f2].indices[k];
                    float xo_v = g_obsv[vid].xo, yo_v = g_obsv[vid].yo, zo_v = g_obsv[vid].zo;
                    float wx = m->verts[vid].x, wy = m->verts[vid].y, wz = m->verts[vid].z;
                    float ax_obs = Ap * xo_v; float by_obs = Bp * yo_v; float cz_obs = Cp * zo_v; float raw_obs = ax_obs + by_obs + cz_obs + Dp; float tempo_obs = V4_arrondi(raw_obs);
                    float ax_w = Ap * wx; float by_w = Bp * wy; float cz_w = Cp * wz; float raw_w = ax_w + by_w + cz_w + Dp; float tempo_world = V4_arrondi(raw_w);
                    fprintf(pvf, "%d,%d,%d,%d,%d,vertex_pq,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                            nbiteration, i, j, f1, f2, vid, xo_v, yo_v, zo_v, wx, wy, wz, ax_obs, by_obs, cz_obs, raw_obs, tempo_obs, ax_w, by_w, cz_w, raw_w, tempo_world);
                }
                // per-vertex tempos for f1 evaluated on plane f2
                for (int k=0; k < m->faces[f1].count; ++k) {
                    int vid = m->faces[f1].indices[k];
                    float xo_v = g_obsv[vid].xo, yo_v = g_obsv[vid].yo, zo_v = g_obsv[vid].zo;
                    float wx = m->verts[vid].x, wy = m->verts[vid].y, wz = m->verts[vid].z;
                    float ax_obs = Aq * xo_v; float by_obs = Bq * yo_v; float cz_obs = Cq * zo_v; float raw_obs = ax_obs + by_obs + cz_obs + Dq; float tempo_obs = V4_arrondi(raw_obs);
                    float ax_w = Aq * wx; float by_w = Bq * wy; float cz_w = Cq * wz; float raw_w = ax_w + by_w + cz_w + Dq; float tempo_world = V4_arrondi(raw_w);
                    fprintf(pvf, "%d,%d,%d,%d,%d,vertex_qp,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                            nbiteration, i, j, f1, f2, vid, xo_v, yo_v, zo_v, wx, wy, wz, ax_obs, by_obs, cz_obs, raw_obs, tempo_obs, ax_w, by_w, cz_w, raw_w, tempo_world);
                }
                // observer-side values
                float ax_p = Ap * xobs_a; float by_p = Bp * yobs_a; float cz_p = Cp * zobs_a; float raw_p = ax_p + by_p + cz_p + Dp; float valobs_p = V4_arrondi(raw_p);
                float ax_q = Aq * xobs_a; float by_q = Bq * yobs_a; float cz_q = Cq * zobs_a; float raw_q = ax_q + by_q + cz_q + Dq; float valobs_q = V4_arrondi(raw_q);
                fprintf(pvf, "%d,%d,%d,%d,%d,valobs,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", nbiteration, i, j, f1, f2, ax_p, by_p, cz_p, raw_p, valobs_p, ax_q, by_q, cz_q, raw_q, valobs_q);
            }
            int dev_pq = eval_side_world(m, f1, f2, xobs_a, yobs_a, zobs_a);
            int dev_qp = eval_side_world(m, f2, f1, xobs_a, yobs_a, zobs_a);
            if (trf) {
                const char* act = "undetermined";
                if (dev_pq == 1) act = "ok_devant_pq";
                else if (dev_pq == -1) act = "swap_devant_pq";
                else if (dev_qp == 1) act = "swap_devant_qp";
                // write trace before taking action
                fprintf(trf, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,\"", nbiteration, i, j, f1, f2, (m->faces[f1].maxx <= m->faces[f2].minx || m->faces[f2].maxx <= m->faces[f1].minx || m->faces[f1].maxy <= m->faces[f2].miny || m->faces[f2].maxy <= m->faces[f1].miny)?1:0, 0, dev_pq, dev_qp, act);
                for (int kk=0; kk<m->face_count; ++kk) fprintf(trf, "%d ", order_out[kk]); fprintf(trf, "\"\n");
            }
            if (dev_pq == 1) continue;
            if (dev_pq == -1) { // swap
                // DoPermut(0)
                int A = f1, B = f2; int found = 0;
                for (int pp=0; pp<tabpermut_count; ++pp) if ((tabpermut[pp].a==A && tabpermut[pp].b==B) || (tabpermut[pp].a==B && tabpermut[pp].b==A)) { found = 1; break; }
                if (found) {
                    // blocked
                    nbpermutb++; tabpermutb = (Pair*)realloc(tabpermutb, sizeof(Pair)*(tabpermutb_count+1)); tabpermutb[tabpermutb_count].a=A; tabpermutb[tabpermutb_count].b=B; tabpermutb_count++;
                    continue;
                } else {
                    tabpermut = (Pair*)realloc(tabpermut, sizeof(Pair)*(tabpermut_count+1)); tabpermut[tabpermut_count].a = (A<B)?A:B; tabpermut[tabpermut_count].b=(A<B)?B:A; tabpermut_count++;
                    // swap in order_out
                    int tmp = order_out[i]; order_out[i] = order_out[j]; order_out[j] = tmp;
                    nbpermut++;
                    continue;
                }
            }
            // symmetric test Devant(Q,p)
            if (dev_qp == 1) continue;
            if (dev_qp == -1) {
                int A = f1, B = f2; int found = 0; for (int pp=0; pp<tabpermut_count; ++pp) if ((tabpermut[pp].a==A && tabpermut[pp].b==B) || (tabpermut[pp].a==B && tabpermut[pp].b==A)) { found = 1; break; }
                if (found) { nbpermutb++; tabpermutb = (Pair*)realloc(tabpermutb, sizeof(Pair)*(tabpermutb_count+1)); tabpermutb[tabpermutb_count].a=A; tabpermutb[tabpermutb_count].b=B; tabpermutb_count++; continue; }
                tabpermut = (Pair*)realloc(tabpermut, sizeof(Pair)*(tabpermut_count+1)); tabpermut[tabpermut_count].a = (A<B)?A:B; tabpermut[tabpermut_count].b=(A<B)?B:A; tabpermut_count++;
                int tmp = order_out[i]; order_out[i] = order_out[j]; order_out[j] = tmp; nbpermut++; continue;
            }
            // if we reach here, indeterminate -> counted implicitly (Pascal increments nbindetermine)
        }
    } while ((nbpermutb > savnbpermutb) || (nbiteration <= 1 && nbpermutb > 0));

    if (tabpermut) free(tabpermut);
    if (tabpermutb) free(tabpermutb);
    if (trf) { fclose(trf); }
    if (pvf) { fclose(pvf); g_v4_pvf = NULL; }
    if (obs) free(obs); g_obsv = NULL; g_model = NULL; return 1;
}

// Dump per-pair debugging info into CSV (used to compare with Pascal TestComplet)
int dump_pairwise_debug(Model* m, int version, const char* outpath) {
    if (!m || !outpath) return 0;
    // Ensure obs vertices and face depths are computed
    ObsVertex* local_obs = malloc(sizeof(ObsVertex)*m->vert_count);
    compute_obs_vertices(m, local_obs);
    g_obsv = local_obs; g_model = m; calculateFaceDepths(m);

    FILE* f = fopen(outpath, "w"); if (!f) { free(local_obs); g_obsv = NULL; g_model = NULL; return 0; }
    fprintf(f, "# pairwise debug,version=%d,faces=%d\n", version, m->face_count);
    fprintf(f, "f1,f2,limiteXY,limiteRegion,Ap,Bp,Cp,Dp,Aq,Bq,Cq,Dq,xobs,yobs,zobs,pos_pq,neg_pq,pos_qp,neg_qp,devant_pq,devant_qp,decision\n");
    for (int i = 0; i < m->face_count - 1; ++i) {
        for (int j = i+1; j < m->face_count; ++j) {
            int f1 = i, f2 = j;
            int limiteXY = (m->faces[f1].maxx <= m->faces[f2].minx || m->faces[f2].maxx <= m->faces[f1].minx || m->faces[f1].maxy <= m->faces[f2].miny || m->faces[f2].maxy <= m->faces[f1].miny) ? 1 : 0;
            int limiteRegion = limite_region_gdi(m, f1, f2);
            // compute world observer arrondi coords used by eval
            float cos_h = cosf(s_angle_h * (M_PI/180.0f)); float sin_h = sinf(s_angle_h * (M_PI/180.0f));
            float cos_v = cosf(s_angle_v * (M_PI/180.0f)); float sin_v = sinf(s_angle_v * (M_PI/180.0f));
            float xobs = V4_arrondi(s_distance * cos_h * cos_v);
            float yobs = V4_arrondi(s_distance * sin_h * cos_v);
            float zobs = V4_arrondi(s_distance * sin_v);
            float Ap,Bp,Cp,Dp; compute_plane_world(m, f1, &Ap,&Bp,&Cp,&Dp);
            float Aq,Bq,Cq,Dq; compute_plane_world(m, f2, &Aq,&Bq,&Cq,&Dq);
            int pos_pq = 0, neg_pq = 0;
            for (int k = 0; k < m->faces[f2].count; ++k) {
                int vid = m->faces[f2].indices[k]; float tempo = V4_arrondi(Ap * g_obsv[vid].xo + Bp * g_obsv[vid].yo + Cp * g_obsv[vid].zo + Dp); if (tempo > 0.0f) pos_pq++; else if (tempo < 0.0f) neg_pq++;
            }
            int pos_qp = 0, neg_qp = 0;
            for (int k = 0; k < m->faces[f1].count; ++k) {
                int vid = m->faces[f1].indices[k]; float tempo = V4_arrondi(Aq * g_obsv[vid].xo + Bq * g_obsv[vid].yo + Cq * g_obsv[vid].zo + Dq); if (tempo > 0.0f) pos_qp++; else if (tempo < 0.0f) neg_qp++;
            }
            int dev_pq = eval_side_world(m, f1, f2, xobs, yobs, zobs);
            int dev_qp = eval_side_world(m, f2, f1, xobs, yobs, zobs);
            const char* decision = "undetermined";
            if (limiteXY) decision = "ok_no_overlap";
            else if (limiteRegion) decision = "ok_no_overlap_region";
            else if (dev_pq == 1) decision = "ok_devant_pq";
            else if (dev_qp == 1) decision = "swap_devant_qp";
            else decision = "undetermined";
            fprintf(f, "%d,%d,%d,%d", f1, f2, limiteXY, limiteRegion);
            fprintf(f, ",%.6f,%.6f,%.6f,%.6f", Ap,Bp,Cp,Dp);
            fprintf(f, ",%.6f,%.6f,%.6f,%.6f", Aq,Bq,Cq,Dq);
            fprintf(f, ",%.6f,%.6f,%.6f", xobs,yobs,zobs);
            fprintf(f, ",%d,%d,%d,%d,%d,%d,%s\n", pos_pq, neg_pq, pos_qp, neg_qp, dev_pq, dev_qp, decision);
            // write per-vertex detailed tempos as commented lines (ignored by CSV reader but useful for diff)
            for (int k = 0; k < m->faces[f2].count; ++k) {
                int vid = m->faces[f2].indices[k]; float xo_v = g_obsv[vid].xo, yo_v = g_obsv[vid].yo, zo_v = g_obsv[vid].zo; float wx = m->verts[vid].x, wy = m->verts[vid].y, wz = m->verts[vid].z; float tempo_obs = V4_arrondi(Ap * xo_v + Bp * yo_v + Cp * zo_v + Dp); float tempo_world = V4_arrondi(Ap * wx + Bp * wy + Cp * wz + Dp);
                fprintf(f, "#DETAIL,pair,%d,%d,vertex_pq,%d,xo=%.6f,yo=%.6f,zo=%.6f,wx=%.6f,wy=%.6f,wz=%.6f,tempo_obs=%.6f,tempo_world=%.6f\n", f1, f2, vid, xo_v, yo_v, zo_v, wx, wy, wz, tempo_obs, tempo_world);
            }
            for (int k = 0; k < m->faces[f1].count; ++k) {
                int vid = m->faces[f1].indices[k]; float xo_v = g_obsv[vid].xo, yo_v = g_obsv[vid].yo, zo_v = g_obsv[vid].zo; float wx = m->verts[vid].x, wy = m->verts[vid].y, wz = m->verts[vid].z; float tempo_obs = V4_arrondi(Aq * xo_v + Bq * yo_v + Cq * zo_v + Dq); float tempo_world = V4_arrondi(Aq * wx + Bq * wy + Cq * wz + Dq);
                fprintf(f, "#DETAIL,pair,%d,%d,vertex_qp,%d,xo=%.6f,yo=%.6f,zo=%.6f,wx=%.6f,wy=%.6f,wz=%.6f,tempo_obs=%.6f,tempo_world=%.6f\n", f1, f2, vid, xo_v, yo_v, zo_v, wx, wy, wz, tempo_obs, tempo_world);
            }
        }
    }
    fclose(f);
    free(local_obs); g_obsv = NULL; g_model = NULL; return 1;
}



int compute_painter_order_V3(Model* m, int* order_out) {
    if (!m || !order_out) return 0;
    g_model = m;
    ObsVertex* obs = malloc(sizeof(ObsVertex)*m->vert_count); g_obsv = obs; compute_obs_vertices(m, obs);
    calculateFaceDepths(m);
    int face_count = m->face_count;
    // Build initial order: visible faces first when display_flag indicates non-visibility
    int visible_count = face_count;
    int idx = 0;
    for (int i=0;i<face_count;i++) {
        if (m->faces[i].display_flag) order_out[idx++] = i;
    }
    visible_count = idx;
    int tail = visible_count;
    for (int i=0;i<face_count;i++) {
        if (!m->faces[i].display_flag) order_out[tail++] = i;
    }
    // Sort only visible faces by z_mean
    qsort(order_out, visible_count, sizeof(int), compar_face_qsort);
    int swapped = 0;
    int ordered_pairs_capacity = face_count * 4;
    typedef struct { int face1; int face2; } OrderedPair;
    OrderedPair* ordered_pairs = NULL; int ordered_pairs_count = 0;
    if (ordered_pairs_capacity > 0) ordered_pairs = (OrderedPair*)malloc(sizeof(OrderedPair) * ordered_pairs_capacity);
    do {
        swapped = 0;
        for (int i=0;i<visible_count-1;i++) {
            int f1 = order_out[i]; int f2 = order_out[i+1];
            int already_ordered = 0;
            for (int p=0;p<ordered_pairs_count;p++) { if (ordered_pairs[p].face1==f1 && ordered_pairs[p].face2==f2) { already_ordered = 1; break; } }
            if (already_ordered) continue;
            if (m->faces[f2].z_max <= m->faces[f1].z_min) continue;
            if (m->faces[f1].z_max <= m->faces[f2].z_min) {
                order_out[i] = f2; order_out[i+1] = f1; swapped = 1;
                if (ordered_pairs != NULL && ordered_pairs_count < ordered_pairs_capacity) { ordered_pairs[ordered_pairs_count].face1 = order_out[i]; ordered_pairs[ordered_pairs_count].face2 = order_out[i+1]; ordered_pairs_count++; }
                continue;
            }
            if (m->faces[f1].maxx <= m->faces[f2].minx || m->faces[f2].maxx <= m->faces[f1].minx) continue;
            if (m->faces[f1].maxy <= m->faces[f2].miny || m->faces[f2].maxy <= m->faces[f1].miny) continue;

            // Tests 4 and 5 only; if inconclusive -> swap and record
            {
                float a1 = m->faces[f1].plane_a;
                float b1 = m->faces[f1].plane_b;
                float c1 = m->faces[f1].plane_c;
                float d1 = m->faces[f1].plane_d;

                float a2 = m->faces[f2].plane_a;
                float b2 = m->faces[f2].plane_b;
                float c2 = m->faces[f2].plane_c;
                float d2 = m->faces[f2].plane_d;

                int n1 = m->faces[f1].count;
                int n2 = m->faces[f2].count;

                // Test 4
                float maxabs = fabsf(d1);
                for (int kk = 0; kk < n2; kk++) { int v = m->faces[f2].indices[kk]; float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1; if (fabsf(tv) > maxabs) maxabs = fabsf(tv); }
                float eps_rel_f = maxabs * 1e-6f; if (eps_rel_f < 0.0001f) eps_rel_f = 0.0001f;
                int obs_side1 = 0; if (d1 > eps_rel_f) obs_side1 = 1; else if (d1 < -eps_rel_f) obs_side1 = -1; else goto skipT4_v3;
                int pos = 0, neg = 0, zero = 0; for (int kk = 0; kk < n2; kk++) { int v = m->faces[f2].indices[kk]; float tv = a1 * g_obsv[v].xo + b1 * g_obsv[v].yo + c1 * g_obsv[v].zo + d1; if (fabsf(tv) <= eps_rel_f) zero++; else if (tv > 0) pos++; else neg++; }
                int thr = (3 * n2 + 3) / 4; if ((obs_side1 == 1 && pos >= thr) || (obs_side1 == -1 && neg >= thr)) { continue; }
            skipT4_v3: ;

                // Test 5
                float maxabs2 = fabsf(d2);
                for (int kk = 0; kk < n1; kk++) { int v = m->faces[f1].indices[kk]; float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2; if (fabsf(tv) > maxabs2) maxabs2 = fabsf(tv); }
                float eps_rel2 = maxabs2 * 1e-6f; if (eps_rel2 < 0.0001f) eps_rel2 = 0.0001f;
                int obs_side2 = 0; if (d2 > eps_rel2) obs_side2 = 1; else if (d2 < -eps_rel2) obs_side2 = -1; else goto skipT5_v3;
                int pos2 = 0, neg2 = 0, zero2 = 0; for (int kk = 0; kk < n1; kk++) { int v = m->faces[f1].indices[kk]; float tv = a2 * g_obsv[v].xo + b2 * g_obsv[v].yo + c2 * g_obsv[v].zo + d2; if (fabsf(tv) <= eps_rel2) zero2++; else if (tv > 0) pos2++; else neg2++; }
                int thr2 = (3 * n1 + 3) / 4; if ((obs_side2 == 1 && neg2 >= thr2) || (obs_side2 == -1 && pos2 >= thr2)) { continue; }
            skipT5_v3: ;

                // inconclusive -> swap and record
                {
                    int a = order_out[i]; int b = order_out[i+1]; order_out[i] = b; order_out[i+1] = a; swapped = 1;
                    if (ordered_pairs != NULL && ordered_pairs_count < ordered_pairs_capacity) { ordered_pairs[ordered_pairs_count].face1 = order_out[i]; ordered_pairs[ordered_pairs_count].face2 = order_out[i+1]; ordered_pairs_count++; }
                }
            }
        }
    } while (swapped);
    if (ordered_pairs) free(ordered_pairs);
    free(obs); g_obsv = NULL; g_model = NULL; return 1;
}

/* Wrapper dispatch: honor runtime selection set by '1'/'2' keys (default = V1) */
int compute_painter_order(Model* m, int* order_out) {
    // Reset inconclusive pairs buffer before each full order computation
    clear_inconclusive_pairs();
    if (g_painter_order_version == 1) return compute_painter_order_V1(m, order_out);
    if (g_painter_order_version == 2) return compute_painter_order_V2(m, order_out);
    if (g_painter_order_version == 3) return compute_painter_order_V3(m, order_out);
    if (g_painter_order_version == 4) return compute_painter_order_V4(m, order_out);
    return compute_painter_order_V1(m, order_out);
}

void dumpFaceEquationsCSV_Model(Model* m) {
    if (!m) return;
    const char* tmp = getenv("TEMP"); char fn[1024]; if (tmp) snprintf(fn, sizeof(fn), "%s\\equ_windows.csv", tmp); else snprintf(fn, sizeof(fn), "equ_windows.csv");
    FILE* out = fopen(fn, "w"); if (!out) return;
    fprintf(out, "# META: angle_h=%d,angle_v=%d,angle_w=%d,distance=%.6f\n", s_angle_h, s_angle_v, s_angle_w, s_distance);
    fprintf(out, "face,a,b,c,d,z_min,z_mean,z_max,minx,maxx,miny,maxy,display_flag,vertex_indices\n");
    for (int i=0;i<m->face_count;i++) {
        Face* f = &m->faces[i]; fprintf(out, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,\"",
            i, f->plane_a, f->plane_b, f->plane_c, f->plane_d, f->z_min, f->z_mean, f->z_max, f->minx, f->maxx, f->miny, f->maxy, f->display_flag);
        for (int k=0;k<f->count;k++) { fprintf(out, "%d ", f->indices[k]); }
        fprintf(out, "\"\n");
    }
    fclose(out);
}

