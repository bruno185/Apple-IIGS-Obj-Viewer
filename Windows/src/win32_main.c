#define UNICODE
#include <windows.h>
#include <commdlg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wchar.h>
#include <ctype.h>
#include "viewer.h"
#include <math.h>
#include <process.h> /* for _beginthreadex */

// Menu command IDs
#define ID_FILE_OPEN  1001
#define ID_FILE_QUIT  1002
#define ID_3D_PARAMS  2001

static Model* g_model = NULL;
static float s_ah = 30.0f;
static float s_av = 20.0f;
static float s_aw = 0.0f;
static float s_dist = 300.0f;
static float s_proj_scale = 100.0f; // angle_h, angle_v, angle_w, distance and projection scale (pixels per unit)

// Diagnostic overlays
static int g_show_inconclusive = 0; // toggled with 'I' key (on Win32)
static int s_dump_on_load = 0; // if set, dump face equations CSV after async model load
static ObsVertex* g_obs = NULL;
static int* g_order = NULL;
static int g_wireframe = 0;

// Forward declarations for functions used by WndProc
static void on_file_open(HWND hwnd);
static void on_show_params(HWND hwnd);
static void load_model_from_path(HWND hwnd, const char* path);
static void render_frame(HWND hwnd);

// Async loader hooks (WM_APP messages and helper types)
#define WM_MODEL_LOADED (WM_APP + 1)
#define WM_MODEL_LOAD_ERROR (WM_APP + 2)
#define WM_MODEL_LOAD_LOG (WM_APP + 3)

typedef struct { Model* m; int* order; } LoadResult;
typedef struct { HWND hwnd; char* path; } LoaderArg;

// Loading debug dialog globals
/* Retain dialog code but keep globals to avoid compile errors (dialog not shown by default) */
static HWND g_dbg_edit = NULL;
static HWND g_dbg_ok = NULL;
static HWND g_load_dbg = NULL;

static DWORD WINAPI loader_thread_fn(LPVOID arg);
static void start_async_load(HWND hwnd, const char* path);

static void close_load_debug(void);

static void compute_projection_and_order(HWND hwnd, int* out_winw, int* out_winh, float* out_cx, float* out_cy, float* out_scale, float* out_pxmin, float* out_pxmax, float* out_pymin, float* out_pymax) {
    RECT r;
    GetClientRect(hwnd, &r);
    int winw = r.right - r.left;
    int winh = r.bottom - r.top;
    *out_winw = winw;
    *out_winh = winh;
    compute_obs_vertices(g_model, g_obs);
    float pxmin = 1e30f;
    float pxmax = -1e30f;
    float pymin = 1e30f;
    float pymax = -1e30f;
    for (int i = 0; i < g_model->vert_count; i++) {
        ObsVertex ov = g_obs[i];
        if (ov.zo == 0.0f) continue;
        float px = ov.xo / ov.zo;
        float py = ov.yo / ov.zo;
        if (px < pxmin) pxmin = px;
        if (px > pxmax) pxmax = px;
        if (py < pymin) pymin = py;
        if (py > pymax) pymax = py;
    }

    // Match GS3Dp behavior: model-space centering applied at load.
    // Use projected-space bbox center so the projected model is centered in the view.
    float cx = 0.0f;
    float cy = 0.0f;
    float scale = s_proj_scale; // user-controlled projection scale
    // Only compute center if we found valid projected bounds
    if (pxmin <= pxmax) {
        cx = (pxmin + pxmax) * 0.5f;
        cy = (pymin + pymax) * 0.5f;
    } else {
        // fallback to origin if nothing valid
        cx = 0.0f; cy = 0.0f;
    }
    // Avoid insanely large scales that push geometry far outside GDI coordinate range
    const float MAX_PROJ_SCALE = 10000.0f;
    if (scale > MAX_PROJ_SCALE) scale = MAX_PROJ_SCALE;

    *out_cx = cx; *out_cy = cy; *out_scale = scale; *out_pxmin = pxmin; *out_pxmax = pxmax; *out_pymin = pymin; *out_pymax = pymax;
    // inform painter module of projection params (so it can compute per-face 2D bbox like GS3Dp)
    extern void set_projection_params(float cx, float cy, float scale);
    set_projection_params(cx, cy, scale);
    compute_painter_order(g_model, g_order);
}

// Apply the simple auto-fit algorithm: distance = 3 * max_dim, then compute projection scale so projected model fits the window
static void apply_auto_fit(HWND hwnd, Model* m, ObsVertex* obs) {
    if (!m || !obs) return;
    // compute bbox dimensions in model space
    float minx = 1e30f, maxx = -1e30f, miny = 1e30f, maxy = -1e30f, minz = 1e30f, maxz = -1e30f;
    for (int i = 0; i < m->vert_count; ++i) {
        float vx = m->verts[i].x, vy = m->verts[i].y, vz = m->verts[i].z;
        if (vx < minx) minx = vx; if (vx > maxx) maxx = vx;
        if (vy < miny) miny = vy; if (vy > maxy) maxy = vy;
        if (vz < minz) minz = vz; if (vz > maxz) maxz = vz;
    }
    float dx = maxx - minx;
    float dy = maxy - miny;
    float dz = maxz - minz;
    float max_dim = dx;
    if (dy > max_dim) max_dim = dy;
    if (dz > max_dim) max_dim = dz;

    // set distance to 3 * max_dim
    s_dist = 3.0f * max_dim;
    set_observer_params(s_ah, s_av, s_aw, s_dist);
    compute_obs_vertices(m, obs);

    // compute projected bounds
    float pxmin = 1e30f, pxmax = -1e30f, pymin = 1e30f, pymax = -1e30f;
    for (int i = 0; i < m->vert_count; ++i) {
        ObsVertex v = obs[i]; if (v.zo == 0.0f) continue;
        float px = v.xo / v.zo; float py = v.yo / v.zo; if (px < pxmin) pxmin = px; if (px > pxmax) pxmax = px; if (py < pymin) pymin = py; if (py > pymax) pymax = py;
    }

    RECT rc; GetClientRect(hwnd, &rc); int winw = rc.right - rc.left; int winh = rc.bottom - rc.top;
    float margin = 0.9f;
    float s1 = (winw * margin) / (pxmax - pxmin);
    float s2 = (winh * margin) / (pymax - pymin);
    // prevent division by zero / extreme scale when bbox is degenerate
    if (!(pxmax > pxmin) || !(pymax > pymin)) {
        s_proj_scale = 200.0f; // default safe scale
    } else {
        s_proj_scale = (s1 < s2) ? s1 : s2;
    }
    // clamp to avoid excessive scaling
    const float MAX_PROJ_SCALE = 10000.0f;
    if (s_proj_scale > MAX_PROJ_SCALE) s_proj_scale = MAX_PROJ_SCALE;

    extern void set_projection_params(float cx, float cy, float scale);
    // center projection on the projected bbox midpoint so apply_auto_fit actually centers the model
    float cx = 0.0f, cy = 0.0f;
    if (pxmax > pxmin) { cx = (pxmin + pxmax) * 0.5f; }
    if (pymax > pymin) { cy = (pymin + pymax) * 0.5f; }
    set_projection_params(cx, cy, s_proj_scale);
    if (g_model && g_order) compute_painter_order(g_model, g_order);
    InvalidateRect(hwnd, NULL, TRUE);
    render_frame(hwnd);
}

static void project_to_screen(float px, float py, int winw, int winh, float scale, float cx, float cy, POINT* out) {
    out->x = (int)(winw*0.5f + (px - cx) * scale);
    out->y = (int)(winh*0.5f - (py - cy) * scale);
}

// Normalize angle to range [0,360)
static float normalize_angle360(float a) {
    float r = fmodf(a, 360.0f);
    if (r < 0.0f) r += 360.0f;
    return r;
} 

static void render_frame(HWND hwnd) {
    if (!g_model) {
        // No model loaded: draw a simple informative message
        HDC hdc = GetDC(hwnd);
        RECT r; GetClientRect(hwnd, &r);
        HBRUSH bg = CreateSolidBrush(RGB(24,24,24)); FillRect(hdc, &r, bg); DeleteObject(bg);
        SetTextColor(hdc, RGB(220,220,220)); SetBkMode(hdc, TRANSPARENT);
        const char* msg = "No model loaded. Use File->Open (Ctrl+O) to load an OBJ.";
        TextOutA(hdc, 20, 20, msg, (int)strlen(msg));
        ReleaseDC(hwnd, hdc);
        return;
    }
    int winw, winh; float cx, cy, scale; float pxmin, pxmax, pymin, pymax; compute_projection_and_order(hwnd, &winw, &winh, &cx, &cy, &scale, &pxmin, &pxmax, &pymin, &pymax);
    // double buffering
    HDC hdc = GetDC(hwnd);
    HDC memdc = CreateCompatibleDC(hdc);
    HBITMAP bmp = CreateCompatibleBitmap(hdc, winw, winh);
    HBITMAP oldbmp = SelectObject(memdc, bmp);

    // background
    HBRUSH bg = CreateSolidBrush(RGB(16,16,16)); FillRect(memdc, &(RECT){0,0,winw,winh}, bg); DeleteObject(bg);



    // diagnostics: write projection bbox/scale, sample vertices and painter order to exe-dir log
    {
        char exe_path[MAX_PATH]; GetModuleFileNameA(NULL, exe_path, MAX_PATH); char exe_dir[MAX_PATH]; strncpy(exe_dir, exe_path, MAX_PATH); char* lastbs = strrchr(exe_dir, '\\'); if (lastbs) *lastbs = '\0'; char logfn[1024]; snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", exe_dir);
        FILE* lf = fopen(logfn, "a");
        if (lf) {
            fprintf(lf, "PROJ: cx=%.6f cy=%.6f pxmin=%.6f pxmax=%.6f pymin=%.6f pymax=%.6f scale=%.6f proj_scale=%.6f dist=%.6f win=%d,%d\n", cx, cy, pxmin, pxmax, pymin, pymax, scale, s_proj_scale, s_dist, winw, winh);
            int n = (g_model->vert_count<6)?g_model->vert_count:6;
            for (int i=0;i<n;i++) {
                ObsVertex v = g_obs[i]; float px = (v.zo==0.0f)?v.xo:(v.xo/v.zo); float py = (v.zo==0.0f)?v.yo:(v.yo/v.zo); fprintf(lf, "VERT[%d] xo=%.6f yo=%.6f zo=%.6f px=%.6f py=%.6f\n", i, v.xo, v.yo, v.zo, px, py);
            }
            // log painter version and full order line
            fprintf(lf, "PAINTER_ORDER: V%d, order:", g_painter_order_version);
            for (int k = 0; k < g_model->face_count; ++k) {
                fprintf(lf, " %d", g_order[k]);
            }
            fprintf(lf, "\r\n");
            fclose(lf);
        }
    }

    // draw faces
    int offscreen_count = 0;
    for (int fi=0; fi<g_model->face_count; fi++) {
        int fidx = g_order[fi]; Face* f = &g_model->faces[fidx]; POINT pts[256];
        for (int k=0;k<f->count;k++) {
            int vi = f->indices[k]; ObsVertex ov = g_obs[vi]; float px = (ov.zo==0.0f)?ov.xo:(ov.xo/ov.zo); float py = (ov.zo==0.0f)?ov.yo:(ov.yo/ov.zo);
            project_to_screen(px, py, winw, winh, scale, cx, cy, &pts[k]);
        }
        // compute face screen bbox and track offscreen faces
        int fminx = pts[0].x, fmaxx = pts[0].x, fminy = pts[0].y, fmaxy = pts[0].y;
        for (int k=1;k<f->count;k++) { if (pts[k].x < fminx) fminx = pts[k].x; if (pts[k].x > fmaxx) fmaxx = pts[k].x; if (pts[k].y < fminy) fminy = pts[k].y; if (pts[k].y > fmaxy) fmaxy = pts[k].y; }
        int visible = !((fmaxx < 0) || (fminx >= winw) || (fmaxy < 0) || (fminy >= winh));
        if (!visible) offscreen_count++;

        // All faces use the same blue color
        COLORREF col = RGB(0, 122, 255);
        if (g_wireframe) {
            // wireframe: draw thick white outline
            HPEN pen = CreatePen(PS_SOLID, 2, RGB(255,255,255)); HPEN oldp = SelectObject(memdc, pen);
            for (int k=0;k<f->count;k++) {
                int j=(k+1)%f->count; MoveToEx(memdc, pts[k].x, pts[k].y, NULL); LineTo(memdc, pts[j].x, pts[j].y);
            }
            SelectObject(memdc, oldp); DeleteObject(pen);
        } else {
            // filled: fill with face color, then draw a white outline with doubled thickness
            HBRUSH brush = CreateSolidBrush(col); HBRUSH oldb = SelectObject(memdc, brush);
            Polygon(memdc, pts, f->count);
            SelectObject(memdc, oldb); DeleteObject(brush);
            HPEN pen = CreatePen(PS_SOLID, 2, RGB(255,255,255)); HPEN oldp = SelectObject(memdc, pen);
            for (int k=0;k<f->count;k++) {
                int j=(k+1)%f->count; MoveToEx(memdc, pts[k].x, pts[k].y, NULL); LineTo(memdc, pts[j].x, pts[j].y);
            }
            SelectObject(memdc, oldp); DeleteObject(pen);
        }
    }

    // Draw inconclusive pairs overlay (if enabled)
    if (g_show_inconclusive && g_model) {
        int cnt = get_inconclusive_pair_count(); int *pairs = get_inconclusive_pairs();
        if (cnt > 0 && pairs) {
            HPEN penI = CreatePen(PS_SOLID, 3, RGB(255,255,255)); HPEN oldpi = SelectObject(memdc, penI);
            for (int pi = 0; pi < cnt; ++pi) {
                int f1 = pairs[2*pi], f2 = pairs[2*pi+1];
                if (f1 >= 0 && f1 < g_model->face_count) {
                    Face* f = &g_model->faces[f1]; POINT pts[256]; for (int k=0;k<f->count;k++) { ObsVertex ov = g_obs[f->indices[k]]; project_to_screen((ov.zo==0.0f)?ov.xo:(ov.xo/ov.zo), (ov.zo==0.0f)?ov.yo:(ov.yo/ov.zo), winw, winh, scale, cx, cy, &pts[k]); }
                    for (int k=0;k<f->count;k++) { int j=(k+1)%f->count; MoveToEx(memdc, pts[k].x, pts[k].y, NULL); LineTo(memdc, pts[j].x, pts[j].y); }
                }
                if (f2 >= 0 && f2 < g_model->face_count) {
                    Face* f = &g_model->faces[f2]; POINT pts[256]; for (int k=0;k<f->count;k++) { ObsVertex ov = g_obs[f->indices[k]]; project_to_screen((ov.zo==0.0f)?ov.xo:(ov.xo/ov.zo), (ov.zo==0.0f)?ov.yo:(ov.yo/ov.zo), winw, winh, scale, cx, cy, &pts[k]); }
                    for (int k=0;k<f->count;k++) { int j=(k+1)%f->count; MoveToEx(memdc, pts[k].x, pts[k].y, NULL); LineTo(memdc, pts[j].x, pts[j].y); }
                }
            }
            SelectObject(memdc, oldpi); DeleteObject(penI);
        }
    }

    // draw parameter overlay at bottom (Angle H, V, W, Distance, Projection scale)
    {
        char buf[256]; snprintf(buf, sizeof(buf), "H=%.1f V=%.1f W=%.1f D=%.3f S=%.1f", s_ah, s_av, s_aw, s_dist, s_proj_scale);
        // Log offscreen face info (help debug missing render)
        if (g_model) {
            char exe_path[MAX_PATH]; GetModuleFileNameA(NULL, exe_path, MAX_PATH); char exe_dir[MAX_PATH]; strncpy(exe_dir, exe_path, MAX_PATH); char* lastbs = strrchr(exe_dir, '\\'); if (lastbs) *lastbs = '\0'; char logfn[1024]; snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", exe_dir);
            FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "DRAW: offscreen_faces=%d total=%d win=%d,%d cx=%.6f cy=%.6f scale=%.6f pxmin=%.6f pxmax=%.6f pymin=%.6f pymax=%.6f\n", offscreen_count, g_model->face_count, winw, winh, cx, cy, scale, pxmin, pxmax, pymin, pymax); fclose(lf); }
        }
        SetTextColor(memdc, RGB(220,220,220)); SetBkMode(memdc, TRANSPARENT);
        TextOutA(memdc, 10, winh - 20, buf, (int)strlen(buf));
    }
    BitBlt(hdc, 0,0,winw,winh, memdc, 0,0, SRCCOPY);
    SelectObject(memdc, oldbmp); DeleteObject(bmp); DeleteDC(memdc); ReleaseDC(hwnd, hdc);
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_CREATE:
        {
            // create menu bar: File (Open, Quit) and 3D (Parameters)
            HMENU hMenu = CreateMenu();
            HMENU hFile = CreatePopupMenu(); AppendMenuW(hFile, MF_STRING, ID_FILE_OPEN, L"Open..."); AppendMenuW(hFile, MF_SEPARATOR, 0, NULL); AppendMenuW(hFile, MF_STRING, ID_FILE_QUIT, L"Quit"); AppendMenuW(hMenu, MF_POPUP, (UINT_PTR)hFile, L"File");
            HMENU h3d = CreatePopupMenu(); AppendMenuW(h3d, MF_STRING, ID_3D_PARAMS, L"Parameters..."); AppendMenuW(hMenu, MF_POPUP, (UINT_PTR)h3d, L"3D");
            SetMenu(hwnd, hMenu);
        }
        return 0;
    case WM_PAINT:
        {
            PAINTSTRUCT ps; HDC hdc = BeginPaint(hwnd, &ps); EndPaint(hwnd, &ps);
            render_frame(hwnd);
        }
        return 0;
    case WM_MODEL_LOAD_LOG:
        {
            // lParam = char* message (allocated by sender)
            char* msg = (char*)lParam; if (!msg) break;
            // append to debug dialog edit if exists
            if (g_dbg_edit) {
                int len = GetWindowTextLengthA(g_dbg_edit);
                SendMessageA(g_dbg_edit, EM_SETSEL, (WPARAM)len, (LPARAM)len);
                SendMessageA(g_dbg_edit, EM_REPLACESEL, 0, (LPARAM)msg);
            }
            // also write to exe-directory log file (viewer_win32.log)
            char exe_path[MAX_PATH]; char logfn[1024]; GetModuleFileNameA(NULL, exe_path, MAX_PATH); char exe_dir[MAX_PATH]; strncpy(exe_dir, exe_path, MAX_PATH); char* lastbs = strrchr(exe_dir, '\\'); if (lastbs) *lastbs = '\0'; snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", exe_dir); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "%s", msg); fclose(lf); }
            free(msg);
        }
        return 0;
    case WM_KEYDOWN:
        {
            int changed = 0;
            if (wParam == VK_LEFT) { s_ah -= 5.0f; changed=1; }
            else if (wParam == VK_RIGHT) { s_ah += 5.0f; changed=1; }
            else if (wParam == VK_UP) { s_av += 5.0f; changed=1; }
            else if (wParam == VK_DOWN) { s_av -= 5.0f; changed=1; }
            else if (wParam == 0x57) { g_wireframe = !g_wireframe; InvalidateRect(hwnd, NULL, TRUE); }
            else if ((wParam == 'O' && (GetKeyState(VK_CONTROL) & 0x8000)) ) { on_file_open(hwnd); }
            else if ((wParam == 'Q' && (GetKeyState(VK_CONTROL) & 0x8000)) ) { PostQuitMessage(0); }
            else if ((wParam == 'D') && (GetKeyState(VK_SHIFT) & 0x8000)) {
                if (g_model) {
                    inspect_faces_before(g_model);
                    PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Inspector: inspect_faces_before executed\r\n"));
                }
                changed = 1; InvalidateRect(hwnd, NULL, TRUE);
            }
            else if (wParam == 'D') {
                if (g_model) {
                    dumpFaceEquationsCSV_Model(g_model);
                    const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "DEBUGDUMP: wrote %s\\equ_windows.csv\n", (tmp?tmp:".") ); fclose(lf); }
                }
            }
            else if (wParam == '1') {
                g_painter_order_version = 1;
                PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Painter order: V1 selected\r\n"));
                if (g_model) {
                    int *ord = (int*)malloc(sizeof(int) * g_model->face_count);
                    if (ord) {
                        if (compute_painter_order(g_model, ord)) {
                            char exe_path[MAX_PATH]; GetModuleFileNameA(NULL, exe_path, MAX_PATH); char exe_dir[MAX_PATH]; strncpy(exe_dir, exe_path, MAX_PATH); char* lastbs = strrchr(exe_dir, '\\'); if (lastbs) *lastbs = '\0'; char ofn[1024]; snprintf(ofn, sizeof(ofn), "%s\\equ_order_runtime_v1.csv", exe_dir);
                            FILE* of = fopen(ofn, "w"); if (of) {
                                fprintf(of, "# version,1,faces,%d\\n", g_model->face_count);
                                for (int ii = 0; ii < g_model->face_count; ++ii) fprintf(of, "%d\\n", ord[ii]);
                                fclose(of);
                            }
                            char logbuf[1024]; snprintf(logbuf, sizeof(logbuf), "Painter order: wrote %s\\r\\n", ofn);
                            PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(logbuf));
                        }
                        free(ord);
                    }
                }
                InvalidateRect(hwnd, NULL, TRUE);
            }
            else if (wParam == 'B') {
                set_cull_back_faces(!get_cull_back_faces());
                char mmsg[128]; snprintf(mmsg, sizeof(mmsg), "Back-face culling: %s\r\n", get_cull_back_faces() ? "ON" : "OFF");
                PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(mmsg));
                InvalidateRect(hwnd, NULL, TRUE);
            }
            else if (wParam == 'S' && (GetKeyState(VK_SHIFT) & 0x8000)) {
                if (g_model) { inspect_faces_after(g_model); PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Inspector: inspect_faces_after executed\r\n")); }
                InvalidateRect(hwnd, NULL, TRUE);
            }
            else if (wParam == 'O' && (GetKeyState(VK_SHIFT) & 0x8000)) {
                if (g_model) { inspect_polygons_overlap(g_model, NULL, ""); PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Inspector: inspect_polygons_overlap executed\r\n")); }
                InvalidateRect(hwnd, NULL, TRUE);
            }
            else if (wParam == 'L') {
                if (g_model) { display_model_face_ids(g_model); PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Inspector: display_model_face_ids executed\r\n")); }
                InvalidateRect(hwnd, NULL, TRUE);
            }
            else if (wParam == 'I') {
                g_show_inconclusive = !g_show_inconclusive;
                char msg[128]; snprintf(msg, sizeof(msg), "Inconclusive pairs overlay: %s\r\n", g_show_inconclusive ? "ON" : "OFF");
                PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(msg));
                InvalidateRect(hwnd, NULL, TRUE);
            }
            else if (wParam == '2') {
                g_painter_order_version = 2;
                PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Painter order: V2 selected\r\n"));
                if (g_model) {
                    int *ord = (int*)malloc(sizeof(int) * g_model->face_count);
                    if (ord) {
                        if (compute_painter_order(g_model, ord)) {
                            char exe_path[MAX_PATH]; GetModuleFileNameA(NULL, exe_path, MAX_PATH); char exe_dir[MAX_PATH]; strncpy(exe_dir, exe_path, MAX_PATH); char* lastbs = strrchr(exe_dir, '\\'); if (lastbs) *lastbs = '\0'; char ofn[1024]; snprintf(ofn, sizeof(ofn), "%s\\equ_order_runtime_v2.csv", exe_dir);
                            FILE* of = fopen(ofn, "w"); if (of) {
                                fprintf(of, "# version,2,faces,%d\\n", g_model->face_count);
                                for (int ii = 0; ii < g_model->face_count; ++ii) fprintf(of, "%d\\n", ord[ii]);
                                fclose(of);
                            }
                            char logbuf[1024]; snprintf(logbuf, sizeof(logbuf), "Painter order: wrote %s\\r\\n", ofn);
                            PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(logbuf));
                        }
                        free(ord);
                    }
                }
                InvalidateRect(hwnd, NULL, TRUE);
            }
            else if (wParam == '3') {
                g_painter_order_version = 3;
                PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Painter order: V3 selected\r\n"));
                if (g_model) {
                    int *ord = (int*)malloc(sizeof(int) * g_model->face_count);
                    if (ord) {
                        if (compute_painter_order(g_model, ord)) {
                            char exe_path[MAX_PATH]; GetModuleFileNameA(NULL, exe_path, MAX_PATH); char exe_dir[MAX_PATH]; strncpy(exe_dir, exe_path, MAX_PATH); char* lastbs = strrchr(exe_dir, '\\'); if (lastbs) *lastbs = '\0'; char ofn[1024]; snprintf(ofn, sizeof(ofn), "%s\\equ_order_runtime_v3.csv", exe_dir);
                            FILE* of = fopen(ofn, "w"); if (of) {
                                fprintf(of, "# version,3,faces,%d\\n", g_model->face_count);
                                for (int ii = 0; ii < g_model->face_count; ++ii) fprintf(of, "%d\\n", ord[ii]);
                                fclose(of);
                            }
                            char logbuf[1024]; snprintf(logbuf, sizeof(logbuf), "Painter order: wrote %s\\r\\n", ofn);
                            PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(logbuf));
                        }
                        free(ord);
                    }
                }
                InvalidateRect(hwnd, NULL, TRUE);
            }
            else if (wParam == '4') {
                g_painter_order_version = 4;
                PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Painter order: V4 selected\r\n"));
                if (g_model) {
                    int *ord = (int*)malloc(sizeof(int) * g_model->face_count);
                    if (ord) {
                        if (compute_painter_order(g_model, ord)) {
                            char exe_path[MAX_PATH]; GetModuleFileNameA(NULL, exe_path, MAX_PATH); char exe_dir[MAX_PATH]; strncpy(exe_dir, exe_path, MAX_PATH); char* lastbs = strrchr(exe_dir, '\\'); if (lastbs) *lastbs = '\0'; char ofn[1024]; snprintf(ofn, sizeof(ofn), "%s\\equ_order_runtime_v4.csv", exe_dir);
                            FILE* of = fopen(ofn, "w"); if (of) {
                                fprintf(of, "# version,4,faces,%d\\n", g_model->face_count);
                                for (int ii = 0; ii < g_model->face_count; ++ii) fprintf(of, "%d\\n", ord[ii]);
                                fclose(of);
                            }
                            char logbuf[1024]; snprintf(logbuf, sizeof(logbuf), "Painter order: wrote %s\\r\\n", ofn);
                            PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(logbuf));
                            // Also produce a pairwise debug CSV to explain decisions (for parity with Pascal TestComplet)
                            char dbgfn[1024]; snprintf(dbgfn, sizeof(dbgfn), "%s\\pair_debug_v4.csv", exe_dir); if (dump_pairwise_debug(g_model, 4, dbgfn)) { char dbglog[1024]; snprintf(dbglog, sizeof(dbglog), "Painter pairwise debug: wrote %s\\r\\n", dbgfn); PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(dbglog)); } else { char dbglog[1024]; snprintf(dbglog, sizeof(dbglog), "Painter pairwise debug: failed to write %s\\r\\n", dbgfn); PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(dbglog)); }
                        }
                        free(ord);
                    }
                }
                InvalidateRect(hwnd, NULL, TRUE);
            }
            // Distance controls: A/a decreases by 10%, Z/z increases by 10%
            else if (wParam == 'A' || wParam == 'a') {
                s_dist *= 0.9f; if (s_dist < 0.01f) s_dist = 0.01f; changed = 1;
            }
            else if (wParam == 'Z' || wParam == 'z') {
                s_dist *= 1.1f; changed = 1;
            }
            // Projection scale controls: handle keypad and OEM keys here; character handling in WM_CHAR
            else if (wParam == VK_ADD || wParam == VK_OEM_PLUS) {
                s_proj_scale *= 1.1f; if (s_proj_scale > 10000.0f) s_proj_scale = 10000.0f; changed = 1;
            }
            else if (wParam == VK_SUBTRACT || wParam == VK_OEM_MINUS) {
                s_proj_scale *= 0.9f; if (s_proj_scale < 1.0f) s_proj_scale = 1.0f; changed = 1;
            }
            // Also attempt to translate this virtual key to a character for layout-independent '+'/'-' detection
            else {
                BYTE kb[256]; if (GetKeyboardState(kb)) {
                    WORD outch = 0; UINT sc = (UINT)((lParam >> 16) & 0xFF);
                    int n = ToAscii((UINT)wParam, sc, kb, &outch, 0);
                    if (n == 1) {
                        char ch = (char)outch;
                        if (ch == '+') { s_proj_scale *= 1.1f; if (s_proj_scale > 10000.0f) s_proj_scale = 10000.0f; changed = 1; }
                        else if (ch == '-') { s_proj_scale *= 0.9f; if (s_proj_scale < 1.0f) s_proj_scale = 1.0f; changed = 1; }
                    }
                }
            }
            if (changed) {
                // keep angles modulo 360
                s_ah = normalize_angle360(s_ah); s_av = normalize_angle360(s_av); s_aw = normalize_angle360(s_aw);
                // apply distance (no distance-scale)
                set_observer_params(s_ah, s_av, s_aw, s_dist);
                if (g_model && g_obs) {
                    compute_obs_vertices(g_model, g_obs);
                    if (g_order) compute_painter_order(g_model, g_order);
                }
                InvalidateRect(hwnd, NULL, TRUE);
                // force immediate render so user sees changes at once
                render_frame(hwnd);
                // log distance/scale change
                const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "UI: dist=%.6f proj_scale=%.6f\n", s_dist, s_proj_scale); fclose(lf); } 
            }
        }
        return 0;
    case WM_CHAR:
        {
            // Catch character-level '+' and '-' (layout-independent for shifted keys)
            if (wParam == '+' ) {
                s_proj_scale *= 1.1f; if (s_proj_scale > 10000.0f) s_proj_scale = 10000.0f;
            } else if (wParam == '-') {
                s_proj_scale *= 0.9f; if (s_proj_scale < 1.0f) s_proj_scale = 1.0f;
            } else if (wParam == 'A' || wParam == 'a') {
                s_dist *= 0.9f; if (s_dist < 0.01f) s_dist = 0.01f;
            } else if (wParam == 'Z' || wParam == 'z') {
                s_dist *= 1.1f;
            } else {
                break; // not handled here
            }
            // apply and render immediately
            set_observer_params(s_ah, s_av, s_aw, s_dist);
            if (g_model && g_obs) {
                compute_obs_vertices(g_model, g_obs);
                if (g_order) compute_painter_order(g_model, g_order);
            }
            InvalidateRect(hwnd, NULL, TRUE);
            render_frame(hwnd);
            const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "UI: dist=%.6f proj_scale=%.6f\n", s_dist, s_proj_scale); fclose(lf); }
        }
        return 0;
    case WM_COMMAND:
        {
            int id = LOWORD(wParam);
            if (id == ID_FILE_OPEN) { on_file_open(hwnd); }
            else if (id == ID_FILE_QUIT) { PostQuitMessage(0); }
            else if (id == ID_3D_PARAMS) { on_show_params(hwnd); }
        }
        return 0;
    case WM_MODEL_LOADED:
        {
            // lParam = LoadResult* allocated by loader thread
            LoadResult* res = (LoadResult*)lParam;
            if (!res) break;
            // Replace current model (free previous)
            if (g_model) { free_model(g_model); g_model = NULL; }
            if (g_obs) { free(g_obs); g_obs = NULL; }
            if (g_order) { free(g_order); g_order = NULL; }
            g_model = res->m; g_order = res->order; g_obs = malloc(sizeof(ObsVertex) * g_model->vert_count);
            set_observer_params(s_ah, s_av, s_aw, s_dist); compute_obs_vertices(g_model, g_obs);
            // apply brute auto-fit: distance = 3 * max_dim; projection scale computed to fit window
            apply_auto_fit(hwnd, g_model, g_obs);
            // log and redraw
            // reset inconclusive pairs between full recompute passes
            clear_inconclusive_pairs();
            const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "ASYNC LOAD: Loaded model verts=%d faces=%d\n", g_model->vert_count, g_model->face_count); fclose(lf); }
            InvalidateRect(hwnd, NULL, TRUE);
            // optionally dump face equations for debugging parity with GS3Dp
            if (s_dump_on_load) { if (g_model) {
                dumpFaceEquationsCSV_Model(g_model);
                const char* tmp2 = getenv("TEMP"); char logfn2[1024]; if (tmp2) snprintf(logfn2, sizeof(logfn2), "%s\\viewer_win32.log", tmp2); else snprintf(logfn2, sizeof(logfn2), "viewer_win32.log"); FILE* lf2 = fopen(logfn2, "a"); if (lf2) { fprintf(lf2, "DEBUGDUMP: wrote %s\\equ_windows.csv\n", (tmp2?tmp2:".")); fclose(lf2); }
                // Additionally dump the computed painter ordering for V1..V4 so we can compare algorithms
                for (int ver = 1; ver <= 4; ++ver) {
                    g_painter_order_version = ver;
                    int *ord = (int*)malloc(sizeof(int) * g_model->face_count);
                    if (ord) {
                        if (compute_painter_order(g_model, ord)) {
                            char ofn[1024]; if (tmp2) snprintf(ofn, sizeof(ofn), "%s\\equ_order_v%d.csv", tmp2, ver); else snprintf(ofn, sizeof(ofn), "equ_order_v%d.csv", ver);
                            FILE* of = fopen(ofn, "w"); if (of) {
                                fprintf(of, "# version,%d,faces,%d\n", ver, g_model->face_count);
                                for (int ii = 0; ii < g_model->face_count; ++ii) fprintf(of, "%d\n", ord[ii]);
                                fclose(of);
                            }
                        }
                        free(ord);
                    }
                }
                // restore default painter version
                g_painter_order_version = 1;
            } }
            // loading finished (no debug dialog): write a log entry
            const char* tmp3 = getenv("TEMP"); char logfn3[1024]; if (tmp3) snprintf(logfn3, sizeof(logfn3), "%s\\viewer_win32.log", tmp3); else snprintf(logfn3, sizeof(logfn3), "viewer_win32.log"); FILE* lf3 = fopen(logfn3, "a"); if (lf3) { fprintf(lf3, "Load complete.\n"); fclose(lf3); }
            free(res);
        }
        return 0;
    case WM_MODEL_LOAD_ERROR:
        {
            // lParam = char* path duplicated by loader thread
            char* bad = (char*)lParam; char msg[2048]; snprintf(msg, sizeof(msg), "Failed to load OBJ: %s", bad ? bad : "(unknown)"); MessageBoxA(hwnd, msg, "Error", MB_OK | MB_ICONERROR); if (bad) free(bad);
            // also notify debug dialog
            if (g_dbg_edit) {
                char buf[1024]; snprintf(buf, sizeof(buf), "Error: failed to load %s\r\n", (bad?bad:"(unknown)")); int len = GetWindowTextLengthA(g_dbg_edit); SendMessageA(g_dbg_edit, EM_SETSEL, (WPARAM)len, (LPARAM)len); SendMessageA(g_dbg_edit, EM_REPLACESEL, 0, (LPARAM)buf);
                if (g_dbg_ok) EnableWindow(g_dbg_ok, TRUE);
            }
        }
        return 0;
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProcW(hwnd, msg, wParam, lParam);
}

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
    int argc; LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);
    // log argv for diagnostics
    const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "argc=%d\n", argc); for (int ai=0; ai<argc; ai++) { char tmpbuf[1024]={0}; wcstombs(tmpbuf, argv[ai], sizeof(tmpbuf)); fprintf(lf, "argv[%d]=%s\n", ai, tmpbuf); } fclose(lf); }

    // Immediate startup diagnostics — log start (no blocking message box)
    {
        FILE* lf2 = fopen(logfn, "a"); if (lf2) { fprintf(lf2, "STARTUP: viewer_win32 started (no messagebox)\n"); fclose(lf2); }
    }
    // If no command-line OBJ path is provided, start GUI without loading a model.
    // User can open a model later via File->Open.

    // Reconstruct potential .obj path even if it was split over argv tokens (paths with spaces)
    char path[1024] = {0}; int last_idx = 1;
    for (int i = 1; i < argc; ++i) {
        char part[512] = {0}; wcstombs(part, argv[i], sizeof(part));
        if (i == 1) strncpy(path, part, sizeof(path)-1);
        else {
            strncat(path, " ", sizeof(path)-strlen(path)-1);
            strncat(path, part, sizeof(path)-strlen(path)-1);
        }
        last_idx = i;
        DWORD attr = GetFileAttributesA(path);
        if (attr != INVALID_FILE_ATTRIBUTES) {
            // found an existing file by progressively joining tokens
            FILE* lf2 = fopen(logfn, "a"); if (lf2) { fprintf(lf2, "Resolved path by joining argv[1..%d] -> %s\n", i, path); fclose(lf2); }
            break;
        }
        // if the current token starts with '-' then probably options follow; stop
        if (part[0] == '-') break;
    }

    // If file still doesn't exist, fall back to raw argv[1] if provided
    if (GetFileAttributesA(path) == INVALID_FILE_ATTRIBUTES) {
        if (argc > 1 && argv[1]) {
            wcstombs(path, argv[1], sizeof(path));
        } else {
            path[0] = '\0';
        }
    }
    // If still not found, search argv tokens for a token containing ".obj" (case-insensitive)
    if (GetFileAttributesA(path) == INVALID_FILE_ATTRIBUTES) {
        const char* tmp2 = getenv("TEMP"); char logfn2[1024]; if (tmp2) snprintf(logfn2, sizeof(logfn2), "%s\\viewer_win32.log", tmp2); else snprintf(logfn2, sizeof(logfn2), "viewer_win32.log");
        for (int i = 1; i < argc; ++i) {
            char orig[1024] = {0}, lower[1024] = {0};
            wcstombs(orig, argv[i], sizeof(orig));
            size_t L = strlen(orig);
            for (size_t k = 0; k < L && k + 1 < sizeof(lower); ++k) lower[k] = (char)tolower((unsigned char)orig[k]);
            lower[L] = '\0';
            if (strstr(lower, ".obj") != NULL) {
                // try this token directly
                DWORD attr = GetFileAttributesA(orig);
                if (attr != INVALID_FILE_ATTRIBUTES) {
                    strncpy(path, orig, sizeof(path)-1); path[sizeof(path)-1] = '\0';
                    FILE* lf = fopen(logfn2, "a"); if (lf) { fprintf(lf, "Found .obj token argv[%d] -> %s\n", i, orig); fclose(lf); }
                    break;
                }
                // try joining with current working dir
                char cwd[1024] = {0}; GetCurrentDirectoryA(sizeof(cwd), cwd);
                char joined[2048] = {0}; snprintf(joined, sizeof(joined), "%s\\%s", cwd, orig);
                if (GetFileAttributesA(joined) != INVALID_FILE_ATTRIBUTES) {
                    strncpy(path, joined, sizeof(path)-1); path[sizeof(path)-1] = '\0';
                    FILE* lf = fopen(logfn2, "a"); if (lf) { fprintf(lf, "Found .obj token argv[%d] -> %s (cwd joined)\n", i, joined); fclose(lf); }
                    break;
                }
            }
        }
    }

    s_ah = 30.0f; s_av = 20.0f; s_aw = 0.0f; s_dist = 300.0f; // defaults matching GS3Dp
    // parse options after the tokens consumed by the path (start at last_idx+1)
    for (int i = last_idx+1; i < argc; ++i) {
        char tmp2[512] = {0}; wcstombs(tmp2, argv[i], sizeof(tmp2));
        if (strcmp(tmp2,"--angle_h")==0 && (i+1)<argc) { char tmp3[128]; wcstombs(tmp3, argv[++i], sizeof(tmp3)); s_ah = (float)atof(tmp3); }
        else if (strcmp(tmp2,"--angle_v")==0 && (i+1)<argc) { char tmp3[128]; wcstombs(tmp3, argv[++i], sizeof(tmp3)); s_av = (float)atof(tmp3); }
        else if (strcmp(tmp2,"--angle_w")==0 && (i+1)<argc) { char tmp3[128]; wcstombs(tmp3, argv[++i], sizeof(tmp3)); s_aw = (float)atof(tmp3); }
        else if (strcmp(tmp2,"--distance")==0 && (i+1)<argc) { char tmp3[128]; wcstombs(tmp3, argv[++i], sizeof(tmp3)); s_dist = (float)atof(tmp3); }
        else if (strcmp(tmp2,"--painter")==0 && (i+1)<argc) { char tmp3[128]; wcstombs(tmp3, argv[++i], sizeof(tmp3)); int pv = atoi(tmp3); if (pv>=1 && pv<=4) { g_painter_order_version = pv; } }
        else if (strcmp(tmp2,"--dump-on-load")==0 || strcmp(tmp2,"--dump")==0) { s_dump_on_load = 1; }
    }

    // normalize angles to [0,360) before applying
    s_ah = normalize_angle360(s_ah); s_av = normalize_angle360(s_av); s_aw = normalize_angle360(s_aw);
    set_observer_params(s_ah, s_av, s_aw, s_dist);

    // If a valid .obj path was provided on the command-line, defer loading until after the window is created
    // so we can show the debug loader window and keep the UI responsive. We will start an async load later.
    if (path[0] != '\0' && GetFileAttributesA(path) != INVALID_FILE_ATTRIBUTES) {
        const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "Startup: deferring async load of path=%s\n", path); fclose(lf); }
    } else {
        // No path found; start GUI with no model
        const char* tmp3 = getenv("TEMP"); char logfn3[1024]; if (tmp3) snprintf(logfn3, sizeof(logfn3), "%s\\viewer_win32.log", tmp3); else snprintf(logfn3, sizeof(logfn3), "viewer_win32.log"); FILE* lf4 = fopen(logfn3, "a"); if (lf4) { fprintf(lf4, "No OBJ path provided; starting GUI without model\n"); fclose(lf4); }
    }

    // register window class
    const wchar_t CLASS_NAME[] = L"GS3DpWin";
    WNDCLASSW wc = {0}; wc.lpfnWndProc = WndProc; wc.hInstance = hInstance; wc.lpszClassName = CLASS_NAME; wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    RegisterClassW(&wc);

    HWND hwnd = CreateWindowExW(0, CLASS_NAME, L"GS3Dp Viewer (Win32 GDI)", WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 1024, 768, NULL, NULL, hInstance, NULL);
    if (!hwnd) { MessageBoxW(NULL, L"CreateWindow failed", L"Error", MB_OK); return 1; }

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    // Ensure main window is visible and focused even when no model was provided
    if (path[0] == '\0' || GetFileAttributesA(path) == INVALID_FILE_ATTRIBUTES) {
        SetWindowTextW(hwnd, L"GS3Dp Viewer (Win32 GDI) - No model loaded");
        ShowWindow(hwnd, SW_SHOWNORMAL);
        BringWindowToTop(hwnd);
        SetForegroundWindow(hwnd);
        // Non-intrusive behavior: do not show a blocking message box when no OBJ was supplied
        const char* tmp3 = getenv("TEMP"); char logfn3[1024]; if (tmp3) snprintf(logfn3, sizeof(logfn3), "%s\\viewer_win32.log", tmp3); else snprintf(logfn3, sizeof(logfn3), "viewer_win32.log"); FILE* lf3 = fopen(logfn3, "a"); if (lf3) { fprintf(lf3, "UI: started with no OBJ path; main window shown (no messagebox)\n"); fclose(lf3); }
    }

    // Additional diagnostics: log visibility and cmdshow for debugging double-click exits
    {
        const char* tmp2 = getenv("TEMP"); char logfn2[1024]; if (tmp2) snprintf(logfn2, sizeof(logfn2), "%s\\viewer_win32.log", tmp2); else snprintf(logfn2, sizeof(logfn2), "viewer_win32.log"); FILE* lf = fopen(logfn2, "a"); if (lf) { fprintf(lf, "UI: ShowWindow called nCmdShow=%d IsWindowVisible=%d\n", nCmdShow, IsWindowVisible(hwnd)); fclose(lf); }
    }

    // initial paint
    InvalidateRect(hwnd, NULL, TRUE);
    // also force an initial render to ensure we draw once
    render_frame(hwnd);

    // If a path was supplied on the command-line, start async load now that the window exists
    if (path[0] != '\0' && GetFileAttributesA(path) != INVALID_FILE_ATTRIBUTES) {
        start_async_load(hwnd, path);
    }
    {
        const char* tmp3 = getenv("TEMP"); char logfn3[1024]; if (tmp3) snprintf(logfn3, sizeof(logfn3), "%s\\viewer_win32.log", tmp3); else snprintf(logfn3, sizeof(logfn3), "viewer_win32.log"); FILE* lf4 = fopen(logfn3, "a"); if (lf4) { fprintf(lf4, "Called render_frame() after window create\n"); fclose(lf4); }
    }

    // message loop
    MSG msg; while (GetMessageW(&msg, NULL, 0,0)) { TranslateMessage(&msg); DispatchMessageW(&msg); }

    free(g_obs); free(g_order); free_model(g_model);
    return 0;
}

// --------- New helpers: File Open and Parameters dialog ----------

static void load_model_from_path(HWND hwnd, const char* path) {
    if (!path) return;
    // Try to load
    Model* nm = load_obj(path);
    if (!nm) { char msg[1024]; snprintf(msg, sizeof(msg), "Failed to load: %s", path); MessageBoxA(hwnd, msg, "Error", MB_OK | MB_ICONERROR); return; }
    // free previous
    if (g_model) { free_model(g_model); g_model = NULL; }
    g_model = nm;
    // allocate observers and order buffers
    if (g_obs) { free(g_obs); g_obs = NULL; }
    if (g_order) { free(g_order); g_order = NULL; }
    g_obs = malloc(sizeof(ObsVertex) * g_model->vert_count);
    g_order = malloc(sizeof(int) * g_model->face_count);
    // recompute
    set_observer_params(s_ah, s_av, s_aw, s_dist);
    compute_obs_vertices(g_model, g_obs);
    compute_painter_order(g_model, g_order);
    // apply brute auto-fit after load
    apply_auto_fit(hwnd, g_model, g_obs);
    InvalidateRect(hwnd, NULL, TRUE);
}

// loader implementation (types & messages declared above)
static DWORD WINAPI loader_thread_fn(LPVOID arg) {
    LoaderArg* a = (LoaderArg*)arg; if (!a) return 0; HWND hwnd = a->hwnd; char* path = a->path;
    // log start
    char tmpbuf[1024]; snprintf(tmpbuf, sizeof(tmpbuf), "Starting background load: %s\r\n", path);
    PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(tmpbuf));

    // also write a direct loader thread debug file (helps when PostMessage doesn't show later)
    const char* tmp = getenv("TEMP"); char loader_logfn[1024]; if (tmp) snprintf(loader_logfn, sizeof(loader_logfn), "%s\\loader_thread_debug.log", tmp); else snprintf(loader_logfn, sizeof(loader_logfn), "loader_thread_debug.log"); FILE* lsf = fopen(loader_logfn, "a"); if (lsf) { fprintf(lsf, "THREAD: starting load for %s\n", path); fclose(lsf); }

    // defensive check: ensure file can be fopen'ed inside the thread
    FILE* tf = fopen(path, "r");
    if (!tf) {
        char errbuf[1024]; snprintf(errbuf, sizeof(errbuf), "loader_thread: fopen failed on %s (errno=%d)\r\n", path, errno);
        PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(errbuf));
        if (lsf) { FILE* lsf2 = fopen(loader_logfn, "a"); if (lsf2) { fprintf(lsf2, "THREAD: fopen FAILED errno=%d\n", errno); fclose(lsf2); } }
        char* pcopy = _strdup(path);
        PostMessageW(hwnd, WM_MODEL_LOAD_ERROR, 0, (LPARAM)pcopy);
        free(a);
        return 0;
    } else {
        fclose(tf);
        PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("loader_thread: fopen succeeded\r\n"));
        if (lsf) { FILE* lsf2 = fopen(loader_logfn, "a"); if (lsf2) { fprintf(lsf2, "THREAD: fopen OK\n"); fclose(lsf2); } }
    }

    // attempt to load the model (runs in background)
    if (lsf) { FILE* lsf3 = fopen(loader_logfn, "a"); if (lsf3) { fprintf(lsf3, "THREAD: calling load_obj()\n"); fclose(lsf3); } }
    Model* nm = load_obj(path);
    if (!nm) {
        // post a message back to the main thread with the failed path so UI can show error
        char errbuf[1024]; snprintf(errbuf, sizeof(errbuf), "Failed to load OBJ: %s (load_obj returned NULL)\r\n", path);
        PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(errbuf));
        if (lsf) { FILE* lsf4 = fopen(loader_logfn, "a"); if (lsf4) { fprintf(lsf4, "THREAD: load_obj returned NULL\n"); fclose(lsf4); } }
        char* pcopy = _strdup(path);
        PostMessageW(hwnd, WM_MODEL_LOAD_ERROR, 0, (LPARAM)pcopy);
        free(a);
        return 0;
    }
    PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("OBJ loaded, parsing complete\r\n"));
    if (lsf) { FILE* lsf5 = fopen(loader_logfn, "a"); if (lsf5) { fprintf(lsf5, "THREAD: load_obj returned model verts=%d faces=%d\n", nm->vert_count, nm->face_count); fclose(lsf5); } }

    // compute painter order in background
    int* order_buf = malloc(sizeof(int) * nm->face_count);
    if (!order_buf) { free_model(nm); char* pcopy = _strdup(path); PostMessageW(hwnd, WM_MODEL_LOAD_ERROR, 0, (LPARAM)pcopy); free(a); return 0; }
    PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Computing painter order...\r\n"));
    if (lsf) { FILE* lsf6 = fopen(loader_logfn, "a"); if (lsf6) { fprintf(lsf6, "THREAD: calling compute_painter_order\n"); fclose(lsf6); } }
    compute_painter_order(nm, order_buf);
    PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Painter order computed.\r\n"));
    if (lsf) { FILE* lsf7 = fopen(loader_logfn, "a"); if (lsf7) { fprintf(lsf7, "THREAD: compute_painter_order done\n"); fclose(lsf7); } }
    // Also attempt to write pairwise debug for V4 on load to help diagnose parity issues (writes to exe dir)
    char exe_path[MAX_PATH]; GetModuleFileNameA(NULL, exe_path, MAX_PATH); char exe_dir[MAX_PATH]; strncpy(exe_dir, exe_path, MAX_PATH); char* lastbs = strrchr(exe_dir, '\\'); if (lastbs) *lastbs = '\0'; char dbgfn[1024]; snprintf(dbgfn, sizeof(dbgfn), "%s\\pair_debug_v4_onload.csv", exe_dir); if (dump_pairwise_debug(nm, 4, dbgfn)) { char dbgl[1024]; snprintf(dbgl, sizeof(dbgl), "Painter pairwise debug: wrote %s\\r\\n", dbgfn); PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(dbgl)); } else { char dbgl[1024]; snprintf(dbgl, sizeof(dbgl), "Painter pairwise debug: failed to write %s\\r\\n", dbgfn); PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(dbgl)); }

    // post result
    LoadResult* res = (LoadResult*)malloc(sizeof(LoadResult)); res->m = nm; res->order = order_buf;
    PostMessageW(hwnd, WM_MODEL_LOADED, 0, (LPARAM)res);
    // signal the UI that loading is finished (no dialog shown); write a log entry instead
    PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Load complete.\r\n"));
    if (lsf) { FILE* lsf8 = fopen(loader_logfn, "a"); if (lsf8) { fprintf(lsf8, "THREAD: posted WM_MODEL_LOADED and finished\n"); fclose(lsf8); } }
    free(a);
    return 0;
}

static void start_async_load(HWND hwnd, const char* path) {
    if (!path) return; 
    // log async load start (no UI dialogs)
    const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf0 = fopen(logfn, "a"); if (lf0) { fprintf(lf0, "UI: start_async_load called for %s\n", path); fclose(lf0); }
    LoaderArg* a = (LoaderArg*)malloc(sizeof(LoaderArg)); if (!a) return;
    a->hwnd = hwnd; a->path = _strdup(path);
    // Use _beginthreadex to ensure CRT is properly initialized in the new thread
    unsigned tidu = 0; uintptr_t th = _beginthreadex(NULL, 0, (unsigned (__stdcall *)(void*))loader_thread_fn, a, 0, &tidu);
    if (!th) { FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "UI: _beginthreadex failed for %s\n", path); fclose(lf); } free(a->path); free(a); }
    else { FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "UI: _beginthreadex started tid=%u for %s\n", (unsigned)tidu, path); fclose(lf); } CloseHandle((HANDLE)th); }
}

static LRESULT CALLBACK LoadDlgProc(HWND dlg, UINT msg, WPARAM wParam, LPARAM lParam) {
    if (msg == WM_CREATE) {
        // create a multiline read-only edit to show logs, and an OK button disabled until done
        g_dbg_edit = CreateWindowExW(0, L"EDIT", L"", WS_CHILD|WS_VISIBLE|WS_BORDER|ES_LEFT|ES_READONLY|ES_MULTILINE|WS_VSCROLL|ES_AUTOVSCROLL, 10,10,360,140,dlg,(HMENU)401,GetModuleHandle(NULL),NULL);
        g_dbg_ok = CreateWindowExW(0, L"BUTTON", L"OK", WS_CHILD|WS_VISIBLE|BS_DEFPUSHBUTTON|WS_DISABLED, 150,160,80,26,dlg,(HMENU)301,GetModuleHandle(NULL),NULL);
        return 0;
    }
    if (msg == WM_COMMAND) {
        int id = LOWORD(wParam);
        if (id == 301) { DestroyWindow(dlg); return 0; }
    }
    if (msg == WM_DESTROY) {
        g_load_dbg = NULL; g_dbg_edit = NULL; g_dbg_ok = NULL; return 0; }
    return DefWindowProcW(dlg, msg, wParam, lParam);
}

static void create_load_debug(HWND parent, const char* path) {
    // log invocation
    const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "UI: create_load_debug called for path=%s\n", path); fclose(lf); }

    if (g_load_dbg) return; // already shown
    const wchar_t CLS[] = L"GS3D_LoadDlg";
    WNDCLASSW wc = {0}; wc.lpfnWndProc = LoadDlgProc; wc.hInstance = GetModuleHandle(NULL); wc.lpszClassName = CLS; wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    RegisterClassW(&wc);
    HWND dlg = CreateWindowExW(0, CLS, L"Loading model...", WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU, CW_USEDEFAULT, CW_USEDEFAULT, 390, 230, parent, NULL, GetModuleHandle(NULL), NULL);
    if (!dlg) { FILE* lf2 = fopen(logfn, "a"); if (lf2) { fprintf(lf2, "ERROR: create_load_debug CreateWindowExW failed\n"); fclose(lf2); } return; }
    // center over parent
    RECT rc; GetWindowRect(parent, &rc); int cx = (rc.left+rc.right)/2 - 195; int cy = (rc.top+rc.bottom)/2 - 115; SetWindowPos(dlg, HWND_TOP, cx, cy, 390, 230, SWP_SHOWWINDOW);
    ShowWindow(dlg, SW_SHOW);
    g_load_dbg = dlg;
    // set initial message with path (post to main window so handler appends to edit)
    char buf[1024]; snprintf(buf, sizeof(buf), "Loading: %s\r\n", path); PostMessageW(parent, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup(buf));
}

static void close_load_debug(void) {
    if (g_load_dbg) { DestroyWindow(g_load_dbg); g_load_dbg = NULL; g_dbg_edit = NULL; g_dbg_ok = NULL; }
}

static void on_file_open(HWND hwnd) {
    WCHAR file[MAX_PATH] = {0};
    OPENFILENAMEW ofn; ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn); ofn.hwndOwner = hwnd; ofn.lpstrFile = file; ofn.nMaxFile = MAX_PATH;
    ofn.lpstrFilter = L"OBJ Files\0*.obj\0All Files\0*.*\0"; ofn.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST;
    ofn.lpstrTitle = L"Open OBJ file";
    if (GetOpenFileNameW(&ofn)) {
        char pathA[MAX_PATH]; wcstombs(pathA, file, sizeof(pathA));
        // log selection
        const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "UI: on_file_open selected=%s\n", pathA); fclose(lf); }
        // start background load so UI stays responsive
        start_async_load(hwnd, pathA);
    }
}

// Parameters dialog implemented as a simple owned window with controls (Angle H, Angle V, Angle W, Distance)
static LRESULT CALLBACK ParamsWndProc(HWND dlg, UINT msg, WPARAM wParam, LPARAM lParam) {
    static HWND hEditAH, hEditAV, hEditAW, hEditDist;
    if (msg == WM_CREATE) {
        wchar_t buf[64];
        swprintf(buf, 64, L"%.0f", s_ah); hEditAH = CreateWindowExW(0, L"EDIT", buf, WS_CHILD|WS_VISIBLE|WS_BORDER|ES_LEFT, 140,20,180,24,dlg,(HMENU)101,GetModuleHandle(NULL),NULL);
        swprintf(buf, 64, L"%.0f", s_av); hEditAV = CreateWindowExW(0, L"EDIT", buf, WS_CHILD|WS_VISIBLE|WS_BORDER|ES_LEFT, 140,60,180,24,dlg,(HMENU)102,GetModuleHandle(NULL),NULL);
        swprintf(buf, 64, L"%.0f", s_aw); hEditAW = CreateWindowExW(0, L"EDIT", buf, WS_CHILD|WS_VISIBLE|WS_BORDER|ES_LEFT, 140,100,180,24,dlg,(HMENU)103,GetModuleHandle(NULL),NULL);
        swprintf(buf, 64, L"%.3f", s_dist); hEditDist = CreateWindowExW(0, L"EDIT", buf, WS_CHILD|WS_VISIBLE|WS_BORDER|ES_LEFT, 140,140,180,24,dlg,(HMENU)104,GetModuleHandle(NULL),NULL);
        swprintf(buf, 64, L"%.1f", s_proj_scale); HWND hEditProjScale = CreateWindowExW(0, L"EDIT", buf, WS_CHILD|WS_VISIBLE|WS_BORDER|ES_LEFT, 140,180,180,24,dlg,(HMENU)106,GetModuleHandle(NULL),NULL);
        // painter version combo
        HWND hComboPainter = CreateWindowExW(0, L"COMBOBOX", NULL, WS_CHILD|WS_VISIBLE|CBS_DROPDOWNLIST|WS_VSCROLL, 140,210,180,120,dlg,(HMENU)107,GetModuleHandle(NULL),NULL);
        SendMessageW(hComboPainter, CB_ADDSTRING, 0, (LPARAM)L"V1 - adjacent swap (original)");
        SendMessageW(hComboPainter, CB_ADDSTRING, 0, (LPARAM)L"V2 - pairwise compare");
        SendMessageW(hComboPainter, CB_ADDSTRING, 0, (LPARAM)L"V3 - V1 without Tests 6/7");
        SendMessageW(hComboPainter, CB_ADDSTRING, 0, (LPARAM)L"V4 - Delphi TestComplet algorithm");
        // future versions can be added here
        int sel = (g_painter_order_version == 2) ? 1 : (g_painter_order_version == 3) ? 2 : (g_painter_order_version == 4) ? 3 : 0;
        SendMessageW(hComboPainter, CB_SETCURSEL, (WPARAM)sel, 0);

        // painter mode combo (GS3Dp FAST / FIXED / FLOAT equivalence)
        HWND hComboPainterMode = CreateWindowExW(0, L"COMBOBOX", NULL, WS_CHILD|WS_VISIBLE|CBS_DROPDOWNLIST|WS_VSCROLL, 140,240,180,120,dlg,(HMENU)108,GetModuleHandle(NULL),NULL);
        SendMessageW(hComboPainterMode, CB_ADDSTRING, 0, (LPARAM)L"FAST - simple sort (fast)");
        SendMessageW(hComboPainterMode, CB_ADDSTRING, 0, (LPARAM)L"FIXED - Fixed32 pipeline");
        SendMessageW(hComboPainterMode, CB_ADDSTRING, 0, (LPARAM)L"FLOAT - float painter");
        int pmode = (g_painter_mode == PAINTER_MODE_FIXED) ? 1 : (g_painter_mode == PAINTER_MODE_FLOAT) ? 2 : 0;
        SendMessageW(hComboPainterMode, CB_SETCURSEL, (WPARAM)pmode, 0);

        // labels
        CreateWindowExW(0, L"STATIC", L"Angle H:", WS_CHILD|WS_VISIBLE, 20,20,110,24,dlg,NULL,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"STATIC", L"Angle V:", WS_CHILD|WS_VISIBLE, 20,60,110,24,dlg,NULL,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"STATIC", L"Angle W:", WS_CHILD|WS_VISIBLE, 20,100,110,24,dlg,NULL,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"STATIC", L"Distance:", WS_CHILD|WS_VISIBLE, 20,140,110,24,dlg,NULL,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"STATIC", L"Projection scale:", WS_CHILD|WS_VISIBLE, 20,180,110,24,dlg,NULL,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"STATIC", L"Painter order:", WS_CHILD|WS_VISIBLE, 20,210,110,24,dlg,NULL,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"STATIC", L"Painter mode:", WS_CHILD|WS_VISIBLE, 20,240,110,24,dlg,NULL,GetModuleHandle(NULL),NULL);
        // buttons (moved down to fit extra combo)
        CreateWindowExW(0, L"BUTTON", L"OK", WS_CHILD|WS_VISIBLE|BS_DEFPUSHBUTTON, 60,300,90,28,dlg,(HMENU)201,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"BUTTON", L"Cancel", WS_CHILD|WS_VISIBLE, 200,300,90,28,dlg,(HMENU)202,GetModuleHandle(NULL),NULL);
        return 0;
    }
    if (msg == WM_COMMAND) {
        int id = LOWORD(wParam);
        if (id == 201) {
            wchar_t buf[64]; GetWindowTextW(hEditAH, buf, 64); s_ah = (float)wcstod(buf, NULL);
            GetWindowTextW(hEditAV, buf, 64); s_av = (float)wcstod(buf, NULL);
            GetWindowTextW(hEditAW, buf, 64); s_aw = (float)wcstod(buf, NULL);
            GetWindowTextW(hEditDist, buf, 64); s_dist = (float)wcstod(buf, NULL);
            GetWindowTextW((HWND)GetDlgItem(dlg, 106), buf, 64); s_proj_scale = (float)wcstod(buf, NULL);
            if (s_proj_scale < 1.0f) s_proj_scale = 1.0f;
            // painter version selection from combo (ID 107)
            int sel = (int)SendMessageW((HWND)GetDlgItem(dlg, 107), CB_GETCURSEL, 0, 0);
            if (sel >= 0) {
                g_painter_order_version = sel + 1; // combo index 0 -> V1, 1 -> V2
                char* msg = _strdup("Painter order selection changed via Parameters dialog\r\n");
                HWND ownerWnd = GetWindow(dlg, GW_OWNER); if (!ownerWnd) ownerWnd = GetParent(dlg); if (!ownerWnd) ownerWnd = GetAncestor(dlg, GA_ROOTOWNER);
                if (ownerWnd) PostMessageW(ownerWnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)msg); else free(msg);
            }
            // painter mode selection from combo (ID 108)
            int psel = (int)SendMessageW((HWND)GetDlgItem(dlg, 108), CB_GETCURSEL, 0, 0);
            if (psel >= 0) {
                // Map combo index 0->FAST, 1->FIXED, 2->FLOAT
                set_painter_mode(psel);
                char* msg2 = _strdup("Painter mode selection changed via Parameters dialog\r\n");
                HWND ownerWnd2 = GetWindow(dlg, GW_OWNER); if (!ownerWnd2) ownerWnd2 = GetParent(dlg); if (!ownerWnd2) ownerWnd2 = GetAncestor(dlg, GA_ROOTOWNER);
                if (ownerWnd2) PostMessageW(ownerWnd2, WM_MODEL_LOAD_LOG, 0, (LPARAM)msg2); else free(msg2);
            }
            // normalize angles to [0,360)
            s_ah = normalize_angle360(s_ah); s_av = normalize_angle360(s_av); s_aw = normalize_angle360(s_aw);
            // apply params and recompute using GS3Dp semantics (distance only)
            set_observer_params(s_ah, s_av, s_aw, s_dist);
            // apply updated projection scale to painter
            extern void set_projection_params(float cx, float cy, float scale);
            set_projection_params(0.0f, 0.0f, s_proj_scale);
            if (g_model && g_obs) {
                compute_obs_vertices(g_model, g_obs);
                if (g_order) compute_painter_order(g_model, g_order);
            }
            HWND owner = GetWindow(dlg, GW_OWNER);
            if (!owner) owner = GetParent(dlg);
            if (!owner) owner = GetAncestor(dlg, GA_ROOTOWNER);
            if (owner) {
                InvalidateRect(owner, NULL, TRUE);
                UpdateWindow(owner);
                render_frame(owner);
            } else {
                // Fallback: invalidate main window class by forcing a paint
                InvalidateRect(NULL, NULL, TRUE);
                // render_frame may be invoked on the current active window if available
                HWND active = GetActiveWindow(); if (active) render_frame(active);
            }
            DestroyWindow(dlg);
            return 0;
        } else if (id == 202) { DestroyWindow(dlg); return 0; }
    }
    return DefWindowProcW(dlg, msg, wParam, lParam);
}

static void on_show_params(HWND hwnd) {
    // create a small dialog window and run a modal loop
    const wchar_t PCLS[] = L"GS3D_ParamsDlg";
    WNDCLASSW wc = {0}; wc.lpfnWndProc = ParamsWndProc; wc.hInstance = GetModuleHandle(NULL); wc.lpszClassName = PCLS; wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    RegisterClassW(&wc);
    // make dialog slightly taller to fit controls cleanly
    HWND dlg = CreateWindowExW(0, PCLS, L"3D Parameters", WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU, CW_USEDEFAULT, CW_USEDEFAULT, 360, 380, hwnd, NULL, GetModuleHandle(NULL), NULL);
    if (!dlg) { // non-blocking: log error, avoid messagebox
        const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "ERROR: Failed to create params dialog\n"); fclose(lf); }
        return; }
    // center over parent (account for increased height)
    RECT rc; GetWindowRect(hwnd, &rc); int cx = (rc.left+rc.right)/2 - 180; int cy = (rc.top+rc.bottom)/2 - 190; SetWindowPos(dlg, HWND_TOP, cx, cy, 360, 380, SWP_SHOWWINDOW);
    ShowWindow(dlg, SW_SHOW);
    // modal loop with Enter handling: pressing Enter will trigger OK (ID 201)
    MSG msg;
    while (IsWindow(dlg) && GetMessageW(&msg, NULL, 0,0)) {
        if (msg.message == WM_KEYDOWN && msg.wParam == VK_RETURN) {
            PostMessageW(dlg, WM_COMMAND, MAKEWPARAM(201,0), 0);
        }
        if (!IsDialogMessage(dlg, &msg)) { TranslateMessage(&msg); DispatchMessageW(&msg); }
        if (!IsWindow(dlg)) break;
    }
}
