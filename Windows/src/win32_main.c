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

static Model* g_model = NULL; static float s_ah = 30.0f, s_av = 20.0f, s_aw = 0.0f, s_dist = 300.0f; // angle_h, angle_v, angle_w, distance (match GS3Dp defaults)
static int s_dump_on_load = 0; // if set, dump face equations CSV after async model load
static ObsVertex* g_obs = NULL;
static int* g_order = NULL;
static int g_wireframe = 0;

// Forward declarations for functions used by WndProc
static void on_file_open(HWND hwnd);
static void on_show_params(HWND hwnd);
static void load_model_from_path(HWND hwnd, const char* path);

// Async loader hooks (WM_APP messages and helper types)
#define WM_MODEL_LOADED (WM_APP + 1)
#define WM_MODEL_LOAD_ERROR (WM_APP + 2)
#define WM_MODEL_LOAD_LOG (WM_APP + 3)

typedef struct { Model* m; int* order; } LoadResult;
typedef struct { HWND hwnd; char* path; } LoaderArg;

// Loading debug dialog globals
static HWND g_load_dbg = NULL; static HWND g_dbg_edit = NULL; static HWND g_dbg_ok = NULL;

static DWORD WINAPI loader_thread_fn(LPVOID arg);
static void start_async_load(HWND hwnd, const char* path);
static void create_load_debug(HWND parent, const char* path);
static void close_load_debug(void);

static void compute_projection_and_order(HWND hwnd, int* out_winw, int* out_winh, float* out_cx, float* out_cy, float* out_scale, float* out_pxmin, float* out_pxmax, float* out_pymin, float* out_pymax) {
    RECT r; GetClientRect(hwnd, &r); int winw = r.right - r.left; int winh = r.bottom - r.top; *out_winw = winw; *out_winh = winh;
    compute_obs_vertices(g_model, g_obs);
    float pxmin=1e30f, pxmax=-1e30f, pymin=1e30f, pymax=-1e30f;
    for (int i=0;i<g_model->vert_count;i++) {
        ObsVertex ov = g_obs[i]; if (ov.zo == 0.0f) continue; float px = ov.xo/ov.zo; float py = ov.yo/ov.zo; if (px<pxmin) pxmin=px; if (px>pxmax) pxmax=px; if (py<pymin) pymin=py; if (py>pymax) pymax=py; }
    float cx = (pxmin + pxmax) * 0.5f; float cy = (pymin + pymax) * 0.5f; float sdx = (pxmax - pxmin); float sdy = (pymax - pymin); float scale = 200.0f; if (sdx>1e-6f && sdy>1e-6f) { float sx = (winw*0.8f) / sdx; float sy = (winh*0.8f) / sdy; scale = fminf(sx, sy); }
    *out_cx = cx; *out_cy = cy; *out_scale = scale; *out_pxmin = pxmin; *out_pxmax = pxmax; *out_pymin = pymin; *out_pymax = pymax;
    // inform painter module of projection params (so it can compute per-face 2D bbox like GS3Dp)
    extern void set_projection_params(float cx, float cy, float scale);
    set_projection_params(cx, cy, scale);
    compute_painter_order(g_model, g_order);
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

    // TEST DRAW: draw simple shapes to verify GDI rendering
    {
        // red rectangle
        HPEN pen = CreatePen(PS_SOLID, 2, RGB(255,0,0)); HPEN oldp = SelectObject(memdc, pen);
        HBRUSH brush = CreateSolidBrush(RGB(255,0,0)); HBRUSH oldb = SelectObject(memdc, brush);
        Rectangle(memdc, 20, 20, 120, 120);
        SelectObject(memdc, oldb); DeleteObject(brush);
        SelectObject(memdc, oldp); DeleteObject(pen);
        // green circle
        HPEN pen2 = CreatePen(PS_SOLID, 2, RGB(0,255,0)); HPEN oldp2 = SelectObject(memdc, pen2);
        HBRUSH brush2 = CreateSolidBrush(RGB(0,255,0)); HBRUSH oldb2 = SelectObject(memdc, brush2);
        Ellipse(memdc, 140, 20, 240, 120);
        SelectObject(memdc, oldb2); DeleteObject(brush2);
        SelectObject(memdc, oldp2); DeleteObject(pen2);
        // text
        SetTextColor(memdc, RGB(255,255,255)); SetBkMode(memdc, TRANSPARENT);
        TextOutA(memdc, 20, 130, "GDI test shapes", 15);
        // log
        const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "TESTDRAW: drew rectangle/circle/text\n"); fclose(lf); }
    }

    // diagnostics: write projection bbox/scale and sample vertices to log
    {
        const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) {
            fprintf(lf, "PROJ: cx=%.6f cy=%.6f pxmin=%.6f pxmax=%.6f pymin=%.6f pymax=%.6f scale=%.6f win=%d,%d\n", cx, cy, pxmin, pxmax, pymin, pymax, scale, winw, winh);
            int n = (g_model->vert_count<6)?g_model->vert_count:6; for (int i=0;i<n;i++) {
                ObsVertex v = g_obs[i]; float px = (v.zo==0.0f)?v.xo:(v.xo/v.zo); float py = (v.zo==0.0f)?v.yo:(v.yo/v.zo); fprintf(lf, "VERT[%d] xo=%.6f yo=%.6f zo=%.6f px=%.6f py=%.6f\n", i, v.xo, v.yo, v.zo, px, py);
            }
            fclose(lf);
        }
    }

    // draw faces
    for (int fi=0; fi<g_model->face_count; fi++) {
        int fidx = g_order[fi]; Face* f = &g_model->faces[fidx]; POINT pts[256];
        for (int k=0;k<f->count;k++) {
            int vi = f->indices[k]; ObsVertex ov = g_obs[vi]; float px = (ov.zo==0.0f)?ov.xo:(ov.xo/ov.zo); float py = (ov.zo==0.0f)?ov.yo:(ov.yo/ov.zo);
            project_to_screen(px, py, winw, winh, scale, cx, cy, &pts[k]);
        }
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

    // draw parameter overlay at bottom (Angle H, V, W, Distance)
    {
        char buf[256]; snprintf(buf, sizeof(buf), "H=%.1f V=%.1f W=%.1f D=%.3f", s_ah, s_av, s_aw, s_dist);
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
            // also write to TEMP log file
            const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "%s", msg); fclose(lf); }
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
            else if (wParam == 'D') { if (g_model) { dumpFaceEquationsCSV_Model(g_model); const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "DEBUGDUMP: wrote %s\\equ_windows.csv\n", (tmp?tmp:".") ); fclose(lf); } }
            }
            if (changed) {
                // keep angles modulo 360
                s_ah = normalize_angle360(s_ah); s_av = normalize_angle360(s_av); s_aw = normalize_angle360(s_aw);
                set_observer_params(s_ah, s_av, s_aw, s_dist);
                if (g_model && g_obs) {
                    compute_obs_vertices(g_model, g_obs);
                    if (g_order) compute_painter_order(g_model, g_order);
                }
                InvalidateRect(hwnd, NULL, TRUE);
                // force immediate render so user sees changes at once
                render_frame(hwnd);
            }
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
            // log and redraw
            const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "ASYNC LOAD: Loaded model verts=%d faces=%d\n", g_model->vert_count, g_model->face_count); fclose(lf); }
            InvalidateRect(hwnd, NULL, TRUE);
            // optionally dump face equations for debugging parity with GS3Dp
            if (s_dump_on_load) { if (g_model) { dumpFaceEquationsCSV_Model(g_model); const char* tmp2 = getenv("TEMP"); char logfn2[1024]; if (tmp2) snprintf(logfn2, sizeof(logfn2), "%s\\viewer_win32.log", tmp2); else snprintf(logfn2, sizeof(logfn2), "viewer_win32.log"); FILE* lf2 = fopen(logfn2, "a"); if (lf2) { fprintf(lf2, "DEBUGDUMP: wrote %s\\equ_windows.csv\n", (tmp2?tmp2:".")); fclose(lf2); } } }
            // notify debug dialog that loading finished
            if (g_dbg_edit) {
                const char* done = "Model loaded successfully.\r\n"; int len = GetWindowTextLengthA(g_dbg_edit); SendMessageA(g_dbg_edit, EM_SETSEL, (WPARAM)len, (LPARAM)len); SendMessageA(g_dbg_edit, EM_REPLACESEL, 0, (LPARAM)done);
                if (g_dbg_ok) EnableWindow(g_dbg_ok, TRUE);
            }
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

    // If file still doesn't exist, fall back to raw argv[1]
    if (GetFileAttributesA(path) == INVALID_FILE_ATTRIBUTES) {
        wcstombs(path, argv[1], sizeof(path));
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

    // post result
    LoadResult* res = (LoadResult*)malloc(sizeof(LoadResult)); res->m = nm; res->order = order_buf;
    PostMessageW(hwnd, WM_MODEL_LOADED, 0, (LPARAM)res);
    // signal the UI that loading is finished (so debug dialog can enable OK)
    PostMessageW(hwnd, WM_MODEL_LOAD_LOG, 0, (LPARAM)_strdup("Load complete. Click OK to close this window.\r\n"));
    if (lsf) { FILE* lsf8 = fopen(loader_logfn, "a"); if (lsf8) { fprintf(lsf8, "THREAD: posted WM_MODEL_LOADED and finished\n"); fclose(lsf8); } }
    free(a);
    return 0;
}

static void start_async_load(HWND hwnd, const char* path) {
    if (!path) return; 
    // log and show debug dialog
    const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf0 = fopen(logfn, "a"); if (lf0) { fprintf(lf0, "UI: start_async_load called for %s\n", path); fclose(lf0); }
    create_load_debug(hwnd, path);
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
        wchar_t buf[64]; swprintf(buf, 64, L"%.0f", s_ah); hEditAH = CreateWindowExW(0, L"EDIT", buf, WS_CHILD|WS_VISIBLE|WS_BORDER|ES_LEFT, 120,20,140,22,dlg,(HMENU)101,GetModuleHandle(NULL),NULL);
        swprintf(buf, 64, L"%.0f", s_av); hEditAV = CreateWindowExW(0, L"EDIT", buf, WS_CHILD|WS_VISIBLE|WS_BORDER|ES_LEFT, 120,55,140,22,dlg,(HMENU)102,GetModuleHandle(NULL),NULL);
        swprintf(buf, 64, L"%.0f", s_aw); hEditAW = CreateWindowExW(0, L"EDIT", buf, WS_CHILD|WS_VISIBLE|WS_BORDER|ES_LEFT, 120,90,140,22,dlg,(HMENU)103,GetModuleHandle(NULL),NULL);
        swprintf(buf, 64, L"%.3f", s_dist); hEditDist = CreateWindowExW(0, L"EDIT", buf, WS_CHILD|WS_VISIBLE|WS_BORDER|ES_LEFT, 120,125,140,22,dlg,(HMENU)104,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"STATIC", L"Angle H:", WS_CHILD|WS_VISIBLE, 20,20,90,22,dlg,NULL,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"STATIC", L"Angle V:", WS_CHILD|WS_VISIBLE, 20,55,90,22,dlg,NULL,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"STATIC", L"Angle W:", WS_CHILD|WS_VISIBLE, 20,90,90,22,dlg,NULL,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"STATIC", L"Distance:", WS_CHILD|WS_VISIBLE, 20,125,90,22,dlg,NULL,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"BUTTON", L"OK", WS_CHILD|WS_VISIBLE|BS_DEFPUSHBUTTON, 40,160,80,26,dlg,(HMENU)201,GetModuleHandle(NULL),NULL);
        CreateWindowExW(0, L"BUTTON", L"Cancel", WS_CHILD|WS_VISIBLE, 160,160,80,26,dlg,(HMENU)202,GetModuleHandle(NULL),NULL);
        return 0;
    }
    if (msg == WM_COMMAND) {
        int id = LOWORD(wParam);
        if (id == 201) {
            wchar_t buf[64]; GetWindowTextW(hEditAH, buf, 64); s_ah = (float)wcstod(buf, NULL);
            GetWindowTextW(hEditAV, buf, 64); s_av = (float)wcstod(buf, NULL);
            GetWindowTextW(hEditAW, buf, 64); s_aw = (float)wcstod(buf, NULL);
            GetWindowTextW(hEditDist, buf, 64); s_dist = (float)wcstod(buf, NULL);
            // normalize angles to [0,360)
            s_ah = normalize_angle360(s_ah); s_av = normalize_angle360(s_av); s_aw = normalize_angle360(s_aw);
            // apply params and recompute using GS3Dp semantics
            set_observer_params(s_ah, s_av, s_aw, s_dist);
            if (g_model && g_obs) {
                compute_obs_vertices(g_model, g_obs);
                if (g_order) compute_painter_order(g_model, g_order);
            }
            HWND owner = GetWindow(dlg, GW_OWNER);
            if (owner) {
                InvalidateRect(owner, NULL, TRUE);
                render_frame(owner);
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
    HWND dlg = CreateWindowExW(0, PCLS, L"3D Parameters", WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU, CW_USEDEFAULT, CW_USEDEFAULT, 320, 240, hwnd, NULL, GetModuleHandle(NULL), NULL);
    if (!dlg) { // non-blocking: log error, avoid messagebox
        const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_win32.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_win32.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "ERROR: Failed to create params dialog\n"); fclose(lf); }
        return; }
    // center over parent
    RECT rc; GetWindowRect(hwnd, &rc); int cx = (rc.left+rc.right)/2 - 160; int cy = (rc.top+rc.bottom)/2 - 120; SetWindowPos(dlg, HWND_TOP, cx, cy, 320, 240, SWP_SHOWWINDOW);
    ShowWindow(dlg, SW_SHOW);
    // modal loop
    MSG msg; while (IsWindow(dlg) && GetMessageW(&msg, NULL, 0,0)) { if (!IsDialogMessage(dlg, &msg)) { TranslateMessage(&msg); DispatchMessageW(&msg); } if (!IsWindow(dlg)) break; }
}
