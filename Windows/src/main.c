#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <SDL.h>
#include "viewer.h"



int main(int argc, char** argv) {
    if (argc < 2) { printf("Usage: viewer <path.obj> [--angle_h <deg>] [--angle_v <deg>] [--angle_w <deg>] [--distance <dist>]\n"); return 1; }
    const char* path = argv[1]; float ah=30.0f, av=20.0f, aw=0.0f, dist=300.0f; float proj_scale = 100.0f; // proj_scale is UI-controlled (keys ',' / '.')
    for (int i=2;i<argc;i++) { if (strcmp(argv[i],"--angle_h")==0 && i+1<argc) ah = atof(argv[++i]); else if (strcmp(argv[i],"--angle_v")==0 && i+1<argc) av = atof(argv[++i]); else if (strcmp(argv[i],"--angle_w")==0 && i+1<argc) aw = atof(argv[++i]); else if (strcmp(argv[i],"--distance")==0 && i+1<argc) dist = atof(argv[++i]); }
    set_observer_params(ah, av, aw, dist); // distance is controlled via the 3D Parameters dialog or CLI --distance option

    Model* m = load_obj(path); if (!m) { printf("Failed to load %s\n", path); return 1; }
    int* order = malloc(sizeof(int)*m->face_count);
    compute_painter_order(m, order);

    ObsVertex* obs = malloc(sizeof(ObsVertex) * m->vert_count);

    // Apply brute auto-fit: distance = 3 * max_dim (computed from model bbox), then compute proj_scale so projected model fits window
    {
        float minx = 1e30f, maxx = -1e30f, miny = 1e30f, maxy = -1e30f, minz = 1e30f, maxz = -1e30f;
        for (int i=0;i<m->vert_count;i++) { float vx = m->verts[i].x, vy = m->verts[i].y, vz = m->verts[i].z; if (vx<minx) minx=vx; if (vx>maxx) maxx=vx; if (vy<miny) miny=vy; if (vy>maxy) maxy=vy; if (vz<minz) minz=vz; if (vz>maxz) maxz=vz; }
        float dx = maxx - minx; float dy = maxy - miny; float dz = maxz - minz; float max_dim = dx; if (dy > max_dim) max_dim = dy; if (dz > max_dim) max_dim = dz;
        dist = 3.0f * max_dim;
        set_observer_params(ah, av, aw, dist);
        compute_obs_vertices(m, obs);
        float pxmin=1e30f, pxmax=-1e30f, pymin=1e30f, pymax=-1e30f;
        for (int i = 0; i < m->vert_count; i++) {
            ObsVertex v = obs[i];
            if (v.yo == 0.0f) continue;
            float px = v.xo / v.yo;
            float py = v.zo / v.yo;
            if (px < pxmin) pxmin = px;
            if (px > pxmax) pxmax = px;
            if (py < pymin) pymin = py;
            if (py > pymax) pymax = py;
        }
        float margin = 0.9f; float s1 = (winw * margin) / (pxmax - pxmin); float s2 = (winh * margin) / (pymax - pymin); proj_scale = (s1 < s2) ? s1 : s2;
        snprintf(titlebuf, sizeof(titlebuf), "GS3Dp Viewer - proj_scale=%.1f", proj_scale); SDL_SetWindowTitle(win, titlebuf);
    }

    if (SDL_Init(SDL_INIT_VIDEO) != 0) { printf("SDL_Init Error: %s\n", SDL_GetError()); return 1; }
    // quick write test to check file permissions
    {
        FILE* tf = fopen("F:\\Bruno\\Dev\\AppleWin\\Projets\\ORCA\\OBJ Viewer Review\\Dony\\viewer_write_test.txt","w");
        if (tf) { fprintf(tf, "write test\n"); fclose(tf); } else { FILE* ef = fopen("%TEMP%\\viewer_write_test_failed.txt","w"); if (ef) { fprintf(ef, "failed\n"); fclose(ef); } }
    }
    SDL_Window* win = SDL_CreateWindow("GS3Dp Viewer", 100,100,1024,768, SDL_WINDOW_SHOWN);
    SDL_Renderer* rend = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED|SDL_RENDERER_PRESENTVSYNC);
    // Show current projection scale in window title
    char titlebuf[128]; snprintf(titlebuf, sizeof(titlebuf), "GS3Dp Viewer - proj_scale=%.1f", proj_scale); SDL_SetWindowTitle(win, titlebuf);

    int winw = 1024, winh = 768;
    int running = 1; SDL_Event ev; int wireframe = 0;

    while (running) {
        int camera_changed = 0;
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_QUIT) { running = 0; }
            else if (ev.type == SDL_KEYDOWN) {
                SDL_Keycode k = ev.key.keysym.sym;
                if (k == SDLK_LEFT) { ah -= 5.0f; camera_changed = 1; }
                else if (k == SDLK_RIGHT) { ah += 5.0f; camera_changed = 1; }
                else if (k == SDLK_UP) { av += 5.0f; camera_changed = 1; }
                else if (k == SDLK_DOWN) { av -= 5.0f; camera_changed = 1; }
                else if (k == SDLK_EQUALS || k == SDLK_PLUS) { dist *= 0.9f; camera_changed = 1; }
                else if (k == SDLK_MINUS) { dist *= 1.1f; camera_changed = 1; }
                /* distance-scale handling removed: distance is controlled via the Distance field in the 3D dialog */
                /* proj_scale handlers remain below */
                else if (k == SDLK_PERIOD) { proj_scale *= 1.1f; if (proj_scale > 2000.0f) proj_scale = 2000.0f; camera_changed = 1; printf("Projection scale increased: %.1f\n", proj_scale); snprintf(titlebuf, sizeof(titlebuf), "GS3Dp Viewer - proj_scale=%.1f", proj_scale); SDL_SetWindowTitle(win, titlebuf); }
                else if (k == SDLK_COMMA) { proj_scale *= 0.9f; if (proj_scale < 10.0f) proj_scale = 10.0f; camera_changed = 1; printf("Projection scale decreased: %.1f\n", proj_scale); snprintf(titlebuf, sizeof(titlebuf), "GS3Dp Viewer - proj_scale=%.1f", proj_scale); SDL_SetWindowTitle(win, titlebuf); }
                else if (k == SDLK_b) { set_cull_back_faces(!get_cull_back_faces()); printf("Back-face culling: %s\n", get_cull_back_faces() ? "ON" : "OFF"); camera_changed = 1; }
                else if (k == SDLK_d) { inspect_faces_before(m); camera_changed = 1; }
                else if (k == SDLK_o) {
                    // Inspect polygon overlap
                    printf("Enter two face IDs (f1 f2) or empty to cancel: "); char buf[128]; if (fgets(buf, sizeof(buf), stdin)) {
                        int f1=-1,f2=-1; if (sscanf(buf, "%d %d", &f1, &f2) == 2) {
                            int ov = projected_polygons_overlap(m, f1, f2);
                            printf("Projected overlap: %s\n", ov ? "YES" : "NO");
                            if (ov) {
                                printf("Press Y to preview highlighting faces, any other key to skip: "); char c = getchar(); if (c=='Y' || c=='y') { int arr[2] = {f1,f2}; set_highlight_faces(arr, 2, 10); wireframe = 1; }
                            }
                        } else { printf("Invalid input\n"); }
                    }
                    camera_changed = 1;
                }
                else if (k == SDLK_l) { display_model_face_ids(m); camera_changed = 1; }
                else if (k == SDLK_w) { wireframe = !wireframe; }
                else if (k == SDLK_s) { /* SHIFT-S will inspect_after, plain 's' saves screenshot */ if (ev.key.keysym.mod & KMOD_SHIFT) { inspect_faces_after(m); camera_changed = 1; } else { save_screenshot(win, rend, "viewer_capture.bmp"); } }
            }
        }
        if (camera_changed) set_observer_params(ah, av, aw, dist);
        compute_obs_vertices(m, obs);
        // compute autoscale / centering
        float pxmin=1e30f, pxmax=-1e30f, pymin=1e30f, pymax=-1e30f;
        for (int i=0;i<m->vert_count;i++) {
            ObsVertex ov = obs[i]; if (ov.yo == 0.0f) continue; float px = ov.xo / ov.yo; float py = ov.zo / ov.yo; if (px<pxmin) pxmin=px; if (px>pxmax) pxmax=px; if (py<pymin) pymin=py; if (py>pymax) pymax=py; }
        // Use projection center (0,0) and user-controlled projection scale (proj_scale) so distance modifies perspective
        float cx = 0.0f; float cy = 0.0f; float scale = proj_scale;
        // keep diagnostics: compute bbox spans for logging but do not use them to autoscale
        float sdx = (pxmax - pxmin); float sdy = (pymax - pymin);
        // Write projection diagnostics when camera changed (or on first frame) for debugging
        static int _viewer_proj_written = 0;
        if (camera_changed || !_viewer_proj_written) {
            _viewer_proj_written = 1;
            FILE* df = fopen("viewer_projection.csv","a");
            if (df) { printf("PROJ: opened viewer_projection.csv (cwd)\n"); }
            if (!df) {
                const char* alt = "F:\\Bruno\\Dev\\AppleWin\\Projets\\ORCA\\OBJ Viewer Review\\Dony\\viewer_projection.csv";
                df = fopen(alt, "a"); if (df) { printf("PROJ: opened alt path %s\n", alt); }
            }
            if (!df) {
                const char* tmp = getenv("TEMP");
                if (tmp) { char tmpfn[1024]; snprintf(tmpfn, sizeof(tmpfn), "%s\\viewer_projection.csv", tmp); df = fopen(tmpfn, "a"); if (df) { printf("PROJ: opened TEMP path %s\n", tmpfn); } }
            }
            if (!df) { printf("PROJ: failed to open any path for viewer_projection.csv\n"); }
            if (df) { fprintf(df, "cx=%.6f,cy=%.6f,pxmin=%.6f,pxmax=%.6f,pymin=%.6f,pymax=%.6f,scale=%.6f,proj_scale=%.6f,dist=%.6f\n", cx, cy, pxmin, pxmax, pymin, pymax, scale, proj_scale, dist); fclose(df); }
        }

        // Update painter projection params and recompute order each frame to reflect camera or projection changes
        extern void set_projection_params(float cx, float cy, float scale);
        set_projection_params(cx, cy, scale);
        compute_painter_order(m, order);

        // draw faces in order
        for (int fi=0; fi<m->face_count; fi++) {
            int fidx = order[fi]; Face* f = &m->faces[fidx]; float xs[256], ys[256];
            for (int k=0;k<f->count;k++) {
                int vi = f->indices[k]; ObsVertex ov = obs[vi]; float px = (ov.yo==0.0f)?ov.xo:(ov.xo/ov.yo); float py = (ov.yo==0.0f)?ov.zo:(ov.zo/ov.yo); int sx, sy; screen_coords_from_proj(px, py, winw, winh, scale, cx, cy, &sx, &sy); xs[k] = (float)sx; ys[k] = (float)sy;
            }
            Uint8 r = (Uint8)((fidx*37)&255), g = (Uint8)((fidx*83)&255), b=(Uint8)((fidx*191)&255);
            if (wireframe) draw_polygon_outline(rend, xs, ys, f->count, r,g,b); else draw_filled_polygon(rend, xs, ys, f->count, r,g,b);
            /* overlay highlights (inspectors) */
            int hcol = face_is_highlighted(fidx);
            if (hcol) {
                Uint8 hr = (hcol==6)?255:0;
                Uint8 hg = (hcol==6)?165:255;
                Uint8 hb = (hcol==6)?0:192;
                draw_polygon_outline(rend, xs, ys, f->count, hr, hg, hb);
            }
}
