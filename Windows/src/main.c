#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <SDL.h>
#include "viewer.h"



int main(int argc, char** argv) {
    if (argc < 2) { printf("Usage: viewer <path.obj> [--angle_h <deg>] [--angle_v <deg>] [--angle_w <deg>] [--distance <dist>]\n"); return 1; }
    const char* path = argv[1]; float ah=30.0f, av=20.0f, aw=0.0f, dist=300.0f;
    for (int i=2;i<argc;i++) { if (strcmp(argv[i],"--angle_h")==0 && i+1<argc) ah = atof(argv[++i]); else if (strcmp(argv[i],"--angle_v")==0 && i+1<argc) av = atof(argv[++i]); else if (strcmp(argv[i],"--angle_w")==0 && i+1<argc) aw = atof(argv[++i]); else if (strcmp(argv[i],"--distance")==0 && i+1<argc) dist = atof(argv[++i]); }
    set_observer_params(ah, av, aw, dist);

    Model* m = load_obj(path); if (!m) { printf("Failed to load %s\n", path); return 1; }
    int* order = malloc(sizeof(int)*m->face_count);
    compute_painter_order(m, order);

    ObsVertex* obs = malloc(sizeof(ObsVertex) * m->vert_count);

    if (SDL_Init(SDL_INIT_VIDEO) != 0) { printf("SDL_Init Error: %s\n", SDL_GetError()); return 1; }
    // quick write test to check file permissions
    {
        FILE* tf = fopen("F:\\Bruno\\Dev\\AppleWin\\Projets\\ORCA\\OBJ Viewer Review\\Dony\\viewer_write_test.txt","w");
        if (tf) { fprintf(tf, "write test\n"); fclose(tf); } else { FILE* ef = fopen("%TEMP%\\viewer_write_test_failed.txt","w"); if (ef) { fprintf(ef, "failed\n"); fclose(ef); } }
    }
    SDL_Window* win = SDL_CreateWindow("GS3Dp Viewer", 100,100,1024,768, SDL_WINDOW_SHOWN);
    SDL_Renderer* rend = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED|SDL_RENDERER_PRESENTVSYNC);

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
                else if (k == SDLK_w) { wireframe = !wireframe; }
                else if (k == SDLK_s) { save_screenshot(win, rend, "viewer_capture.bmp"); }
            }
        }
        if (camera_changed) set_observer_params(h, ah, av, dist);
        compute_obs_vertices(m, obs);
        // compute autoscale / centering
        float pxmin=1e30f, pxmax=-1e30f, pymin=1e30f, pymax=-1e30f;
        for (int i=0;i<m->vert_count;i++) {
            ObsVertex ov = obs[i]; if (ov.yo == 0.0f) continue; float px = ov.xo / ov.yo; float py = ov.zo / ov.yo; if (px<pxmin) pxmin=px; if (px>pxmax) pxmax=px; if (py<pymin) pymin=py; if (py>pymax) pymax=py; }
        float cx = (pxmin + pxmax) * 0.5f; float cy = (pymin + pymax) * 0.5f; float sdx = (pxmax - pxmin); float sdy = (pymax - pymin); float scale = 200.0f; if (sdx>1e-6f && sdy>1e-6f) { float sx = (winw*0.8f) / sdx; float sy = (winh*0.8f) / sdy; scale = fminf(sx, sy); }
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
            if (df) { fprintf(df, "cx=%.6f,cy=%.6f,pxmin=%.6f,pxmax=%.6f,pymin=%.6f,pymax=%.6f,scale=%.6f\n", cx, cy, pxmin, pxmax, pymin, pymax, scale); fclose(df); }
        }

        SDL_SetRenderDrawColor(rend, 16,16,16,255); SDL_RenderClear(rend);
        // recompute order each frame to reflect camera changes
        compute_painter_order(m, order);
        // draw faces in order
        for (int fi=0; fi<m->face_count; fi++) {
            int fidx = order[fi]; Face* f = &m->faces[fidx]; float xs[256], ys[256];
            for (int k=0;k<f->count;k++) {
                int vi = f->indices[k]; ObsVertex ov = obs[vi]; float px = (ov.yo==0.0f)?ov.xo:(ov.xo/ov.yo); float py = (ov.yo==0.0f)?ov.zo:(ov.zo/ov.yo);
                int sx, sy; screen_coords_from_proj(px, py, winw, winh, scale, cx, cy, &sx, &sy); xs[k] = (float)sx; ys[k] = (float)sy;
            }
            // pick color
            Uint8 r = (Uint8)((fidx*37)&255), g = (Uint8)((fidx*83)&255), b=(Uint8)((fidx*191)&255);
            if (wireframe) draw_polygon_outline(rend, xs, ys, f->count, r,g,b); else draw_filled_polygon(rend, xs, ys, f->count, r,g,b);
        }
        SDL_RenderPresent(rend);
        SDL_Delay(16);
    }
    SDL_DestroyRenderer(rend); SDL_DestroyWindow(win); SDL_Quit(); free(order); free_model(m); free(obs);
    return 0;
}
