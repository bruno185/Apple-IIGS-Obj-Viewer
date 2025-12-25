#include <SDL.h>
#include <stdio.h>
#include <stdlib.h>

// Simple scanline polygon fill for convex polygons; vertices are in screen coords (float)

void draw_filled_polygon(SDL_Renderer* rend, float* xs, float* ys, int n, Uint8 r, Uint8 g, Uint8 b) {
    if (n < 3) return; // compute ymin/ymax
    float ymin = ys[0], ymax = ys[0]; for (int i=1;i<n;i++) { if (ys[i] < ymin) ymin = ys[i]; if (ys[i] > ymax) ymax = ys[i]; }
    int y0 = (int)floorf(ymin), y1 = (int)ceilf(ymax);
    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(rend, r, g, b, 255);
    for (int y = y0; y <= y1; y++) {
        float py = (float)y + 0.5f;
        // find intersections with edges
        int cnt = 0; float inter[128];
        for (int i=0;i<n;i++) {
            int j = (i+1)%n; float y0f = ys[i], y1f = ys[j]; if ((py < fminf(y0f,y1f)) || (py > fmaxf(y0f,y1f))) continue; float x0f = xs[i], x1f = xs[j]; if (fabsf(y1f-y0f) < 1e-6f) continue; float t = (py - y0f)/(y1f - y0f); float xi = x0f + t*(x1f - x0f); inter[cnt++] = xi; if (cnt >= 128) break; }
        if (cnt < 2) continue; // sort
        for (int i=0;i<cnt-1;i++) for (int j=i+1;j<cnt;j++) if (inter[j] < inter[i]) { float tmp = inter[i]; inter[i]=inter[j]; inter[j]=tmp; }
        for (int i=0;i<cnt;i+=2) {
            int xstart = (int)floorf(inter[i]); int xend = (int)ceilf(inter[i+1]); SDL_RenderDrawLine(rend, xstart, y, xend, y);
        }
    }
}

void screen_coords_from_proj(float px, float py, int winw, int winh, float scale, float centerx, float centery, int* sx, int* sy) {
    *sx = (int)(winw*0.5f + (px - centerx) * scale);
    *sy = (int)(winh*0.5f - (py - centery) * scale);
}

void draw_polygon_outline(SDL_Renderer* rend, float* xs, float* ys, int n, Uint8 r, Uint8 g, Uint8 b) {
    if (n < 2) return;
    SDL_SetRenderDrawColor(rend, r, g, b, 255);
    for (int i=0;i<n;i++) {
        int j = (i+1)%n;
        SDL_RenderDrawLine(rend, (int)xs[i], (int)ys[i], (int)xs[j], (int)ys[j]);
    }
}

void save_screenshot(SDL_Window* win, SDL_Renderer* rend, const char* path) {
    int w,h; SDL_GetRendererOutputSize(rend, &w, &h);
    SDL_Surface* surf = SDL_CreateRGBSurfaceWithFormat(0, w, h, 32, SDL_PIXELFORMAT_ARGB8888);
    if (!surf) return;
    SDL_RenderReadPixels(rend, NULL, SDL_PIXELFORMAT_ARGB8888, surf->pixels, surf->pitch);
    SDL_SaveBMP(surf, path);
    SDL_FreeSurface(surf);
}
