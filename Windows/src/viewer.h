#ifndef VIEWER_H
#define VIEWER_H

#include "fixed32.h"

typedef struct { float x,y,z; } Vertex;
typedef struct { int *indices; int count; float z_min,z_mean,z_max; float plane_a, plane_b, plane_c, plane_d; int minx,maxx,miny,maxy; int display_flag; } Face;
typedef struct { Vertex *verts; int vert_count; Face *faces; int face_count; } Model;

Model* load_obj(const char* path);
void free_model(Model* m);

void set_observer_params(float ah, float av, float aw, float dist);
int compute_painter_order(Model* m, int* order_out);
void dumpFaceEquationsCSV_Model(Model* m);
typedef struct { float xo, yo, zo; } ObsVertex;
void compute_obs_vertices(Model* m, ObsVertex* out); // out = ObsVertex*

void screen_coords_from_proj(float px, float py, int winw, int winh, float scale, float centerx, float centery, int* sx, int* sy);

#ifdef USE_SDL
#include <SDL.h>
void draw_filled_polygon(SDL_Renderer* rend, float* xs, float* ys, int n, unsigned char r, unsigned char g, unsigned char b);
void draw_polygon_outline(SDL_Renderer* rend, float* xs, float* ys, int n, Uint8 r, Uint8 g, Uint8 b);
void save_screenshot(SDL_Window* win, SDL_Renderer* rend, const char* path);
#endif

#endif // VIEWER_H
