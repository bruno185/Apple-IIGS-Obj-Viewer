#ifndef VIEWER_H
#define VIEWER_H


typedef struct { float x,y,z; } Vertex;
typedef struct { int *indices; int count; float z_min,z_mean,z_max; float plane_a, plane_b, plane_c, plane_d; int minx,maxx,miny,maxy; int display_flag; } Face;
typedef struct { Vertex *verts; int vert_count; Face *faces; int face_count; } Model;

Model* load_obj(const char* path);
void free_model(Model* m);

void set_observer_params(float ah, float av, float aw, float dist);
int compute_painter_order(Model* m, int* order_out);

/* Inconclusive pairs accessors (diagnostic) */
int get_inconclusive_pair_count(void);
int* get_inconclusive_pairs(void);
void clear_inconclusive_pairs(void);
int dump_pairwise_debug(Model* m, int version, const char* outpath);
extern int g_painter_order_version; /* 1 = V1 (default), 2 = V2 */
void dumpFaceEquationsCSV_Model(Model* m);
typedef struct { float xo, yo, zo; } ObsVertex;
void compute_obs_vertices(Model* m, ObsVertex* out); // out = ObsVertex*

// Painter modes (match GS3Dp): FAST, FIXED, FLOAT
#define PAINTER_MODE_FAST 0
#define PAINTER_MODE_FIXED 1
#define PAINTER_MODE_FLOAT 2
extern int g_painter_mode;
void set_painter_mode(int mode);

void screen_coords_from_proj(float px, float py, int winw, int winh, float scale, float centerx, float centery, int* sx, int* sy);

/* Windows-specific parity helpers added to align with GS3Dp features */
void set_cull_back_faces(int v);
int get_cull_back_faces(void);
int projected_polygons_overlap(Model* m, int f1, int f2);
void inspect_polygons_overlap(Model* m, void* params, const char* filename);
void inspect_faces_before(Model* m);
void inspect_faces_after(Model* m);
void display_model_face_ids(Model* m);
void set_highlight_faces(int* faces, int count, int color);
void clear_highlight_faces(void);
int face_is_highlighted(int f);

#ifdef USE_SDL
#include <SDL.h>
void draw_filled_polygon(SDL_Renderer* rend, float* xs, float* ys, int n, unsigned char r, unsigned char g, unsigned char b);
void draw_polygon_outline(SDL_Renderer* rend, float* xs, float* ys, int n, Uint8 r, Uint8 g, Uint8 b);
void save_screenshot(SDL_Window* win, SDL_Renderer* rend, const char* path);
#endif

#endif // VIEWER_H
