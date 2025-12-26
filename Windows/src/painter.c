#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <string.h>
#include "fixed32.h"

// Minimal structures (mirrors those in obj.c)
typedef struct { float x,y,z; } Vertex;
// Extended Face struct to match GS3Dp: add plane coeffs, bbox and display flag
typedef struct { int *indices; int count; float z_min,z_mean,z_max; float plane_a, plane_b, plane_c, plane_d; int minx, maxx, miny, maxy; int display_flag; } Face;
typedef struct { Vertex* verts; int vert_count; Face* faces; int face_count; } Model;

// Observer params (match GS3Dp semantics: angle_h, angle_v, angle_w, distance)
static int s_angle_h = 30; static int s_angle_v = 20; static int s_angle_w = 0; static float s_distance = 300.0f;

void set_observer_params(float ah, float av, float aw, float dist) { s_angle_h = (int)roundf(ah); s_angle_v = (int)roundf(av); s_angle_w = (int)roundf(aw); s_distance = dist; }

// Projection parameters set by win32_main to ensure bbox uses same cx/cy/scale
static float s_proj_cx = 0.0f, s_proj_cy = 0.0f, s_proj_scale = 200.0f;
void set_projection_params(float cx, float cy, float scale) { s_proj_cx = cx; s_proj_cy = cy; s_proj_scale = scale; }

// Transform model into observer (GS3Dp float port of processModelFast per-vertex transform)
typedef struct { float xo, yo, zo; } ObsVertex;

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
    float distance = s_distance;
    for (int i = 0; i < m->vert_count; ++i) {
        float x = m->verts[i].x;
        float y = m->verts[i].y;
        float z = m->verts[i].z;
        float term1 = x * cos_h_cos_v;
        float term2 = y * sin_h_cos_v;
        float term3 = z * sin_v;
        float zo = -term1 - term2 - term3 + distance;
        out[i].zo = zo;
        if (zo > 0.0f) {
            float xo = -x * sin_h + y * cos_h;
            float yo = -x * cos_h_sin_v - y * sin_h_sin_v + z * cos_v;
            out[i].xo = xo; out[i].yo = yo;
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

// Pair comparator helpers and GS3Dp-like calculateFaceDepths + insertion-based Newell/Sancha
static Model* g_model = NULL; static ObsVertex* g_obsv = NULL; // for comparator

static int point_in_triangle(float px, float py, float Ax, float Ay, float Bx, float By, float Cx, float Cy) {
    float v0x = Cx - Ax, v0y = Cy - Ay; float v1x = Bx - Ax, v1y = By - Ay; float v2x = px - Ax, v2y = py - Ay;
    float den = v0x * v1y - v1x * v0y; if (fabsf(den) < 1e-9f) return 0;
    float u = (v2x * v1y - v1x * v2y) / den; float v = (v0x * v2y - v2x * v0y) / den; return (u>=0 && v>=0 && (u+v)<=1.0f);
}

// Quick comparator for initial qsort by z_mean descending (stable tie-breaker by index)
static int compar_face_qsort(const void* pa, const void* pb) {
    int a = *(const int*)pa; int b = *(const int*)pb;
    float za = g_model->faces[a].z_mean; float zb = g_model->faces[b].z_mean;
    if (za > zb) return -1; if (za < zb) return 1; if (a < b) return -1; if (a > b) return 1; return 0;
}

static int face_depth_at_point_tri(int face_idx, float px, float py, float* depth) {
    Face* f = &g_model->faces[face_idx]; if (f->count < 3) return 0;
    for (int k=1;k<f->count-1;k++) {
        int a=f->indices[0], b=f->indices[k], c=f->indices[k+1]; ObsVertex A=g_obsv[a], B=g_obsv[b], C=g_obsv[c]; float Ax,Ay,Bx,By,Cx,Cy; project_vertex(&A,&Ax,&Ay); project_vertex(&B,&Bx,&By); project_vertex(&C,&Cx,&Cy); if (!point_in_triangle(px,py,Ax,Ay,Bx,By,Cx,Cy)) continue;
        float denom = (By - Cy)*(Ax - Cx) + (Cx - Bx)*(Ay - Cy); if (fabsf(denom) < 1e-9f) return 0; float wA = ((By - Cy)*(px - Cx) + (Cx - Bx)*(py - Cy))/denom; float wB = ((Cy - Ay)*(px - Cx) + (Ax - Cx)*(py - Cy))/denom; float wC = 1.0f - wA - wB; *depth = wA*A.yo + wB*B.yo + wC*C.yo; return 1; }
    return 0;
}

static void calculateFaceDepths(Model* model) {
    if (!model) return;
    for (int i = 0; i < model->face_count; ++i) {
        Face* face = &model->faces[i]; float zmin = 1e30f, zmax = -1e30f, sum = 0.0f; int n = face->count; int minx = 999999, maxx = -999999, miny = 999999, maxy = -999999; int display_flag = 1;
        for (int k = 0; k < n; ++k) {
            int vid = face->indices[k]; if (vid < 0 || vid >= model->vert_count) continue;
            ObsVertex ov = g_obsv[vid]; float zo = ov.zo; if (zo < 0.0f) display_flag = 0; if (zo < zmin) zmin = zo; if (zo > zmax) zmax = zo; sum += zo;
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
                ObsVertex A = g_obsv[idx0], B = g_obsv[idx1], C = g_obsv[idx2]; float x1=A.xo, y1=A.yo, z1=A.zo; float x2=B.xo, y2=B.yo, z2=B.zo; float x3=C.xo, y3=C.yo, z3=C.zo;
                float a = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2);
                float b = -x1*(z2 - z3) + x2*(z1 - z3) - x3*(z1 - z2);
                float c = x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2);
                float t1 = y2*z3 - y3*z2; float t2 = y1*z3 - y3*z1; float t3 = y1*z2 - y2*z1;
                float d = -x1 * t1 + x2 * t2 - x3 * t3;
                face->plane_a = a; face->plane_b = b; face->plane_c = c; face->plane_d = d;
            }
        }
        face->z_min = (n>0)?zmin:0.0f; face->z_max = (n>0)?zmax:0.0f; face->z_mean = (n>0)?(sum/n):0.0f; face->minx = (n>0)?minx:0; face->maxx = (n>0)?maxx:0; face->miny = (n>0)?miny:0; face->maxy = (n>0)?maxy:0; face->display_flag = display_flag;
    }
}

static int pair_should_swap(Model* model, int f1, int f2) {
    Face* fA = &model->faces[f1]; Face* fB = &model->faces[f2];
    if (fB->z_max <= fA->z_min) return -1; if (fA->z_max <= fB->z_min) return 1;
    if (fA->maxx <= fB->minx || fB->maxx <= fA->minx) return -1; if (fA->maxy <= fB->miny || fB->maxy <= fA->miny) return -1;
    float a1=fA->plane_a, b1=fA->plane_b, c1=fA->plane_c, d1=fA->plane_d;
    float a2=fB->plane_a, b2=fB->plane_b, c2=fB->plane_c, d2=fB->plane_d;
    int n1 = fA->count, n2 = fB->count;
    float maxabs1 = fabsf(d1);
    for (int k=0;k<n2;k++) { int v = fB->indices[k]; float tv = a1*g_obsv[v].xo + b1*g_obsv[v].yo + c1*g_obsv[v].zo + d1; if (fabsf(tv) > maxabs1) maxabs1 = fabsf(tv); }
    float eps_rel1 = maxabs1 * 1e-6f; if (eps_rel1 < 0.0001f) eps_rel1 = 0.0001f; int obs1 = 0; if (d1 > eps_rel1) obs1 = 1; else if (d1 < -eps_rel1) obs1 = -1;
    if (obs1 != 0) {
        int pos=0, neg=0; for (int k=0;k<n2;k++) { int v=fB->indices[k]; float tv = a1*g_obsv[v].xo + b1*g_obsv[v].yo + c1*g_obsv[v].zo + d1; if (fabsf(tv) <= eps_rel1) ; else if (tv>0) pos++; else neg++; }
        int thr = (3*n2 + 3)/4; if ((obs1==1 && pos>=thr) || (obs1==-1 && neg>=thr)) return -1;
    }
    float maxabs2 = fabsf(d2);
    for (int k=0;k<n1;k++) { int v = fA->indices[k]; float tv = a2*g_obsv[v].xo + b2*g_obsv[v].yo + c2*g_obsv[v].zo + d2; if (fabsf(tv) > maxabs2) maxabs2 = fabsf(tv); }
    float eps_rel2 = maxabs2 * 1e-6f; if (eps_rel2 < 0.0001f) eps_rel2 = 0.0001f; int obs2=0; if (d2>eps_rel2) obs2=1; else if (d2 < -eps_rel2) obs2 = -1;
    if (obs2 != 0) { int pos2=0, neg2=0; for (int k=0;k<n1;k++){ int v=fA->indices[k]; float tv = a2*g_obsv[v].xo + b2*g_obsv[v].yo + c2*g_obsv[v].zo + d2; if (fabsf(tv) <= eps_rel2) ; else if (tv>0) pos2++; else neg2++; } int thr2=(3*n1+3)/4; if ((obs2==1 && neg2>=thr2) || (obs2==-1 && pos2>=thr2)) return -1; }
    if (obs1 == 0) {
        int pos=0, neg=0; for (int k=0;k<n2;k++){ int v=fB->indices[k]; float tv=a1*g_obsv[v].xo + b1*g_obsv[v].yo + c1*g_obsv[v].zo + d1; if (fabsf(tv)<=eps_rel1) ; else if (tv>0) pos++; else neg++; } int thr=(3*n2+3)/4; if ((obs1==1 && neg>=thr) || (obs1==-1 && pos>=thr)) return 1;
    }
    if (obs2 == 0) { int pos2=0, neg2=0; for (int k=0;k<n1;k++){ int v=fA->indices[k]; float tv=a2*g_obsv[v].xo + b2*g_obsv[v].yo + c2*g_obsv[v].zo + d2; if (fabsf(tv)<=eps_rel2) ; else if (tv>0) pos2++; else neg2++; } int thr2=(3*n1+3)/4; if ((obs2==1 && pos2>=thr2) || (obs2==-1 && neg2>=thr2)) return 1; }
    return 0;
}

int compute_painter_order(Model* m, int* order_out) {
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
    for (int i=0;i<m->face_count;i++) order_out[i] = i;
    qsort(order_out, m->face_count, sizeof(int), compar_face_qsort);
    int face_count = m->face_count;
    int swapped = 0;
    int ordered_pairs_capacity = face_count * 4;
    typedef struct { int face1; int face2; } OrderedPair;
    OrderedPair* ordered_pairs = NULL; int ordered_pairs_count = 0;
    if (ordered_pairs_capacity > 0) ordered_pairs = (OrderedPair*)malloc(sizeof(OrderedPair) * ordered_pairs_capacity);
    do {
        swapped = 0;
        for (int i=0;i<face_count-1;i++) {
            int f1 = order_out[i]; int f2 = order_out[i+1];
            int already_ordered = 0;
            for (int p=0;p<ordered_pairs_count;p++) { if (ordered_pairs[p].face1==f1 && ordered_pairs[p].face2==f2) { already_ordered = 1; break; } }
            if (already_ordered) continue;
            if (m->faces[f2].z_max <= m->faces[f1].z_min) continue;
            if (m->faces[f1].z_max <= m->faces[f2].z_min) {
                order_out[i] = f2; order_out[i+1] = f1; swapped = 1;
                if (ordered_pairs && ordered_pairs_count + 1 < ordered_pairs_capacity) { ordered_pairs[ordered_pairs_count].face1 = order_out[i]; ordered_pairs[ordered_pairs_count].face2 = order_out[i+1]; ordered_pairs_count++; }
                continue;
            }
            if (m->faces[f1].maxx <= m->faces[f2].minx || m->faces[f2].maxx <= m->faces[f1].minx) continue;
            if (m->faces[f1].maxy <= m->faces[f2].miny || m->faces[f2].maxy <= m->faces[f1].miny) continue;
            int res = pair_should_swap(m, f1, f2);
            if (res == 1) {
                int cand = order_out[i+1];
                for (int t = i+1; t < face_count-1; t++) order_out[t] = order_out[t+1];
                int j = i;
                while (j > 0) {
                    int prev = order_out[j-1];
                    int cmp = pair_should_swap(m, prev, cand);
                    if (cmp == 1) j--; else break;
                }
                for (int t = face_count-1; t > j; t--) order_out[t] = order_out[t-1];
                order_out[j] = cand;
                swapped = 1;
                if (ordered_pairs && ordered_pairs_count + 1 < ordered_pairs_capacity) { ordered_pairs[ordered_pairs_count].face1 = order_out[j]; ordered_pairs[ordered_pairs_count].face2 = order_out[j+1]; ordered_pairs_count++; }
            }
            if (ordered_pairs && ordered_pairs_count + 1 < ordered_pairs_capacity) { ordered_pairs[ordered_pairs_count].face1 = f2; ordered_pairs[ordered_pairs_count].face2 = f1; ordered_pairs_count++; }
        }
    } while (swapped);
    if (ordered_pairs) free(ordered_pairs);
    free(obs); g_obsv = NULL; g_model = NULL; return 1;
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

