#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fixed32.h"

// Simple OBJ loader for v x y z and f i j k ... (1-based indices)

typedef struct {
    float x,y,z;
} Vertex;

typedef struct {
    int *indices; // 0-based
    int count;
    float z_min, z_mean, z_max;
    // Cached planar coefficients and bbox to match GS3Dp
    float plane_a, plane_b, plane_c, plane_d;
    int minx, maxx, miny, maxy;
    int display_flag; // 1 = visible, 0 = hidden
} Face;

typedef struct {
    Vertex *verts; int vert_count;
    Face *faces; int face_count;
} Model;

static int ensure_verts(Model* m, int n) {
    if (m->verts == NULL) { m->verts = malloc(sizeof(Vertex)*n); m->vert_count = 0; return 1; }
    return 1;
}

Model* load_obj(const char* path) {
    // entry log
    const char* ttmp = getenv("TEMP"); char debugfn[1024]; if (ttmp) snprintf(debugfn, sizeof(debugfn), "%s\\loader_thread_debug.log", ttmp); else snprintf(debugfn, sizeof(debugfn), "loader_thread_debug.log"); FILE* dbg = fopen(debugfn, "a"); if (dbg) { fprintf(dbg, "load_obj: start path=%s\n", path); fflush(dbg); fclose(dbg); }

    FILE* f = fopen(path, "r");
    if (!f) {
        FILE* dbg2 = fopen(debugfn, "a"); if (dbg2) { fprintf(dbg2, "load_obj: fopen failed for %s (errno=%d)\n", path, errno); fclose(dbg2); }
        return NULL;
    }
    Model* m = calloc(1, sizeof(Model));
    char line[4096];
    int verts_cap = 0, faces_cap = 0;
    int vlines=0, flines=0;
    int lineno = 0;
    while (fgets(line, sizeof(line), f)) {
        lineno++;
        // debug: log the raw line prefix to help diagnose parsing
        { char preview[128]; int L = (int)strlen(line); int plen = (L<64)?L:64; memcpy(preview, line, plen); preview[plen] = '\0'; FILE* dbgline = fopen(debugfn, "a"); if (dbgline) { fprintf(dbgline, "load_obj: line %d: '%s'\n", lineno, preview); // also print hex of first bytes
            fprintf(dbgline, "load_obj: line %d hex:", lineno); for (int hi=0; hi<plen && hi<32; ++hi) fprintf(dbgline, " %02X", (unsigned char)preview[hi]); fprintf(dbgline, "\n"); fclose(dbgline); } }
        if (line[0] == 'v' && line[1] == ' ') {
            vlines++;
            float x,y,z; if (sscanf(line+2, "%f %f %f", &x, &y, &z) == 3) {
                if (m->vert_count+1 > verts_cap) { verts_cap = (verts_cap==0)?256:verts_cap*2; m->verts = realloc(m->verts, sizeof(Vertex)*verts_cap); }
                // Match GS3Dp: OBJ Z is up => swap Y/Z into internal coords
                // GS3Dp stores: x, y := OBJ z, z := OBJ y
                m->verts[m->vert_count].x = x; m->verts[m->vert_count].y = z; m->verts[m->vert_count].z = y; m->vert_count++;
                // debug log per-vertex (note swap)
                FILE* dbgv = fopen(debugfn, "a"); if (dbgv) { fprintf(dbgv, "load_obj: parsed v[%d]=%.6f %.6f %.6f (swapped Y<->Z)\n", m->vert_count-1, x, z, y); fclose(dbgv); }
            }
        }
        else {
            // more robust face detection: skip leading whitespace and check for 'f' token
            char *s = line; while (*s == ' ' || *s == '\t') s++;
            if (*s == 'f' && (s[1] == ' ' || s[1] == '\t')) {
                flines++;
                // face indices
                int tmp[64]; int c=0; char *p = s+1; while (*p) {
                    while (*p==' '||*p=='\t' || *p=='\r' || *p=='\n') p++; if (!*p) break;
                    int idx=0; int rv = sscanf(p, "%d", &idx); if (rv == 1) { tmp[c++] = idx-1; FILE* dbgtok = fopen(debugfn, "a"); if (dbgtok) { fprintf(dbgtok, "load_obj: token parsed idx=%d (token start='%c')\n", idx, *p); fclose(dbgtok); } }
                    while (*p && *p!=' ' && *p!='\t' && *p!='\n') p++; if (c>=60) break;
                }
                if (c>=3) {
                    if (m->face_count+1 > faces_cap) { faces_cap = (faces_cap==0)?256:faces_cap*2; m->faces = realloc(m->faces, sizeof(Face)*faces_cap); }
                    Face* face = &m->faces[m->face_count++]; face->count = c; face->indices = malloc(sizeof(int)*c); memcpy(face->indices, tmp, sizeof(int)*c);
                    // debug log per-face
                    FILE* dbgf = fopen(debugfn, "a"); if (dbgf) { fprintf(dbgf, "load_obj: parsed f[%d] count=%d\n", m->face_count-1, c); fclose(dbgf); }
                }
            }
        }
        if ((lineno & 1023) == 0) { FILE* dbg3 = fopen(debugfn, "a"); if (dbg3) { fprintf(dbg3, "load_obj: parsed %d lines, v=%d f=%d\n", lineno, vlines, flines); fclose(dbg3); } }
    }
    fclose(f);
    // log that file reading completed
    {
        FILE* dbg4 = fopen(debugfn, "a"); if (dbg4) { fprintf(dbg4, "load_obj: finished reading file; vlines=%d flines=%d\n", vlines, flines); fclose(dbg4); }
    }
    // compute basic z stats per face
    for (int i=0;i<m->face_count;i++) {
        Face* face = &m->faces[i]; float zmin = 1e30f, zmax = -1e30f, sum = 0; for (int k=0;k<face->count;k++) { float z = m->verts[face->indices[k]].z; if (z<zmin) zmin=z; if (z>zmax) zmax=z; sum+=z; } face->z_min = zmin; face->z_max = zmax; face->z_mean = sum / face->count; }
    // log counts
    {
        const char* tmp = getenv("TEMP"); char logfn[1024]; if (tmp) snprintf(logfn, sizeof(logfn), "%s\\viewer_obj.log", tmp); else snprintf(logfn, sizeof(logfn), "viewer_obj.log"); FILE* lf = fopen(logfn, "a"); if (lf) { fprintf(lf, "OBJ: parsed vlines=%d flines=%d vert_count=%d face_count=%d\n", vlines, flines, m->vert_count, m->face_count); fclose(lf); }
    }
    // exit log
    FILE* dbgfin = fopen(debugfn, "a"); if (dbgfin) { fprintf(dbgfin, "load_obj: finished v=%d f=%d\n", m->vert_count, m->face_count); fclose(dbgfin); }
    return m;
}

void free_model(Model* m) {
    if (!m) return; for (int i=0;i<m->face_count;i++) free(m->faces[i].indices); free(m->faces); free(m->verts); free(m);
}

