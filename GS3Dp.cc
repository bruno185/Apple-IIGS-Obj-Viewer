/*
 * ============================================================================
 *                              GS3Df.CC - Fixed32 Optimized Version
 * ============================================================================
 * 
 * 3D rendering program for Apple IIGS with ORCA/C - FIXED POINT VERSION
 * Version 0.4 Fixed-Point 64-bit - Ultra-Optimized Fixed32 Arithmetic
 * 
 * DESCRIPTION:
 *   This program reads 3D model files in simplified OBJ format,
 *   applies 3D geometric transformations (rotation, translation),
 *   projects the points on a 2D screen, and draws the resulting polygons
 *   using QuickDraw.
 *   
 *   This is the OPTIMIZED VERSION using Fixed32 16.16 arithmetic with
 *   lookup tables and 64-bit overflow protection. 
 *   Significantly faster performance compared to SANE Extended reference version.
 * 
 * FEATURES:
 *   - Reading OBJ files (vertices "v" and faces "f")
 *   - 3D transformations using Fixed32 16.16 arithmetic
 *   - 64-bit arithmetic for overflow-safe multiplication/division
 *   - Perspective projection on 2D screen
 *   - Graphic rendering with QuickDraw
 *   - Interactive interface with corrected parameter display
 *   - Performance measurements for comparison with SANE version
 * 
 * ARITHMETIC:
 *   - Fixed32 16.16 format (16 bits integer + 16 bits fractional)
 *   - FIXED_SCALE = 65536, optimized FIXED_MUL_64/FIXED_DIV_64
 *   - Pre-computed trigonometric lookup tables
 *   - Overflow-safe 64-bit intermediate calculations
 * 
 * PERFORMANCE OPTIMIZATIONS:
 *   - Direct lookup table access: rad = deg_to_rad_table[degree]
 *   - Combined transformation+projection pipeline
 *   - Pre-calculated trigonometric products
 *   - Eliminated function calls in critical loops
 * 
 * AUTHOR: Bruno
 * DATE: 2025
 * PLATFORM: Apple IIGS - ORCA/C with Fixed32 arithmetic
 * ============================================================================
 */

// ============================================================================
//                           HEADER INCLUDES
// ============================================================================

#include <stdio.h>      // Standard input/output (printf, fgets, etc.)
#include <asm.h>        // ORCA specific assembler functions
#include <string.h>     // String manipulation (strlen, strcmp, etc.)
#include <misctool.h>   // ORCA misc tools
#include <stdlib.h>     // Standard functions (malloc, free, atof, etc.)
#include <math.h>       // Math functions (cos, sin, sqrt, etc.)
#include <quickdraw.h>  // Apple IIGS QuickDraw graphics API
#include <event.h>      // System event management
#include <memory.h>     // Advanced memory management (NewHandle, etc.)
#include <window.h>     // Window management
#include <orca.h>       // ORCA specific functions (startgraph, etc.)


#pragma memorymodel 1

segment "data";

// ============================================================================
// --- Variable globale pour la vérification des indices de sommet dans readFaces ---
// (Ajouté pour garantir la cohérence OBJ)
int readVertices_last_count = 0;

// --- Global persistent polygon handle for drawing ---
// Allocated once, HLocked only during use, prevents repeated NewHandle/DisposeHandle
static Handle globalPolyHandle = NULL;
static int poly_handle_locked = 0;  // Track lock state
static int framePolyOnly = 0; // Toggle: 1 = frame-only, 0 = fill+frame (default: filled polygons)
static int painterFastMode = 1; // Toggle: 1 = fast painter (tests 1-3 only) (default: ON)

// ============================================================================
//                            FIXED POINT DEFINITIONS
// ============================================================================

/**
 * FIXED POINT ARITHMETIC - 64-bit safe version
 * =============================================
 * 
 * Cette version utilise l'arithmétique en virgule fixe 16.16 avec
 * des calculs intermédiaires 64-bit pour éviter les débordements.
 * 
 * Format: 16.16 (16 bits entiers, 16 bits fractionnaires)
 * Plage: -32768.0 à +32767.99998 avec précision de 1/65536
 */

// Basic fixed-point definitions
typedef long Fixed32;           // 32-bit fixed point number (16.16)
typedef long long Fixed64;      // 64-bit for intermediate calculations

#define FIXED_SHIFT     16                    // Number of fractional bits
#define FIXED_SCALE     (1L << FIXED_SHIFT)  // 65536
#define FIXED_MASK      (FIXED_SCALE - 1)    // 0xFFFF
#define FIXED_HALF      (FIXED_SCALE >> 1)   // 32768 (for rounding)

// Mathematical constants in fixed point
#define FIXED_PI        205887L               // PI ≈ 3.14159265 in 16.16
#define FIXED_2PI       411775L               // 2*PI ≈ 6.28318530 in 16.16
#define FIXED_PI_2      102944L               // PI/2 ≈ 1.57079632 in 16.16
#define FIXED_ONE       FIXED_SCALE           // 1.0 in 16.16
#define FIXED_PI_180    1143LL                // PI/180 ≈ 0.017453293 in 16.16 (64-bit)

// Lookup tables: sin(0..359) and cos(0..359) in Fixed32 16.16
static const Fixed32 sin_table[360] = {
    0, 1144, 2287, 3430, 4572, 5712, 6850, 7987, 9121, 10252,
    11380, 12505, 13626, 14742, 15855, 16962, 18064, 19161, 20252, 21336,
    22415, 23486, 24550, 25607, 26656, 27697, 28729, 29753, 30767, 31772,
    32768, 33754, 34729, 35693, 36647, 37590, 38521, 39441, 40348, 41243,
    42126, 42995, 43852, 44695, 45525, 46341, 47143, 47930, 48703, 49461,
    50203, 50931, 51643, 52339, 53020, 53684, 54332, 54963, 55578, 56175,
    56756, 57319, 57865, 58393, 58903, 59396, 59870, 60326, 60764, 61183,
    61584, 61966, 62328, 62672, 62997, 63303, 63589, 63856, 64104, 64332,
    64540, 64729, 64898, 65048, 65177, 65287, 65376, 65446, 65496, 65526,
    65536, 65526, 65496, 65446, 65376, 65287, 65177, 65048, 64898, 64729,
    64540, 64332, 64104, 63856, 63589, 63303, 62997, 62672, 62328, 61966,
    61584, 61183, 60764, 60326, 59870, 59396, 58903, 58393, 57865, 57319,
    56756, 56175, 55578, 54963, 54332, 53684, 53020, 52339, 51643, 50931,
    50203, 49461, 48703, 47930, 47143, 46341, 45525, 44695, 43852, 42995,
    42126, 41243, 40348, 39441, 38521, 37590, 36647, 35693, 34729, 33754,
    32768, 31772, 30767, 29753, 28729, 27697, 26656, 25607, 24550, 23486,
    22415, 21336, 20252, 19161, 18064, 16962, 15855, 14742, 13626, 12505,
    11380, 10252, 9121, 7987, 6850, 5712, 4572, 3430, 2287, 1144,
    0, -1144, -2287, -3430, -4572, -5712, -6850, -7987, -9121, -10252,
    -11380, -12505, -13626, -14742, -15855, -16962, -18064, -19161, -20252, -21336,
    -22415, -23486, -24550, -25607, -26656, -27697, -28729, -29753, -30767, -31772,
    -32768, -33754, -34729, -35693, -36647, -37590, -38521, -39441, -40348, -41243,
    -42126, -42995, -43852, -44695, -45525, -46341, -47143, -47930, -48703, -49461,
    -50203, -50931, -51643, -52339, -53020, -53684, -54332, -54963, -55578, -56175,
    -56756, -57319, -57865, -58393, -58903, -59396, -59870, -60326, -60764, -61183,
    -61584, -61966, -62328, -62672, -62997, -63303, -63589, -63856, -64104, -64332,
    -64540, -64729, -64898, -65048, -65177, -65287, -65376, -65446, -65496, -65526,
    -65536, -65526, -65496, -65446, -65376, -65287, -65177, -65048, -64898, -64729,
    -64540, -64332, -64104, -63856, -63589, -63303, -62997, -62672, -62328, -61966,
    -61584, -61183, -60764, -60326, -59870, -59396, -58903, -58393, -57865, -57319,
    -56756, -56175, -55578, -54963, -54332, -53684, -53020, -52339, -51643, -50931,
    -50203, -49461, -48703, -47930, -47143, -46341, -45525, -44695, -43852, -42995,
    -42126, -41243, -40348, -39441, -38521, -37590, -36647, -35693, -34729, -33754,
    -32768, -31772, -30767, -29753, -28729, -27697, -26656, -25607, -24550, -23486,
    -22415, -21336, -20252, -19161, -18064, -16962, -15855, -14742, -13626, -12505,
    -11380, -10252, -9121, -7987, -6850, -5712, -4572, -3430, -2287, -1144,
};
static const Fixed32 cos_table[360] = {
    65536, 65526, 65496, 65446, 65376, 65287, 65177, 65048, 64898, 64729,
    64540, 64332, 64104, 63856, 63589, 63303, 62997, 62672, 62328, 61966,
    61584, 61183, 60764, 60326, 59870, 59396, 58903, 58393, 57865, 57319,
    56756, 56175, 55578, 54963, 54332, 53684, 53020, 52339, 51643, 50931,
    50203, 49461, 48703, 47930, 47143, 46341, 45525, 44695, 43852, 42995,
    42126, 41243, 40348, 39441, 38521, 37590, 36647, 35693, 34729, 33754,
    32768, 31772, 30767, 29753, 28729, 27697, 26656, 25607, 24550, 23486,
    22415, 21336, 20252, 19161, 18064, 16962, 15855, 14742, 13626, 12505,
    11380, 10252, 9121, 7987, 6850, 5712, 4572, 3430, 2287, 1144,
    0, -1144, -2287, -3430, -4572, -5712, -6850, -7987, -9121, -10252,
    -11380, -12505, -13626, -14742, -15855, -16962, -18064, -19161, -20252, -21336,
    -22415, -23486, -24550, -25607, -26656, -27697, -28729, -29753, -30767, -31772,
    -32768, -33754, -34729, -35693, -36647, -37590, -38521, -39441, -40348, -41243,
    -42126, -42995, -43852, -44695, -45525, -46341, -47143, -47930, -48703, -49461,
    -50203, -50931, -51643, -52339, -53020, -53684, -54332, -54963, -55578, -56175,
    -56756, -57319, -57865, -58393, -58903, -59396, -59870, -60326, -60764, -61183,
    -61584, -61966, -62328, -62672, -62997, -63303, -63589, -63856, -64104, -64332,
    -64540, -64729, -64898, -65048, -65177, -65287, -65376, -65446, -65496, -65526,
    -65536, -65526, -65496, -65446, -65376, -65287, -65177, -65048, -64898, -64729,
    -64540, -64332, -64104, -63856, -63589, -63303, -62997, -62672, -62328, -61966,
    -61584, -61183, -60764, -60326, -59870, -59396, -58903, -58393, -57865, -57319,
    -56756, -56175, -55578, -54963, -54332, -53684, -53020, -52339, -51643, -50931,
    -50203, -49461, -48703, -47930, -47143, -46341, -45525, -44695, -43852, -42995,
    -42126, -41243, -40348, -39441, -38521, -37590, -36647, -35693, -34729, -33754,
    -32768, -31772, -30767, -29753, -28729, -27697, -26656, -25607, -24550, -23486,
    -22415, -21336, -20252, -19161, -18064, -16962, -15855, -14742, -13626, -12505,
    -11380, -10252, -9121, -7987, -6850, -5712, -4572, -3430, -2287, -1144,
    0, 1144, 2287, 3430, 4572, 5712, 6850, 7987, 9121, 10252,
    11380, 12505, 13626, 14742, 15855, 16962, 18064, 19161, 20252, 21336,
    22415, 23486, 24550, 25607, 26656, 27697, 28729, 29753, 30767, 31772,
    32768, 33754, 34729, 35693, 36647, 37590, 38521, 39441, 40348, 41243,
    42126, 42995, 43852, 44695, 45525, 46341, 47143, 47930, 48703, 49461,
    50203, 50931, 51643, 52339, 53020, 53684, 54332, 54963, 55578, 56175,
    56756, 57319, 57865, 58393, 58903, 59396, 59870, 60326, 60764, 61183,
    61584, 61966, 62328, 62672, 62997, 63303, 63589, 63856, 64104, 64332,
    64540, 64729, 64898, 65048, 65177, 65287, 65376, 65446, 65496, 65526,
};

// Integer-degree sin/cos helpers (direct table lookup)
static inline Fixed32 sin_deg_int(int deg) {
    deg %= 360;
    if (deg < 0) deg += 360;
    return sin_table[deg];
}
static inline Fixed32 cos_deg_int(int deg) {
    deg %= 360;
    if (deg < 0) deg += 360;
    return cos_table[deg];
}


// Conversion macros
#define INT_TO_FIXED(x)     ((Fixed32)(x) << FIXED_SHIFT)
#define FIXED_TO_INT(x)     ((int)((x) >> FIXED_SHIFT))
// Round fixed to int (nearest) handling negatives
static inline int FIXED_ROUND_TO_INT(Fixed32 x) {
    if (x >= 0) return (int)(((x) + FIXED_HALF) >> FIXED_SHIFT);
    else return (int)(((x) - FIXED_HALF) >> FIXED_SHIFT);
}
#define FLOAT_TO_FIXED(x)   ((Fixed32)((x) * FIXED_SCALE))
#define FIXED_TO_FLOAT(x)   ((float)(x) / (float)FIXED_SCALE)

// Arithmetic operations
#define FIXED_ADD(a, b)     ((a) + (b))
#define FIXED_SUB(a, b)     ((a) - (b))
#define FIXED_NEG(x)        (-(x))
#define FIXED_ABS(x)        ((x) >= 0 ? (x) : -(x))
#define FIXED_FRAC(x)       ((x) & FIXED_MASK)

// Simple multiplication and division for ORCA/C
#define FIXED_MUL(a, b)     (((long)(a) * (long)(b)) >> FIXED_SHIFT)
#define FIXED_DIV(a, b)     (((long)(a) << FIXED_SHIFT) / (long)(b))

// 64-bit safe multiplication and division for critical calculations
#define FIXED_MUL_64(a, b)  ((Fixed32)(((Fixed64)(a) * (Fixed64)(b)) >> FIXED_SHIFT))
#define FIXED_DIV_64(a, b)  ((Fixed32)(((Fixed64)(a) << FIXED_SHIFT) / (Fixed64)(b)))
#define FIXED64_TO_32(x)    ((Fixed32)(x))

// NOTE: Fixed32 angle normalization removed — code now uses integer-degree
// normalization via `normalize_deg(int)`. This avoids accidental usage of
// Fixed32-based normalization when angles are stored as integer degrees.

// Integer degree normalization
static inline int normalize_deg(int deg) {
    deg %= 360;
    if (deg < 0) deg += 360;
    return deg;
}

// ============================================================================
//                            GLOBAL CONSTANTS
// ============================================================================

// Performance and debug configuration
#define ENABLE_DEBUG_SAVE 0     // 1 = Enable debug save (SLOW!), 0 = Disable
//#define PERFORMANCE_MODE 0      // 1 = Optimized performance mode, 0 = Debug mode
// OPTIMIZATION: Performance mode - disable printf
#define PERFORMANCE_MODE 1      // 1 = no printf, 0 = normal printf

#define MAX_LINE_LENGTH 256     // Maximum file line size
#define MAX_VERTICES 6000       // Maximum vertices in a 3D model
#define MAX_FACES 6000          // Maximum faces in a 3D model (using parallel arrays)
#define MAX_FACE_VERTICES 6     // Maximum vertices per face (triangles/quads/hexagons)
#define CENTRE_X 160            // Screen center in X (320/2)
#define CENTRE_Y 100            // Screen center in Y (200/2)
//#define mode 640               // Graphics mode 640x200 pixels
#define mode 320               // Graphics mode 320x200 pixels

// ============================================================================
//                          DATA STRUCTURES
// ============================================================================

/**
 * Structure Vertex3D
 * 
 * DESCRIPTION:
 *   Represents a point in 3D space with its different representations
 *   throughout the 3D rendering pipeline.
 * 
 * FIELDS:
 *   x, y, z    : Original coordinates read from OBJ file
 *   xo, yo, zo : Transformed coordinates in the observer system
 *                (after applying rotations and translation)
 *   x2d, y2d   : Final projected coordinates on 2D screen
 * 
 * USAGE:
 *   This structure preserves all transformation steps to
 *   allow debugging and recalculations without rereading the file.
 */


// Parallel arrays for vertex data (to break 32K/64K struct limit)
typedef struct {
    Handle xHandle, yHandle, zHandle;
    Handle xoHandle, yoHandle, zoHandle;
    Handle x2dHandle, y2dHandle;
    Fixed32 *x, *y, *z;
    Fixed32 *xo, *yo, *zo;
    int *x2d, *y2d;
    int vertex_count;
} VertexArrays3D;

/**
 * Structure FaceArrays3D - Compact dynamic face storage with depth-sorted rendering
 * Each face stores ONLY the vertices it needs:
 * - vertex_count: How many vertices this face has (3 for tri, 4 for quad, etc.)
 * - vertex_indices_buffer: ONE packed buffer with all indices (NO WASTED SLOTS!)
 * - vertex_indices_ptr: Offset array pointing to each face's slice in the buffer
 * - sorted_face_indices: Array of face indices SORTED by depth (for painter's algorithm)
 * - z_max: Depth for sorting
 * - display_flag: Culling flag
 * 
 * MEMORY LAYOUT:
 * Instead of 4 arrays of 6000 elements each, we use ONE packed buffer.
 * Triangles (1538 faces × 3 indices) + Quads (2504 faces × 4 indices) = packed linearly
 * Saves ~40-60% memory vs fixed 4 vertices/face
 * 
 * DEPTH SORTING STRATEGY:
 * Instead of moving data around (complex with variable-length indices), we maintain
 * sorted_face_indices[] which contains face numbers in depth order (farthest first).
 * Drawing loop: for(i=0; i<face_count; i++) { int face_id = sorted_face_indices[i]; ... }
 * This keeps the buffer untouched while providing correct rendering order.
 */
typedef struct {
    Handle vertex_countHandle;           // 1 array: face_count × 4 bytes
    Handle vertex_indicesBufferHandle;   // 1 buffer: all indices packed (NO WASTED SLOTS!)
    Handle vertex_indicesPtrHandle;      // 1 array: offset to each face's indices
    Handle z_maxHandle;                  // 1 array: face_count × 4 bytes
    Handle display_flagHandle;           // 1 array: face_count × 4 bytes
    Handle sorted_face_indicesHandle;    // 1 array: face numbers sorted by depth
    
    int *vertex_count;                   // Points to: [3, 3, 4, 3, 4, ...]
    int *vertex_indices_buffer;          // Points to: [v1, v2, v3, v1, v2, v3, v4, v1, v2, v3, ...]
    int *vertex_indices_ptr;             // Points to: [offset0, offset3, offset6, offset10, ...]
    Fixed32 *z_min;
    Fixed32 *z_max;
    Fixed32 *z_mean;                     // Mean zo per face (computed in calculateFaceDepths)
    Fixed32 *plane_a;                     // per-face normalized normal X (a)
    Fixed32 *plane_b;                     // per-face normalized normal Y (b)
    Fixed32 *plane_c;                     // per-face normalized normal Z (c)
    Fixed32 *plane_d;                     // per-face D term for plane equation
    int *minx;                            // cached 2D bounding box (projected x/y)
    int *maxx;
    int *miny;
    int *maxy;
    int *display_flag;
    int *sorted_face_indices;            // Points to: [face_id1, face_id2, ...] sorted by z_max
    int face_count;                      // Actual number of loaded faces
    int total_indices;                   // Total indices across all faces (sum of all vertex_counts)
} FaceArrays3D;

/**
 * Structure Face3D
 * 
 * DESCRIPTION:
 *   Represents a face (polygon) of a 3D object. A face is defined
 *   by a list of indices pointing to vertices in the model's
 *   vertex array.
 * 
 * FIELDS:
 *   vertex_count    : Number of vertices composing this face (3+ for polygon)
 *   vertex_indices  : Array of vertex indices (1-based numbering as in OBJ format)
 * 
 * NOTES:
 *   - Indices are stored in base 1 (first vertex = index 1)
 *   - Conversion to base 0 needed to access the C array
 *   - Maximum MAX_FACE_VERTICES vertices per face (now 6 for triangles/quads/hexagons)
 *   - LEGACY STRUCTURE: Now replaced by FaceArrays3D for parallel array storage
 */
typedef struct {
    int vertex_count;                           // Number of vertices in the face
    int vertex_indices[MAX_FACE_VERTICES];     // Vertex indices (base 1, max 6 for polygons)
    Fixed32 z_max;                             // Maximum depth of the face (for sorting, Fixed Point)
    int display_flag;                          // 1 = display face, 0 = don't display (behind camera)
} Face3D;

/**
 * Structure DynamicPolygon
 * 
 * DESCRIPTION:
 *   Structure compatible with QuickDraw for drawing polygons.
 *   This structure must be dynamically allocated because its size
 *   varies according to the number of points in the polygon.
 * 
 * FIELDS:
 *   polySize    : Total size of the structure in bytes
 *   polyBBox    : Polygon bounding box rectangle
 *   polyPoints  : Array of polygon points in screen coordinates
 * 
 * QUICKDRAW FORMAT:
 *   QuickDraw expects a structure with header (size + bbox) followed
 *   by points. The size must include the header + all points.
 */
typedef struct {
    int polySize;                               // Total structure size (bytes)
    Rect polyBBox;                             // Bounding box rectangle
    Point polyPoints[MAX_FACE_VERTICES];       // Polygon points (screen coordinates)
} DynamicPolygon;

/**
 * Structure ObserverParams
 * 
 * DESCRIPTION:
 *   Contains all parameters defining the position and orientation
 *   of the observer (camera) in 3D space, as well as projection
 *   parameters.
 * 
 * FIELDS:
 *   angle_h  : Horizontal rotation angle of the observer (degrees)
 *              Rotation around Y-axis (left/right)
 *   angle_v  : Vertical rotation angle of the observer (degrees)
 *              Rotation around X-axis (up/down)
 *   angle_w  : Screen projection rotation angle (degrees)
 *              Rotation in the final 2D plane
 *   distance : Distance from observer to model center
 *              Larger = smaller object, smaller = larger object
 * 
 * MATHEMATICAL NOTES:
 *   - Angles are in degrees (converted to radians for calculations)
 *   - Distance affects perspective and apparent size
 *   - angle_w allows final rotation to adjust orientation
 */
typedef struct {
    int angle_h;   // Observer horizontal angle (degrees)
    int angle_v;   // Observer vertical angle (degrees)
    int angle_w;   // 2D projection rotation angle (degrees)
    Fixed32 distance;  // Observer-object distance (perspective, Fixed Point)
} ObserverParams;

/**
 * Structure Model3D
 * 
 * DESCRIPTION:
 *   Main structure containing all data of a 3D model.
 *   It groups vertices, faces, and associated counters.
 * 
 * FIELDS:
 *   vertices      : Pointer to dynamic vertex array
 *   faces         : Pointer to dynamic face array
 *   vertex_count  : Actual number of loaded vertices
 *   face_count    : Actual number of loaded faces
 * 
 * MEMORY MANAGEMENT:
 *   - Arrays are dynamically allocated (malloc)
 *   - Allows exceeding Apple IIGS stack limits
 *   - Mandatory cleanup with destroyModel3D()
 * 
 * USAGE:
 *   Model3D* model = createModel3D();
 *   loadModel3D(model, "file.obj");
 *   // ... usage ...
 *   destroyModel3D(model);
 */
typedef struct {
    VertexArrays3D vertices;          // Parallel arrays for all vertex data
    FaceArrays3D faces;               // Parallel arrays for all face data

    /* Auto-scale metadata (non-destructive) */
    Fixed32 auto_scale;               // Fixed32 scale factor applied (FIXED_ONE if none)
    Fixed32 auto_center_x, auto_center_y, auto_center_z; // center used during scaling
    int auto_scaled;                  // 1 if auto-scaling has been applied
    int auto_centered;                // 1 if auto-scale used centering

    /* Optional backup of original coordinates to allow exact revert and dynamic scale adjustments */
    Fixed32 *orig_x;                   // NULL if no backup
    Fixed32 *orig_y;
    Fixed32 *orig_z;
    float *radius_buf;            // scratch buffer for radii (squared distances) (float for cache)
    int radius_buf_capacity;        // capacity of radius_buf
    float *coord_buf;               // scratch buffer for converted float coordinates (x,y,z interleaved)
    int coord_buf_capacity;         // capacity in number of vertices for coord_buf

    /* Bounding sphere (computed once at load) */
    float bs_cx;   // center x (float)
    float bs_cy;   // center y
    float bs_cz;   // center z
    float bs_r;    // radius
    int bs_valid;  // 0 = not computed, 1 = valid
} Model3D;

// ============================================================================
//                       FUNCTION DECLARATIONS
// ============================================================================

/**
 * FIXED POINT MATHEMATICAL FUNCTIONS
 * ===================================
 */
/**
 * OBJ FILE READING FUNCTIONS
 * ===========================
/**
 * readVertices
 * 
 * DESCRIPTION:
 *   Reads vertices (3D points) from an OBJ format file.
 *   Searches for lines starting with "v " and extracts X,Y,Z coordinates.
 * 
 * PARAMETERS:
 *   filename     : OBJ filename to read
 *   vertices     : Destination array to store vertices
 *   max_vertices : Maximum array size (overflow protection)
 * 
 * RETURN:
 *   Number of successfully read vertices, or -1 on error
 * 
 * OBJ FORMAT:
 *   v 1.234 5.678 9.012
 *   v -2.5 0.0 3.14
 */
int readVertices(const char* filename, VertexArrays3D* vtx, int max_vertices);

/**
 * readFaces
 * 
 * DESCRIPTION:
 *   Reads faces (polygons) from an OBJ format file.
 *   Searches for lines starting with "f " and extracts vertex indices.
 * 
 * PARAMETERS:
 *   filename   : OBJ filename to read
 *   faces      : Destination array to store faces
 *   max_faces  : Maximum array size (overflow protection)
 * 
 * RETURN:
 *   Number of successfully read faces, or -1 on error
 * 
 * OBJ FORMAT:
 *   f 1 2 3        (triangle with vertices 1, 2, 3)
 *   f 4 5 6 7      (quadrilateral with vertices 4, 5, 6, 7)
 */

/**
 * 3D GEOMETRIC TRANSFORMATION FUNCTIONS
 * =====================================
 */

/**
 * projectTo2D
 * 
 * DESCRIPTION:
 *   Projects the transformed 3D coordinates onto a 2D screen using
 *   perspective projection. Also applies a final rotation
 *   in the 2D plane.
 * 
 * PARAMETERS:
 *   vertices     : Array of vertices to project
 *   vertex_count : Number of vertices in the array

/**
 * projectTo2D
 * 
 * DESCRIPTION:
 *   Projects the transformed 3D coordinates onto a 2D screen using
 *   perspective projection. Also applies a final rotation
 *   in the 2D plane.
 * 
 * PARAMETERS:
 *   vertices     : Array of vertices to project
 *   vertex_count : Number of vertices in the array
 *   angle_w      : Rotation angle in the 2D plane (degrees)
 * 
 * ALGORITHM:
 *   1. Perspective projection: x2d = (xo * scale) / zo + center_x
 *   2. Same for y2d with Y-axis inversion
 *   3. Final rotation in the 2D plane according to angle_w
 *   4. Points behind observer (zo <= 0) marked invisible
 * 
 * RESULTING COORDINATES:
 *   The x2d, y2d fields of the vertices contain the final screen coordinates.
 */
void projectTo2D(VertexArrays3D* vtx, int angle_w_deg);

/* Bounding sphere helpers */
void computeModelBoundingSphere(Model3D* model);
Fixed32 computeDistanceFromBoundingSphere(Model3D* model, float margin);

/* DEPRECATED: computeDistanceToFit - vertex-based O(n) distance fit.
 * Historically used to compute observer distance by scanning all vertices.
 * Replaced by:
 *   - computeModelBoundingSphere() (computed once at model load)
 *   - computeDistanceFromBoundingSphere() (fast O(1) estimate based on sphere)
 * Kept here as a fallback and for historical reference; consider removal after validation.
 */
Fixed32 computeDistanceToFit(VertexArrays3D* vtx, float margin);
void getObserverParams(ObserverParams* params, Model3D* model);

/* Auto-scaling helpers (non-destructive): allow automatic scaling on import with rollback */
void autoScaleModel(Model3D* model, float target_max_dim, float min_scale, float max_scale, int center_flag);
void revertAutoScaleModel(Model3D* model);

// Fit model to view using sphere metric with percentile trimming (non-destructive)
void fitModelToView(Model3D* model, ObserverParams* params, float target_max_dim, float margin, float percentile, int center_flag);

// Fast distance adjust: try to approximate change in distance by scaling 2D coords and updating depths/plane_d.
// Returns 1 if applied safely, 0 if caller should fall back to full recompute.


// Internal helpers: backup/restore original vertex arrays
void backupModelCoords(Model3D* model);
void freeBackupModelCoords(Model3D* model);


/**
 * GRAPHIC RENDERING FUNCTIONS
 * ===========================
 */

/**
 * drawPolygons
 * 
 * DESCRIPTION:
 *   Draws all polygons (faces) of the 3D model on screen using
 *   Apple IIGS QuickDraw API. Each face is rendered with a different
 *   color for visualization.
 * 
 * PARAMETERS:
 *   vertices   : Array of vertices with calculated 2D coordinates
 *   faces      : Array of faces to draw
 *   face_count : Number of faces in the array
 * 
 * ALGORITHM:
 *   1. QuickDraw graphics mode initialization
 *   2. For each face:
 *      - Check vertex visibility
 *      - Create QuickDraw polygon structure
 *      - Calculate bounding box
 *      - Dynamic memory allocation
 *      - Drawing with PaintPoly()
 *      - Memory cleanup
 * 
 * COLOR MANAGEMENT:
 *   Cyclic colors based on face index (i % 15 + 1)
 * 
 * OPTIMIZATIONS:
 *   - Faces with less than 3 visible vertices ignored
 *   - Off-screen vertices handled correctly
 */
void drawPolygons(Model3D* model, int* vertex_count, int face_count, int vertex_count_total);
void calculateFaceDepths(Model3D* model, Face3D* faces, int face_count);


/**
 * PAINTER'S ALGORITHM - NEWELL/SANCHA
 * ===================================
 * Tri Z décroissant, puis test de recouvrement 2D (bounding box)
 * Corrige l'ordre si recouvrement ambigu (sans split)
 *
 * Appel : painter_newell_sancha(model, face_count);
 */

// Comparator support for qsort in painter_newell_sancha
static FaceArrays3D* qsort_faces_ptr_for_cmp = NULL;
static int cmp_faces_by_zmean(const void* pa, const void* pb) {
    int a = *(const int*)pa;
    int b = *(const int*)pb;
    Fixed32 za = qsort_faces_ptr_for_cmp->z_mean[a];
    Fixed32 zb = qsort_faces_ptr_for_cmp->z_mean[b];
    if (za > zb) return -1;   // larger z_mean first (descending)
    if (za < zb) return 1;
    if (a < b) return -1;     // tie-breaker: smaller index first
    if (a > b) return 1;
    return 0;
}

/**
 * FAST VERSION: Only performs Test 1 (depth overlap), Test 2 (X bbox), Test 3 (Y bbox)
 * No plane coefficients, no pair caching. Much faster but less robust.
 */
void painter_newell_sancha_fast(Model3D* model, int face_count) {
    FaceArrays3D* faces = &model->faces;
    if (!faces->z_mean) return; // safety

    // Initial stable sort by z_mean (descending) with tie-breaker on index
    qsort_faces_ptr_for_cmp = faces;
    qsort(faces->sorted_face_indices, face_count, sizeof(int), cmp_faces_by_zmean);
    qsort_faces_ptr_for_cmp = NULL;

/*    2. Correction stricte d'ordre avec les trois tests de plan (Newell/Sancha)
    int swapped;
    int pass = 0;
    const int max_passes = face_count * 2; // safety cap
    do {
        swapped = 0;
        int i;
        for (i = 0; i < face_count - 1; i++) {
            int f1 = faces->sorted_face_indices[i];
            int f2 = faces->sorted_face_indices[i+1];

            // Test 1: Depth overlap
            if (faces->z_max[f2] <= faces->z_min[f1]) continue;
            if (faces->z_max[f1] < faces->z_min[f2]) {
                // Swap needed
                faces->sorted_face_indices[i] = f2;
                faces->sorted_face_indices[i+1] = f1;
                swapped = 1;
                continue;
            }

            // Test 2: X overlap
            int minx1 = faces->minx[f1], maxx1 = faces->maxx[f1];
            int minx2 = faces->minx[f2], maxx2 = faces->maxx[f2];
            if (maxx1 <= minx2 || maxx2 <= minx1) continue;

            // Test 3: Y overlap
            int miny1 = faces->miny[f1], maxy1 = faces->maxy[f1];
            int miny2 = faces->miny[f2], maxy2 = faces->maxy[f2];
            if (maxy1 <= miny2 || maxy2 <= miny1) continue;


        }
        pass++;
    } while (swapped && pass < max_passes); */
}

void painter_newell_sancha(Model3D* model, int face_count) {
    // ...existing code...
    FaceArrays3D* faces = &model->faces;
    VertexArrays3D* vtx = &model->vertices;
    int i, j;
    Fixed32* face_zmean = faces->z_mean;
    if (!face_zmean) return; // sécurité


    long t_start = GetTick();
    // Etape 1 : Tri Z décroissant STABLE sur la moyenne, avec tie-breaker sur l'indice d'origine
    for (i = 0; i < face_count; i++) faces->sorted_face_indices[i] = i;
    // Use qsort for O(n log n) sorting while preserving the exact tie-breaker
    qsort_faces_ptr_for_cmp = faces;
    qsort(faces->sorted_face_indices, face_count, sizeof(int), cmp_faces_by_zmean);
    qsort_faces_ptr_for_cmp = NULL;
    long t_end = GetTick();
    if (!PERFORMANCE_MODE)
    {
        long elapsed = t_end - t_start;
        double ms = ((double)elapsed * 1000.0) / 60.0; // 60 ticks per second
        printf("[TIMING] initial sort (qsort): %ld ticks (%.2f ms)\n", elapsed, ms);
    }
    
    // 3. Correction stricte d'ordre avec les deux tests de plan (Newell/Sancha)
    
    int swap_count = 0;
    int swapped = 0; // flag utilisé par la boucle de correction


    // Structure pour stocker les paires de faces ordonnées
    typedef struct {
        int face1;  // Face qui doit être avant
        int face2;  // Face qui doit être après
    } OrderedPair;
    
    // Préallocation unique : meilleure performance en évitant realloc fréquents.
    // On préalloue une capacité basée sur face_count * 4 (choix empirique).
    int ordered_pairs_capacity = face_count * 4;
    OrderedPair* ordered_pairs = NULL;
    if (ordered_pairs_capacity > 0) {
        ordered_pairs = (OrderedPair*)malloc(ordered_pairs_capacity * sizeof(OrderedPair));
        if (!ordered_pairs) {
            // Si l'allocation échoue, revenir au mode dynamique par réallocation (capacity = 0)
            ordered_pairs_capacity = 0;
        }
    }
    int ordered_pairs_count = 0;

    // Oscillation / runaway protection
    int pass_number = 0;
    const int max_passes = 2000; // safety cap to avoid infinite loops
    unsigned int last_checksum = 0;
    unsigned int prev_checksum = 0;

    do {
        swapped = 0;
        int t1 = 0;
        int t2 = 0;
        int t3 = 0;
        int t4 = 0;
        int t5 = 0;
        int t6 = 0;
        int t7 = 0;

        for (i = 0; i < face_count-1; i++) {
            int f1 = faces->sorted_face_indices[i];
            int f2 = faces->sorted_face_indices[i+1];
            
            
            // Vérifier si cette paire a déjà été ordonnée définitivement
            int already_ordered = 0;
            int p;
            for (p = 0; p < ordered_pairs_count; p++) {
                if (ordered_pairs[p].face1 == f1 && ordered_pairs[p].face2 == f2) {
                    already_ordered = 1;
                    break;
                }
            }
            if (already_ordered) {
                continue;
            }

            t1++;
            // Test 1 : Depth overlap
            // Depth quick test: if farthest point of P2 is in front of nearest point of P1, order is respected
            if (faces->z_max[f2] <= faces->z_min[f1]) continue;
            if (faces->z_max[f1] <= faces->z_min[f2]) goto do_swap;
            
            // Bounding Box 2D overlap (split into X and Y parts)
            // Use cached bounding boxes (computed in calculateFaceDepths)

            t2++;
            // Test 2 : X overlap only
            int minx1 = faces->minx[f1], maxx1 = faces->maxx[f1], miny1 = faces->miny[f1], maxy1 = faces->maxy[f1];
            int minx2 = faces->minx[f2], maxx2 = faces->maxx[f2], miny2 = faces->miny[f2], maxy2 = faces->maxy[f2];

            if (maxx1 <= minx2 || maxx2 <= minx1) continue;
            
            t3++;
            // Test 3 : Y overlap only
            if (maxy1 <= miny2 || maxy2 <= miny1) continue;

            // Use cached plane normals and d terms computed in calculateFaceDepths
            int n1 = faces->vertex_count[f1];
            int n2 = faces->vertex_count[f2];
            int offset1 = faces->vertex_indices_ptr[f1];
            int offset2 = faces->vertex_indices_ptr[f2];
            int k;
            Fixed32 a1 = faces->plane_a[f1];
            Fixed32 b1 = faces->plane_b[f1];
            Fixed32 c1 = faces->plane_c[f1];
            Fixed32 d1 = faces->plane_d[f1];
            Fixed32 a2 = faces->plane_a[f2];
            Fixed32 b2 = faces->plane_b[f2];
            Fixed32 c2 = faces->plane_c[f2];
            Fixed32 d2 = faces->plane_d[f2];
            Fixed32 epsilon = FLOAT_TO_FIXED(0.000001f);

            int obs_side1 = 0; // côté de l'observateur par rapport au plan de f1 : +1, -1 ou 0 (inconclusive)
            int obs_side2 = 0; // côté de l'observateur par rapport au plan de f2 : +1, -1 ou 0 (inconclusive)
            int side;           // coté du vertex.
            int all_same_side; // flag pour indiquer si tous les vertex sont du même coté
            int all_opposite_side; // flag pour indiquer si tous les vertex sont du coté opposé
            Fixed32 test_value;

            t4++;
            // Test 4 : Test si f2 est du même côté que l'observatur par rapport au plan de f1. 
            // Si oui, f2 est bien devant f1, pas d'échange.
            if (ENABLE_DEBUG_SAVE) {
            printf("Test 4 : Testing faces %d and %d\n", f1, f2);
            }
            //keypress();
            obs_side1 = 0; // sign of d1: +1, -1 or 0 (inconclusive)
            if (d1 > epsilon) obs_side1 = 1; 
            else if (d1 < -epsilon) obs_side1 = -1;
            else goto skipT4; // si l'observateur est sur le plan, on ne peut rien conclure, il faut faire d'autres tests
            all_same_side = 1;
            for (k=0; k<n2; k++) {
                    int v = faces->vertex_indices_buffer[offset2+k]-1;
                    // test_value = a1*vtx->xo[v] + b1*vtx->yo[v] + c1*vtx->zo[v] + d1;
                    test_value = a2*vtx->xo[v] + b2*vtx->yo[v] + c2*vtx->zo[v] + d2;
                    if  (test_value > epsilon) side = 1;
                    else if (test_value < -epsilon) side = -1;
                    if (obs_side1 != side) { 
                        // si un vertex est de l'autre coté, on sort de la boucle
                        // et on met le flag à 0 pour indiquer que le test a échoué (et passer au test suivant)
                        all_same_side = 0;   
                        break; 
                    }
            }
            if (all_same_side) continue; // faces are ordered correctly, move to next pair

            skipT4:

            t5++;
            // Test 5 : Test si f1 est du coté opposé de l'observateur par rapport au plan de f2.      
            // Si oui, f1 est devant f2, pas d'échange
            if (ENABLE_DEBUG_SAVE) {
            printf("Test 5 : Testing faces %d and %d\n", f1, f2);
          }
            obs_side2 = 0; // sign of d1: +1, -1 or 0 (inconclusive)
            if (d2 > epsilon) obs_side2 = 1; 
            else if (d2 < -epsilon) obs_side2 = -1;
            else goto skipT5; // si l'observateur est sur le plan, on ne peut rien conclure, il faut faire d'autres tests
            all_opposite_side = 1;
            for (k=0; k<n1; k++) {
                int v = faces->vertex_indices_buffer[offset1+k]-1;
                // test_value = a2*vtx->xo[v] + b2*vtx->yo[v] + c2*vtx->zo[v] + d2;
                test_value = a1*vtx->xo[v] + b1*vtx->yo[v] + c1*vtx->zo[v] + d1;
                if  (test_value > epsilon) side = 1;
                else if (test_value < -epsilon) side = -1;
                if (obs_side2 == side) {
                    // si un vertex est du même coté, on sort de la boucle
                    // et on met le flag à 0 pour indiquer que le test a échoué (et passer au test suivant)
                    all_opposite_side = 0; 
                    break; }
                }
                if (all_opposite_side) continue; // faces are ordered correctly, move to next pair
            
            skipT5:

            t6++;
            // Test 6 : Test si f2 est du  côté opposé de l'observateur par rapport au plan de f1. 
            // Si oui, f2 est derrière f1, on doit échanger l'ordre
                
            if (ENABLE_DEBUG_SAVE) {
            printf("Test 6 : Testing faces %d and %d\n", f1, f2);
                }
            obs_side1 = 0; // sign of d1: +1, -1 or 0 (inconclusive)
            if (d1 > epsilon) obs_side1 = 1; 
            else if (d1 < -epsilon) obs_side1 = -1;
            else goto skipT6; // si l'observateur est sur le plan, on ne peut rien conclure, il faut faire d'autres tests

                all_opposite_side = 1;
                for (k=0; k<n2; k++) {
                    int v = faces->vertex_indices_buffer[offset2+k]-1;
                    int side;
                    test_value = a2*vtx->xo[v] + b2*vtx->yo[v] + c2*vtx->zo[v] + d2;
                    // test_value = a1*vtx->xo[v] + b1*vtx->yo[v] + c1*vtx->zo[v] + d1;
                    if  (test_value > epsilon) side = 1;
                    else side = -1;
                    if (obs_side1 == side) { 
                        all_opposite_side = 0; 
                        break; 
                        }
                }
                if (all_opposite_side == 0) continue;
                // f2 n'est pas du coté opposé de l'observateur, donc f2 n'est pas derrière f1

                // Si on arrive ici, f2 est du même côté que l'observateur, donc f2 est devant f1
                // on peut donc inverser l'ordre des faces
                else {
                    goto do_swap;
                }

            skipT6: ;

            t7++;
            // Test 7 : Test si f1 est du même côté de l'observateur par rapport au plan de f2. 
            // Si oui, f1 est devant f2, on doit échanger l'ordre

            if (ENABLE_DEBUG_SAVE) {
            printf("Test 7 : Testing faces %d and %d\n", f1, f2);
                }
            obs_side2 = 0; // sign of d1: +1, -1 or 0 (inconclusive)
            if (d2 > epsilon) obs_side2 = 1; 
            else if (d2 < -epsilon) obs_side2 = -1;
            else goto skipT7; // si l'observateur est sur le plan, on ne peut rien conclure, il faut faire d'autres tests
            all_same_side = 1;
            for (k=0; k<n1; k++) {
                int v = faces->vertex_indices_buffer[offset1+k]-1;
                int side;
                // test_value = a2*vtx->xo[v] + b2*vtx->yo[v] + c2*vtx->zo[v] + d2;
                test_value = a1*vtx->xo[v] + b1*vtx->yo[v] + c1*vtx->zo[v] + d1;
                if  (test_value > epsilon) side = 1;
                else side = -1;
                if (obs_side2 != side) { 
                    all_same_side = 0; 
                    break; 
                    }
            }
                if (all_same_side == 0) goto skipT7;
                // f1 n'est pas du même côté de l'observateur, donc f1 n'est pas devant f2
                // on ne doit pas échanger l'ordre des faces
                else {
                    goto do_swap;
                }

            do_swap: {

                if (ENABLE_DEBUG_SAVE) {
                printf("Swapping faces %d and %d\n", f1, f2);
                // removed blocking keypress() to avoid hangs in GS runtime
                }

                int tmp = faces->sorted_face_indices[i];
                faces->sorted_face_indices[i] = faces->sorted_face_indices[i+1];
                faces->sorted_face_indices[i+1] = tmp;
                swapped = 1;
                swap_count++;
                
                // Ajouter cette paire à la liste des paires ordonnées
                // Après l'échange, f2 est maintenant avant f1 dans le tableau
                // Utiliser uniquement le buffer préalloué (pas de realloc) : si on dépasse la capacité, on ignore la paire
                if (ordered_pairs != NULL && ordered_pairs_count < ordered_pairs_capacity) {
                    ordered_pairs[ordered_pairs_count].face1 = f2;
                    ordered_pairs[ordered_pairs_count].face2 = f1;
                    ordered_pairs_count++;
                }
            }

        skipT7: ;
        if (ENABLE_DEBUG_SAVE){
                printf("NON CONCLUTANT POUR LES FACES %d ET %d\n", f1, f2);
        }
        // on les met dans la liste des paires ordonnées pour ne plus les tester
        if (ordered_pairs != NULL && ordered_pairs_count < ordered_pairs_capacity) {
                ordered_pairs[ordered_pairs_count].face1 = f2;
                ordered_pairs[ordered_pairs_count].face2 = f1;
                ordered_pairs_count++;
        }
        // keypress();
        // Ici, on devrait découper f1 par f2 (ou inversement), mais on ne le fait pas pour l'instant
        }
        if (ENABLE_DEBUG_SAVE) {printf("Pass completed, swaps this pass: %d\n", swap_count);
                // removed blocking keypress();
        if (swapped) {
                printf("swapped = %d\n", swapped);
                keypress();
                }
        }

        if (ENABLE_DEBUG_SAVE) {
        printf("t1=%d t2=%d t3=%d t4=%d t5=%d t6=%d t7=%d\n", t1, t2, t3, t4, t5, t6, t7);
        keypress();
        }

    } while (swapped);
    
    // Libérer la mémoire de la liste des paires ordonnées
    if (ordered_pairs) {
        free(ordered_pairs);
    }    
}
/**
 * UTILITY FUNCTIONS
 * ==================
 */
// ============================================================================
//                    3D MODEL MANAGEMENT FUNCTIONS
// ============================================================================
/**
 * CREATING A NEW 3D MODEL
 * ========================
 * 
 * This function dynamically allocates all structures necessary
 * for a 3D model. Dynamic allocation is crucial on Apple IIGS
 * because the stack is limited and cannot contain large arrays.
 * 
 * ALLOCATION STRATEGY:
 * 1. Main Model3D structure allocation
 * 2. Vertex array allocation (MAX_VERTICES elements)
 * 3. Face array allocation (MAX_FACES elements)
 * 4. On failure: cleanup of previous allocations
 * 
 * ERROR HANDLING:
 * - Check each allocation
 * - Automatic cascade cleanup on partial failure
 * - Return NULL if unable to allocate
 */
Model3D* createModel3D(void) {
    // Step 1: Main structure allocation
    Model3D* model = (Model3D*)malloc(sizeof(Model3D));
    if (model == NULL) {
        return NULL;
    }
    int n = MAX_VERTICES;
    model->vertices.vertex_count = n;

    // Initialize auto-scale metadata
    model->auto_scale = FIXED_ONE;
    model->auto_center_x = 0;
    model->auto_center_y = 0;
    model->auto_center_z = 0;
    model->auto_scaled = 0;
    model->auto_centered = 0;
    model->orig_x = NULL;
    model->orig_y = NULL;
    model->orig_z = NULL;
    model->radius_buf = NULL;
    model->radius_buf_capacity = 0;
    model->coord_buf = NULL;
    model->coord_buf_capacity = 0;
    model->bs_cx = 0.0f;
    model->bs_cy = 0.0f;
    model->bs_cz = 0.0f;
    model->bs_r = 0.0f;
    model->bs_valid = 0;
    
    // Step 2: Allocate vertex arrays using malloc (handles bank crossing better)
    // Note: malloc() should handle bank boundaries better than NewHandle()
    model->vertices.x = (Fixed32*)malloc(n * sizeof(Fixed32));
    model->vertices.y = (Fixed32*)malloc(n * sizeof(Fixed32));
    model->vertices.z = (Fixed32*)malloc(n * sizeof(Fixed32));
    model->vertices.xo = (Fixed32*)malloc(n * sizeof(Fixed32));
    model->vertices.yo = (Fixed32*)malloc(n * sizeof(Fixed32));
    model->vertices.zo = (Fixed32*)malloc(n * sizeof(Fixed32));
    model->vertices.x2d = (int*)malloc(n * sizeof(int));
    model->vertices.y2d = (int*)malloc(n * sizeof(int));
    
    if (!model->vertices.x || !model->vertices.y || !model->vertices.z ||
        !model->vertices.xo || !model->vertices.yo || !model->vertices.zo ||
        !model->vertices.x2d || !model->vertices.y2d) {
        printf("Error: Unable to allocate memory for vertex arrays\n");
        keypress();
        // Allocation failed, cleanup
        if (model->vertices.x) free(model->vertices.x);
        if (model->vertices.y) free(model->vertices.y);
        if (model->vertices.z) free(model->vertices.z);
        if (model->vertices.xo) free(model->vertices.xo);
        if (model->vertices.yo) free(model->vertices.yo);
        if (model->vertices.zo) free(model->vertices.zo);
        if (model->vertices.x2d) free(model->vertices.x2d);
        if (model->vertices.y2d) free(model->vertices.y2d);
        free(model);
        return NULL;
    }
    
    // Set dummy handles to NULL (not used with malloc)
    model->vertices.xHandle = NULL;
    model->vertices.yHandle = NULL;
    model->vertices.zHandle = NULL;
    model->vertices.xoHandle = NULL;
    model->vertices.yoHandle = NULL;
    model->vertices.zoHandle = NULL;
    model->vertices.x2dHandle = NULL;
    model->vertices.y2dHandle = NULL;
    
    // Step 3: Face array allocation using parallel arrays (like vertices)
    // Each element stored separately to fit 32KB limit per allocation
    int nf = MAX_FACES;
    
    // Allocate vertex count array: nf * 4 bytes = 24KB
    model->faces.vertex_count = (int*)malloc(nf * sizeof(int));
    if (!model->faces.vertex_count) {
        printf("Error: Unable to allocate memory for face vertex_count array\n");
        keypress();
        free(model->vertices.x);
        free(model->vertices.y);
        free(model->vertices.z);
        free(model->vertices.xo);
        free(model->vertices.yo);
        free(model->vertices.zo);
        free(model->vertices.x2d);
        free(model->vertices.y2d);
        free(model);
        return NULL;
    }
    
    // Allocate SINGLE packed buffer for all vertex indices
    // Estimate: average 3.5 indices per face (mix of triangles and quads)
    // For 6000 faces: ~21KB. We allocate conservatively at 5 per face = 120KB max
    int estimated_total_indices = nf * 5;
    model->faces.vertex_indices_buffer = (int*)malloc(estimated_total_indices * sizeof(int));
    if (!model->faces.vertex_indices_buffer) {
        printf("Error: Unable to allocate memory for vertex_indices_buffer\n");
        keypress();
        free(model->vertices.x);
        free(model->vertices.y);
        free(model->vertices.z);
        free(model->vertices.xo);
        free(model->vertices.yo);
        free(model->vertices.zo);
        free(model->vertices.x2d);
        free(model->vertices.y2d);
        free(model->faces.vertex_count);
        free(model);
        return NULL;
    }
    
    // Allocate offset array: one offset per face into the packed buffer
    model->faces.vertex_indices_ptr = (int*)malloc(nf * sizeof(int));
    if (!model->faces.vertex_indices_ptr) {
        printf("Error: Unable to allocate memory for vertex_indices_ptr array\n");
        keypress();
        free(model->vertices.x);
        free(model->vertices.y);
        free(model->vertices.z);
        free(model->vertices.xo);
        free(model->vertices.yo);
        free(model->vertices.zo);
        free(model->vertices.x2d);
        free(model->vertices.y2d);
        free(model->faces.vertex_count);
        free(model->faces.vertex_indices_buffer);
        free(model);
        return NULL;
    }
    
    // Allocate z_min and z_max arrays (min and max depth per face)
    model->faces.z_min = (Fixed32*)malloc(nf * sizeof(Fixed32));
    model->faces.z_max = (Fixed32*)malloc(nf * sizeof(Fixed32));
    if (!model->faces.z_min || !model->faces.z_max) {
        printf("Error: Unable to allocate memory for face depth arrays\n");
        keypress();
        if (model->faces.z_min) free(model->faces.z_min);
        if (model->faces.z_max) free(model->faces.z_max);
        free(model->vertices.x);
        free(model->vertices.y);
        free(model->vertices.z);
        free(model->vertices.xo);
        free(model->vertices.yo);
        free(model->vertices.zo);
        free(model->vertices.x2d);
        free(model->vertices.y2d);
        free(model->faces.vertex_count);
        free(model->faces.vertex_indices_buffer);
        free(model->faces.vertex_indices_ptr);
        free(model);
        return NULL;
    }
    // Allocate z_mean array (mean depth per face) - used by painter_newell_sancha
    model->faces.z_mean = (Fixed32*)malloc(nf * sizeof(Fixed32));
    if (!model->faces.z_mean) {
        printf("Error: Unable to allocate memory for face z_mean array\n");
        keypress();
        if (model->faces.z_mean) free(model->faces.z_mean);
        free(model->vertices.x);
        free(model->vertices.y);
        free(model->vertices.z);
        free(model->vertices.xo);
        free(model->vertices.yo);
        free(model->vertices.zo);
        free(model->vertices.x2d);
        free(model->vertices.y2d);
        free(model->faces.vertex_count);
        free(model->faces.vertex_indices_buffer);
        free(model->faces.vertex_indices_ptr);
        if (model->faces.z_min) free(model->faces.z_min);
        if (model->faces.z_max) free(model->faces.z_max);
        free(model);
        return NULL;
    }

    // Allocate plane coefficient arrays (normalized normals + d term)
    model->faces.plane_a = (Fixed32*)malloc(nf * sizeof(Fixed32));
    model->faces.plane_b = (Fixed32*)malloc(nf * sizeof(Fixed32));
    model->faces.plane_c = (Fixed32*)malloc(nf * sizeof(Fixed32));
    model->faces.plane_d = (Fixed32*)malloc(nf * sizeof(Fixed32));
    if (!model->faces.plane_a || !model->faces.plane_b || !model->faces.plane_c || !model->faces.plane_d) {
        printf("Error: Unable to allocate memory for face plane arrays\n");
        keypress();
        if (model->faces.plane_a) free(model->faces.plane_a);
        if (model->faces.plane_b) free(model->faces.plane_b);
        if (model->faces.plane_c) free(model->faces.plane_c);
        if (model->faces.plane_d) free(model->faces.plane_d);
        if (model->faces.z_mean) free(model->faces.z_mean);
        if (model->faces.z_min) free(model->faces.z_min);
        if (model->faces.z_max) free(model->faces.z_max);
        free(model->vertices.x);
        free(model->vertices.y);
        free(model->vertices.z);
        free(model->vertices.xo);
        free(model->vertices.yo);
        free(model->vertices.zo);
        free(model->vertices.x2d);
        free(model->vertices.y2d);
        free(model->faces.vertex_count);
        free(model->faces.vertex_indices_buffer);
        free(model->faces.vertex_indices_ptr);
        free(model->faces.z_max);
        free(model);
        return NULL;
    }

    // Allocate cached 2D bounding boxes arrays
    model->faces.minx = (int*)malloc(nf * sizeof(int));
    model->faces.maxx = (int*)malloc(nf * sizeof(int));
    model->faces.miny = (int*)malloc(nf * sizeof(int));
    model->faces.maxy = (int*)malloc(nf * sizeof(int));
    if (!model->faces.minx || !model->faces.maxx || !model->faces.miny || !model->faces.maxy) {
        printf("Error: Unable to allocate memory for face bounding box arrays\n");
        keypress();
        if (model->faces.minx) free(model->faces.minx);
        if (model->faces.maxx) free(model->faces.maxx);
        if (model->faces.miny) free(model->faces.miny);
        if (model->faces.maxy) free(model->faces.maxy);
        if (model->faces.plane_a) free(model->faces.plane_a);
        if (model->faces.plane_b) free(model->faces.plane_b);
        if (model->faces.plane_c) free(model->faces.plane_c);
        if (model->faces.plane_d) free(model->faces.plane_d);
        if (model->faces.z_mean) free(model->faces.z_mean);
        if (model->faces.z_min) free(model->faces.z_min);
        if (model->faces.z_max) free(model->faces.z_max);
        free(model->vertices.x);
        free(model->vertices.y);
        free(model->vertices.z);
        free(model->vertices.xo);
        free(model->vertices.yo);
        free(model->vertices.zo);
        free(model->vertices.x2d);
        free(model->vertices.y2d);
        free(model->faces.vertex_count);
        free(model->faces.vertex_indices_buffer);
        free(model->faces.vertex_indices_ptr);
        free(model->faces.z_max);
        free(model);
        return NULL;
    }
    
    // Allocate display_flag array: nf * 4 bytes = 24KB
    model->faces.display_flag = (int*)malloc(nf * sizeof(int));
    if (!model->faces.display_flag) {
        printf("Error: Unable to allocate memory for face display_flag array\n");
        keypress();
        free(model->vertices.x);
        free(model->vertices.y);
        free(model->vertices.z);
        free(model->vertices.xo);
        free(model->vertices.yo);
        free(model->vertices.zo);
        free(model->vertices.x2d);
        free(model->vertices.y2d);
        free(model->faces.vertex_count);
        free(model->faces.vertex_indices_buffer);
        free(model->faces.vertex_indices_ptr);
        free(model->faces.z_max);
        if (model->faces.z_mean) free(model->faces.z_mean);
        free(model);
        return NULL;
    }
    
    // Initialize structure
    model->faces.vertex_countHandle = NULL;
    model->faces.vertex_indicesBufferHandle = NULL;
    model->faces.vertex_indicesPtrHandle = NULL;
    model->faces.z_maxHandle = NULL;
    model->faces.display_flagHandle = NULL;
    model->faces.sorted_face_indicesHandle = NULL;
    
    // Allocate sorted_face_indices array: nf * 4 bytes = 24KB max
    model->faces.sorted_face_indices = (int*)malloc(nf * sizeof(int));
    if (!model->faces.sorted_face_indices) {
        printf("Error: Unable to allocate memory for sorted_face_indices array\n");
        keypress();
        free(model->vertices.x);
        free(model->vertices.y);
        free(model->vertices.z);
        free(model->vertices.xo);
        free(model->vertices.yo);
        free(model->vertices.zo);
        free(model->vertices.x2d);
        free(model->vertices.y2d);
        free(model->faces.vertex_count);
        free(model->faces.vertex_indices_buffer);
        free(model->faces.vertex_indices_ptr);
        if (model->faces.z_min) free(model->faces.z_min);
        if (model->faces.z_max) free(model->faces.z_max);
        if (model->faces.z_mean) free(model->faces.z_mean);
        free(model->faces.display_flag);
        free(model);
        return NULL;
    }
    
    // Initialize structure
    model->faces.vertex_countHandle = NULL;
    model->faces.vertex_indicesBufferHandle = NULL;
    model->faces.vertex_indicesPtrHandle = NULL;
    model->faces.z_maxHandle = NULL;
    model->faces.display_flagHandle = NULL;
    model->faces.sorted_face_indicesHandle = NULL;
    model->faces.total_indices = 0;
    
    return model;
}

/**
 * MEMORY CLEANUP - DESTROY MODEL 3D
 * ==================================
 * 
 * This function frees all memory allocated by createModel3D().
 * Must be called when the model is no longer needed to prevent memory leaks.
 * 
 * CLEANUP SEQUENCE:
 * 1. Free all vertex arrays (x, y, z, xo, yo, zo, x2d, y2d)
 * 2. Free all face arrays (vertex_count, vertex_indices_buffer, vertex_indices_ptr, z_max, display_flag, sorted_face_indices)
 * 3. Free main Model3D structure
 */
void destroyModel3D(Model3D* model) {
    if (model != NULL) {
        // Free all vertex arrays
        if (model->vertices.x) free(model->vertices.x);
        if (model->vertices.y) free(model->vertices.y);
        if (model->vertices.z) free(model->vertices.z);
        if (model->vertices.xo) free(model->vertices.xo);
        if (model->vertices.yo) free(model->vertices.yo);
        if (model->vertices.zo) free(model->vertices.zo);
        if (model->vertices.x2d) free(model->vertices.x2d);
        if (model->vertices.y2d) free(model->vertices.y2d);
        
        // Free all face arrays (now simplified with packed buffer)
        if (model->faces.vertex_count) free(model->faces.vertex_count);
        if (model->faces.vertex_indices_buffer) free(model->faces.vertex_indices_buffer);
        if (model->faces.vertex_indices_ptr) free(model->faces.vertex_indices_ptr);
        if (model->faces.z_max) free(model->faces.z_max);
        if (model->faces.z_mean) free(model->faces.z_mean);
        if (model->faces.z_min) free(model->faces.z_min);
        if (model->faces.plane_a) free(model->faces.plane_a);
        if (model->faces.plane_b) free(model->faces.plane_b);
        if (model->faces.plane_c) free(model->faces.plane_c);
        if (model->faces.plane_d) free(model->faces.plane_d);
        if (model->faces.minx) free(model->faces.minx);
        if (model->faces.maxx) free(model->faces.maxx);
        if (model->faces.miny) free(model->faces.miny);
        if (model->faces.maxy) free(model->faces.maxy);
        if (model->faces.display_flag) free(model->faces.display_flag);
        if (model->faces.sorted_face_indices) free(model->faces.sorted_face_indices);

        // Free optional buffers and backups
        if (model->orig_x) free(model->orig_x);
        if (model->orig_y) free(model->orig_y);
        if (model->orig_z) free(model->orig_z);
        if (model->radius_buf) free(model->radius_buf);
        if (model->coord_buf) free(model->coord_buf);
        
        // Free main structure
        free(model);
    }
}

segment "code";
/**
 * COMPLETE 3D MODEL LOADING
 * ==========================
 * 
 * This function coordinates the complete loading of an OBJ file
 * by successively calling the vertex and face reading functions.
 * 
 * LOADING PIPELINE:
 * 1. Input parameter validation
 * 2. Read vertices from file
 * 3. Read faces from file  
 * 4. Update counters in structure
 * 
 * ERROR HANDLING:
 * - Vertex reading failure: immediate stop
 * - Face reading failure: warning but continue
 *   (vertices-only model remains usable)
 */
int loadModel3D(Model3D* model, const char* filename) {
    // Input parameter validation
    if (model == NULL || filename == NULL) {
        return -1;  // Invalid parameters
    }
    
    // Step 1: Read vertices from OBJ file
    // --- MAJ du compteur global pour la vérification des indices de faces ---
    readVertices_last_count = model->vertices.vertex_count;
    
    int vcount = readVertices(filename, &model->vertices, MAX_VERTICES);
    if (vcount < 0) {
        return -1;  // Critical failure: unable to read vertices
    }
    model->vertices.vertex_count = vcount;
    
    // Step 2: Read faces from OBJ file (using chunked allocation)
    // This function handles reading faces into 2 chunks transparently
    int fcount = readFaces_model(filename, model);
    if (fcount < 0) {
        // Critical failure: unable to read faces
        printf("\nWarning: Unable to read faces\n");
        model->faces.face_count = 0;  // No faces available
    } else {
        model->faces.face_count = fcount;
    }
    
    // After loading vertices, compute bounding sphere for fast auto-fit
    computeModelBoundingSphere(model);
    return 0;  // Success: model loaded (with or without faces)
}

// Compute bounding sphere (centroid + max radius) - O(n) once at load
void computeModelBoundingSphere(Model3D* model) {
    if (!model) return;
    VertexArrays3D* vtx = &model->vertices;
    int n = vtx->vertex_count;
    if (n <= 0) { model->bs_valid = 0; return; }
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (int i = 0; i < n; ++i) {
        cx += FIXED_TO_FLOAT(vtx->x[i]);
        cy += FIXED_TO_FLOAT(vtx->y[i]);
        cz += FIXED_TO_FLOAT(vtx->z[i]);
    }
    cx /= (double)n; cy /= (double)n; cz /= (double)n;
    double max_r2 = 0.0;
    for (int i = 0; i < n; ++i) {
        double dx = FIXED_TO_FLOAT(vtx->x[i]) - cx;
        double dy = FIXED_TO_FLOAT(vtx->y[i]) - cy;
        double dz = FIXED_TO_FLOAT(vtx->z[i]) - cz;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 > max_r2) max_r2 = r2;
    }
    model->bs_cx = (float)cx;
    model->bs_cy = (float)cy;
    model->bs_cz = (float)cz;
    model->bs_r = (float)sqrt(max_r2);
    model->bs_valid = 1;
}

// Compute distance estimate directly from bounding sphere (cheap, O(1))
Fixed32 computeDistanceFromBoundingSphere(Model3D* model, float margin) {
    if (!model || !model->bs_valid) return FLOAT_TO_FIXED(30.0f);
    float r = model->bs_r;
    const float scale = 100.0f;
    float screenHalfX = (float)CENTRE_X * margin;
    float screenHalfY = (float)CENTRE_Y * margin;
    float minScreenHalf = fminf(screenHalfX, screenHalfY);
    float d_proj = (minScreenHalf > 0.0f) ? (r * scale) / minScreenHalf : 1.0f;
    float d = d_proj + r + 1.0f;
    d *= 1.05f;
    if (!isfinite(d) || d < 1.0f) d = 1.0f;
    return FLOAT_TO_FIXED(d);
}

// ============================================================================
//                    USER INTERFACE FUNCTIONS
// ============================================================================

/**
 * OBSERVER PARAMETER INPUT
 * ========================
 * 
 * This function presents a text interface to allow the
 * user to specify 3D visualization parameters.
 * 
 * REQUESTED PARAMETERS:
 * - Horizontal angle: rotation around Y-axis (left/right view)
 * - Vertical angle: rotation around X-axis (up/down view)
 * - Distance: observer distance (zoom)
 * - Screen rotation angle: final rotation in 2D plane
 * 
 * ERROR HANDLING:
 * - Default values if input failure
 * - Automatic string->Fixed32 conversion with atof() then FLOAT_TO_FIXED
 */
void getObserverParams(ObserverParams* params, Model3D* model) {
    char input[50];  // Buffer for user input
    
    // Display section header
    printf("\nObserver parameters:\n");
    printf("============================\n");
    printf("(Press ENTER to use default values - Distance will auto-fit the model)\n");
    
    // Input horizontal angle (rotation around Y)
    printf("Horizontal angle (degrees, default 30): ");
    if (fgets(input, sizeof(input), stdin) != NULL) {
        // Remove newline
        input[strcspn(input, "\n")] = 0;
        if (strlen(input) == 0) {
            params->angle_h = 30;     // Default value if ENTER (degrees)
        } else {
            params->angle_h = atoi(input);  // Parse integer degrees from input
        }
    } else {
        params->angle_h = 30;         // Default value if error (degrees)
    }
    
    // Input vertical angle (rotation around X)  
    printf("Vertical angle (degrees, default 20): ");
    if (fgets(input, sizeof(input), stdin) != NULL) {
        // Remove newline
        input[strcspn(input, "\n")] = 0;
        if (strlen(input) == 0) {
            params->angle_v = 20;     // Default value if ENTER (degrees)
        } else {
            params->angle_v = atoi(input);  // Parse integer degrees from input
        }
    } else {
        params->angle_v = 20;         // Default value if error (degrees)
    }
    

    // Input screen rotation angle (final 2D rotation)
    printf("Screen rotation angle (degrees, default 0): ");
    if (fgets(input, sizeof(input), stdin) != NULL) {
        // Remove newline
        input[strcspn(input, "\n")] = 0;
        if (strlen(input) == 0) {
            params->angle_w = 0;      // Default value if ENTER (degrees)
        } else {
            params->angle_w = atoi(input);  // Parse integer degrees from input
        }
    } else {
        params->angle_w = 0;          // No rotation by default (degrees)
    }

    // Input observation distance (zoom/perspective).
    // Press ENTER = auto-scale + center using sphere-based fit (default target),
    // or enter a numeric value to use that distance directly (no scaling).
    printf("Distance (ENTER = auto-scale & center, or enter a value): ");

    if (fgets(input, sizeof(input), stdin) != NULL) {
        // Remove newline
        input[strcspn(input, "\n")] = 0;
        if (strlen(input) == 0) {
            // User pressed ENTER: perform full auto-fit + centering using default target
            if (model != NULL) {
                float default_target = 200.0f;
                fitModelToView(model, params, default_target, 0.92f, 0.99f, 1);

                printf("Auto-fit/sphere applied (factor %.4f). Press 'r' to revert or +/- to adjust.\n", FIXED_TO_FLOAT(model->auto_scale));
            } else {
                params->distance = FLOAT_TO_FIXED(30.0);
            }
        } else {
            // User provided a distance value -> use it directly (no auto-scale)
            params->distance = FLOAT_TO_FIXED(atof(input)); // String->Fixed32 conversion
        }
    } else {
        // fgets failed; default behavior: if we have a model, auto-fit, else fallback distance
        if (model != NULL) {
            float default_target = 200.0f;
            fitModelToView(model, params, default_target, 0.92f, 0.99f, 1);
            printf("Auto-fit/sphere applied (factor %.4f). Press 'r' to revert or +/- to adjust.\n", FIXED_TO_FLOAT(model->auto_scale));
        } else {
            params->distance = FLOAT_TO_FIXED(30.0);        // Default distance: balanced view (Fixed Point)
        }
    }

    // Debug: show parsed observer angles in degrees
    if (!PERFORMANCE_MODE)
    {
    printf("Observer angles (degrees) - H: %d, V: %d, W: %d\n", params->angle_h, params->angle_v, params->angle_w);
    }
}

/**
 * ULTRA-FAST FUNCTION: Combined Transformation + Projection
 * ==========================================================
 */
void processModelFast(Model3D* model, ObserverParams* params, const char* filename) {
    int i;

    Fixed32 cos_h, sin_h, cos_v, sin_v, cos_w, sin_w;
    Fixed32 x, y, z, zo, xo, yo;
    Fixed32 inv_zo, x2d_temp, y2d_temp;
    
    // Direct table access - ultra-fast! (no function calls)
    cos_h = cos_deg_int(params->angle_h);
    sin_h = sin_deg_int(params->angle_h);
    cos_v = cos_deg_int(params->angle_v);
    sin_v = sin_deg_int(params->angle_v);
    cos_w = cos_deg_int(params->angle_w);
    sin_w = sin_deg_int(params->angle_w);

    // Pre-calculate all trigonometric products in Fixed32 - using 64-bit multiply
    const Fixed32 cos_h_cos_v = FIXED_MUL_64(cos_h, cos_v);
    const Fixed32 sin_h_cos_v = FIXED_MUL_64(sin_h, cos_v);
    const Fixed32 cos_h_sin_v = FIXED_MUL_64(cos_h, sin_v);
    const Fixed32 sin_h_sin_v = FIXED_MUL_64(sin_h, sin_v);
    const Fixed32 scale = INT_TO_FIXED(100); // avoid float conversion
    const Fixed32 centre_x_f = INT_TO_FIXED(CENTRE_X);
    const Fixed32 centre_y_f = INT_TO_FIXED(CENTRE_Y);
    Fixed32 distance = params->distance;
    // Apply auto-scale if any
    //distance = FIXED_MUL_64(distance, INT_TO_FIXED(4));
    
    // 100% Fixed32 loop - ZERO conversions, maximum speed!
    VertexArrays3D* vtx = &model->vertices;
    Fixed32 *x_arr = vtx->x, *y_arr = vtx->y, *z_arr = vtx->z;
    Fixed32 *xo_arr = vtx->xo, *yo_arr = vtx->yo, *zo_arr = vtx->zo;
    int *x2d_arr = vtx->x2d, *y2d_arr = vtx->y2d;
    int vcount = vtx->vertex_count;

    long t_loop_start = GetTick();
    for (i = 0; i < vcount; i++) {
        x = x_arr[i];
        y = y_arr[i];
        z = z_arr[i];
        // 3D transformation in pure Fixed32 (64-bit multiply)
        Fixed32 term1 = FIXED_MUL_64(x, cos_h_cos_v);
        Fixed32 term2 = FIXED_MUL_64(y, sin_h_cos_v);
        Fixed32 term3 = FIXED_MUL_64(z, sin_v);
        zo = FIXED_ADD(FIXED_SUB(FIXED_SUB(FIXED_NEG(term1), term2), term3), distance);
        if (zo > 0) {
            xo = FIXED_ADD(FIXED_NEG(FIXED_MUL_64(x, sin_h)), FIXED_MUL_64(y, cos_h));
            yo = FIXED_ADD(FIXED_SUB(FIXED_NEG(FIXED_MUL_64(x, cos_h_sin_v)), FIXED_MUL_64(y, sin_h_sin_v)), FIXED_MUL_64(z, cos_v));
            zo_arr[i] = zo;
            xo_arr[i] = xo;
            yo_arr[i] = yo;
            inv_zo = FIXED_DIV_64(scale, zo);
            x2d_temp = FIXED_ADD(FIXED_MUL_64(xo, inv_zo), centre_x_f);
            y2d_temp = FIXED_SUB(centre_y_f, FIXED_MUL_64(yo, inv_zo));
            x2d_arr[i] = FIXED_ROUND_TO_INT(FIXED_ADD(FIXED_SUB(FIXED_MUL_64(cos_w, FIXED_SUB(x2d_temp, centre_x_f)), FIXED_MUL_64(sin_w, FIXED_SUB(centre_y_f, y2d_temp))), centre_x_f));
            y2d_arr[i] = FIXED_ROUND_TO_INT(FIXED_SUB(centre_y_f, FIXED_ADD(FIXED_MUL_64(sin_w, FIXED_SUB(x2d_temp, centre_x_f)), FIXED_MUL_64(cos_w, FIXED_SUB(centre_y_f, y2d_temp)))));
                // XXX
        //     x2d_arr[i] = (x2d_arr[i]- 160)*4 + 160; // stretch to full resolution
        //     y2d_arr[i] = (y2d_arr[i]- 100)*4 + 100;
        } else {
            zo_arr[i] = zo;
            xo_arr[i] = 0;
            yo_arr[i] = 0;
            x2d_arr[i] = -1;
            y2d_arr[i] = -1;
        }
    }
    long t_loop_end = GetTick();
    if (!PERFORMANCE_MODE) {
        long elapsed_loop = t_loop_end - t_loop_start;
        double ms_loop = ((double)elapsed_loop * 1000.0) / 60.0; // 60 ticks per second
        printf("[TIMING] transform+project loop: %ld ticks (%.2f ms)\n", elapsed_loop, ms_loop);
    }
    
    // Face sorting after transformation
    long t_start, t_end;

    t_start = GetTick();
    calculateFaceDepths(model, NULL, model->faces.face_count);
    t_end = GetTick();
    if (!PERFORMANCE_MODE)
    {
        long elapsed = t_end - t_start;
        double ms = ((double)elapsed * 1000.0) / 60.0; // 60 ticks per second
        printf("[TIMING] calculateFaceDepths: %ld ticks (%.2f ms)\n", elapsed, ms);
    }

    
    // CRITICAL: Reset sorted_face_indices before each sort to prevent corruption
    for (i = 0; i < model->faces.face_count; i++) {
        model->faces.sorted_face_indices[i] = i;
    }

    // painter_newell_sancha (remplace sortFacesByDepth)
    t_start = GetTick();
    if (painterFastMode) painter_newell_sancha_fast(model, model->faces.face_count);
    else painter_newell_sancha(model, model->faces.face_count);
    t_end = GetTick();

    if (!PERFORMANCE_MODE)
    {
        long elapsed = t_end - t_start;
        double ms = ((double)elapsed * 1000.0) / 60.0; // 60 ticks per second
        printf("[TIMING] %s: %ld ticks (%.2f ms)\n", painterFastMode ? "painter_newell_sancha_fast" : "painter_newell_sancha", elapsed, ms);
        keypress();
    }
}

// Lightweight wireframe processing: only transform & project vertices, set face visibility
// No per-face depth calculations or sorting performed here for maximum speed in wireframe mode
void processModelWireframe(Model3D* model, ObserverParams* params, const char* filename) {
    int i, j;
    Fixed32 cos_h, sin_h, cos_v, sin_v, cos_w, sin_w;
    Fixed32 x, y, z, zo, xo, yo;
    Fixed32 inv_zo, x2d_temp, y2d_temp;

    cos_h = cos_deg_int(params->angle_h);
    sin_h = sin_deg_int(params->angle_h);
    cos_v = cos_deg_int(params->angle_v);
    sin_v = sin_deg_int(params->angle_v);
    cos_w = cos_deg_int(params->angle_w);
    sin_w = sin_deg_int(params->angle_w);

    const Fixed32 cos_h_cos_v = FIXED_MUL_64(cos_h, cos_v);
    const Fixed32 sin_h_cos_v = FIXED_MUL_64(sin_h, cos_v);
    const Fixed32 cos_h_sin_v = FIXED_MUL_64(cos_h, sin_v);
    const Fixed32 sin_h_sin_v = FIXED_MUL_64(sin_h, sin_v);
    const Fixed32 scale = INT_TO_FIXED(100);
    const Fixed32 centre_x_f = INT_TO_FIXED(CENTRE_X);
    const Fixed32 centre_y_f = INT_TO_FIXED(CENTRE_Y);
    const Fixed32 distance = params->distance;

    VertexArrays3D* vtx = &model->vertices;
    Fixed32 *x_arr = vtx->x, *y_arr = vtx->y, *z_arr = vtx->z;
    Fixed32 *xo_arr = vtx->xo, *yo_arr = vtx->yo, *zo_arr = vtx->zo;
    int *x2d_arr = vtx->x2d, *y2d_arr = vtx->y2d;
    int vcount = vtx->vertex_count;

    // Local copies for speed
    FaceArrays3D* faces = &model->faces;
    int *vertex_indices_buffer = faces->vertex_indices_buffer;
    int *vertex_indices_ptr = faces->vertex_indices_ptr;
    int *face_vertex_count = faces->vertex_count;

    for (i = 0; i < vcount; i++) {
        x = x_arr[i];
        y = y_arr[i];
        z = z_arr[i];
        Fixed32 term1 = FIXED_MUL_64(x, cos_h_cos_v);
        Fixed32 term2 = FIXED_MUL_64(y, sin_h_cos_v);
        Fixed32 term3 = FIXED_MUL_64(z, sin_v);
        zo = FIXED_ADD(FIXED_SUB(FIXED_SUB(FIXED_NEG(term1), term2), term3), distance);
        if (zo > 0) {
            // compute projected xy directly into x2d/y2d and store intermediate observer coords for depth tests
            Fixed32 xo_local = FIXED_ADD(FIXED_NEG(FIXED_MUL_64(x, sin_h)), FIXED_MUL_64(y, cos_h));
            Fixed32 yo_local = FIXED_ADD(FIXED_SUB(FIXED_NEG(FIXED_MUL_64(x, cos_h_sin_v)), FIXED_MUL_64(y, sin_h_sin_v)), FIXED_MUL_64(z, cos_v));
            inv_zo = FIXED_DIV_64(scale, zo);
            Fixed32 tmp_x = FIXED_ADD(FIXED_MUL_64(xo_local, inv_zo), centre_x_f);
            Fixed32 tmp_y = FIXED_SUB(centre_y_f, FIXED_MUL_64(yo_local, inv_zo));
            // apply screen rotation and round
            x2d_arr[i] = FIXED_ROUND_TO_INT(FIXED_ADD(FIXED_SUB(FIXED_MUL_64(cos_w, FIXED_SUB(tmp_x, centre_x_f)), FIXED_MUL_64(sin_w, FIXED_SUB(centre_y_f, tmp_y))), centre_x_f));
            y2d_arr[i] = FIXED_ROUND_TO_INT(FIXED_SUB(centre_y_f, FIXED_ADD(FIXED_MUL_64(sin_w, FIXED_SUB(tmp_x, centre_x_f)), FIXED_MUL_64(cos_w, FIXED_SUB(centre_y_f, tmp_y)))));
            // Store observer-space coordinates so subsequent face tests see valid values
            zo_arr[i] = zo;
            xo_arr[i] = xo_local;
            yo_arr[i] = yo_local;
        } else {
            // negative zo (behind camera) — mark as invalid projection and store zo<=0
            zo_arr[i] = zo;
            xo_arr[i] = 0;
            yo_arr[i] = 0;
            x2d_arr[i] = -1;
            y2d_arr[i] = -1;
        }
    }

    // Set simple visibility flag per face: visible if any vertex projected on-screen (x2d != -1)
    for (i = 0; i < faces->face_count; ++i) {
        int offset = vertex_indices_ptr[i];
        int vcount_face = face_vertex_count[i];
        int *indices_base = &vertex_indices_buffer[offset];
        int visible = 0;
        for (j = 0; j < vcount_face; ++j) {
            int vi = indices_base[j] - 1;
            if (vi >= 0 && vi < vcount && x2d_arr[vi] != -1) { visible = 1; break; }
        }
        faces->display_flag[i] = visible;
        faces->sorted_face_indices[i] = i; // identity order; no sorting required
    }
}

// ============================================================================
//                    BASIC FUNCTION IMPLEMENTATIONS
// ============================================================================

/**
 * VERTEX READING FROM OBJ FILE
 * ============================
 * 
 * This function parses an OBJ format file to extract
 * vertices (3D points). It searches for lines starting with "v "
 * and extracts X, Y, Z coordinates.
 * 
 * OBJ FORMAT FOR VERTICES:
 *   v 1.234 5.678 9.012
 *   v -2.5 0.0 3.14159
 * 
 * ALGORITHM:
 * 1. Open file in read mode
 * 2. Read line by line with fgets()
 * 3. Detect "v " lines with character verification
 * 4. Extract coordinates with sscanf()
 * 5. Store in array with bounds checking
 * 6. Progressive display for user feedback
 * 
 * ERROR HANDLING:
 * - File opening verification
 * - Array overflow protection
 * - Coordinate format validation
 */
int readVertices(const char* filename, VertexArrays3D* vtx, int max_vertices) {
    FILE *file;
    char line[MAX_LINE_LENGTH];
    int line_number = 1;
    int vertex_count = 0;
    
    // Open file in read mode
    file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error: Unable to open file '%s'\n", filename);
        printf("Check that the file exists.\n\n");
        return -1;  // Return -1 on error
    }
    
    printf("\nReading vertices from file...'%s':\n", filename);
    
    // Read file line by line
    while (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] == 'v' && line[1] == ' ') {
            if (vertex_count < max_vertices) {
                float x, y, z;  // Temporary reading in float
                if (sscanf(line + 2, "%f %f %f", &x, &y, &z) == 3) {
                    vtx->x[vertex_count] = FLOAT_TO_FIXED(x);
                    // OBJ files often use Z as up; convert to viewer convention by swapping axes
                    vtx->y[vertex_count] = FLOAT_TO_FIXED(z);   // treat OBJ Z as up
                    vtx->z[vertex_count] = FLOAT_TO_FIXED(y);  // rotate -90° around X, corrected orientation
                    vertex_count++;
                    if (vertex_count % 10 == 0) printf("..");

                } else {
                    printf("[\nDEBUG] readVertices: sscanf failed at line %d: %s\n", line_number, line);
                    keypress();
                }
            } else {
                printf("\n[DEBUG] readVertices: vertex limit reached (%d)\n", max_vertices);
                keypress();
            }
        }
        line_number++;
    }
    printf("\n");
    printf("Reading vertices finished : %d vertices read.\n", vertex_count);

    // Close file
    fclose(file);
    
    // printf("\n\nAnalyse terminee. %d lignes lues.\n", line_number - 1);
    return vertex_count;  // Return the number of vertices read
}

// Function to read faces into parallel arrays in FaceArrays3D structure
int readFaces_model(const char* filename, Model3D* model) {
    FILE *file;
    char line[MAX_LINE_LENGTH];
    int line_number = 1;
    int face_count = 0;
    int i;
    
    // Validate model structure
    if (model == NULL || model->faces.vertex_count == NULL) {
        printf("Error: Invalid model structure for readFaces_model\n");
        return -1;
    }
    
    // Open file in read mode
    file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error: Unable to open file '%s' to read faces\n", filename);
        return -1;
    }
    
    printf("\nReading faces from file '%s' :\n", filename);
    
    int buffer_pos = 0;  // Current position in the packed buffer
    
    // Read file line by line
    while (fgets(line, sizeof(line), file) != NULL) {
        // Check if line starts with "f " (face)
        if (line[0] == 'f' && line[1] == ' ') {
            if (face_count < MAX_FACES) {
                // Initialize face data
                model->faces.vertex_count[face_count] = 0;
                model->faces.display_flag[face_count] = 1;  // Displayable by default
                model->faces.vertex_indices_ptr[face_count] = buffer_pos;  // Store offset to this face's indices
                
                // Parse vertices from this face
                char *ptr = line + 2;  // Start after "f "
                int temp_indices[MAX_FACE_VERTICES];
                int temp_vertex_count = 0;
                int invalid_index_found = 0;
                
                // Parse character by character
                while (*ptr != '\0' && *ptr != '\n' && temp_vertex_count < MAX_FACE_VERTICES) {
                    // Skip spaces and tabs
                    while (*ptr == ' ' || *ptr == '\t') ptr++;
                    
                    if (*ptr == '\0' || *ptr == '\n') break;
                    
                    // Read the number
                    int vertex_index = 0;
                    while (*ptr >= '0' && *ptr <= '9') {
                        vertex_index = vertex_index * 10 + (*ptr - '0');
                        ptr++;
                    }
                    
                    // Skip texture/normal data (after /)
                    while (*ptr != '\0' && *ptr != ' ' && *ptr != '\t' && *ptr != '\n') {
                        ptr++;
                    }
                    
                    // Validate vertex index
                    if (vertex_index >= 1) {
                        if (vertex_index > readVertices_last_count) {
                            // Index out of bounds
                            invalid_index_found = 1;
                        }
                        temp_indices[temp_vertex_count] = vertex_index;
                        temp_vertex_count++;
                    }
                }
                
                // Check for errors
                if (invalid_index_found) {
                    printf("\nERROR: Face at line %d references vertex index > %d vertices\n", 
                           line_number, readVertices_last_count);
                    fclose(file);
                    return -1;
                } else {
                    // Store valid indices into the packed buffer
                    for (i = 0; i < temp_vertex_count; i++) {
                        model->faces.vertex_indices_buffer[buffer_pos++] = temp_indices[i];
                    }
                    model->faces.vertex_count[face_count] = temp_vertex_count;
                    model->faces.total_indices += temp_vertex_count;
                    
                    // Skip faces with zero vertices
                    if (model->faces.vertex_count[face_count] > 0) {
                        face_count++;
                        if (face_count % 10 == 0) {printf(".");}
                    } else {
                        printf("     -> WARNING: Face without valid vertices ignored\n");
                    }
                }
            } else {
                printf("     -> WARNING: Face limit reached (%d)\n", MAX_FACES);
            }
        }
        
        line_number++;
    }
    
    // Close file
    fclose(file);
    
    model->faces.face_count = face_count;
    
    // Initialize sorted_face_indices with identity mapping (will be sorted later)
    for (i = 0; i < face_count; i++) {
        model->faces.sorted_face_indices[i] = i;
    }
    
    printf("\nReading faces finished : %d faces read.\n", face_count);
    return face_count;
}

// Function to project 3D coordinates onto 2D screen - FIXED POINT VERSION
void projectTo2D(VertexArrays3D* vtx, int angle_w_deg) {
    int i;
    Fixed32 cos_w, sin_w;
    Fixed32 x2d_temp, y2d_temp;
    cos_w = cos_deg_int(angle_w_deg);
    sin_w = sin_deg_int(angle_w_deg);
    const Fixed32 scale = INT_TO_FIXED(100);
    const Fixed32 centre_x_f = INT_TO_FIXED(CENTRE_X);
    const Fixed32 centre_y_f = INT_TO_FIXED(CENTRE_Y);

    for (i = 0; i < vtx->vertex_count; i++) {
        if (vtx->zo[i] > 0) {
            Fixed32 xo = vtx->xo[i];
            Fixed32 yo = vtx->yo[i];
            Fixed32 inv_zo = FIXED_DIV_64(scale, vtx->zo[i]);
            x2d_temp = FIXED_ADD(FIXED_MUL_64(xo, inv_zo), centre_x_f);
            y2d_temp = FIXED_SUB(centre_y_f, FIXED_MUL_64(yo, inv_zo));
            vtx->x2d[i] = FIXED_ROUND_TO_INT(FIXED_ADD(FIXED_SUB(FIXED_MUL_64(cos_w, FIXED_SUB(x2d_temp, centre_x_f)), FIXED_MUL_64(sin_w, FIXED_SUB(centre_y_f, y2d_temp))), centre_x_f));
            vtx->y2d[i] = FIXED_ROUND_TO_INT(FIXED_SUB(centre_y_f, FIXED_ADD(FIXED_MUL_64(sin_w, FIXED_SUB(x2d_temp, centre_x_f)), FIXED_MUL_64(cos_w, FIXED_SUB(centre_y_f, y2d_temp)))));
        } else {
            vtx->x2d[i] = -1;
            vtx->y2d[i] = -1;
        }
    }
}

/**
 * CALCULATING MINIMUM FACE DEPTHS AND VISIBILITY FLAGS
 * =====================================================
 * 
 * This function calculates for each face:
 * 1. The minimum depth (z_min) of all its vertices in the observer coordinate system
 * 2. The display visibility flag based on vertex positions relative to camera
 * 
 * The z_min value is used for face sorting during rendering (painter's algorithm).

 * We use minimum (closest point) for correct occlusion in the painter's algorithm.
 * The display_flag is used to cull faces that have vertices behind the camera.
 * 
 * PARAMETERS:
 *   vertices   : Array of vertices with coordinates in observer system
 *   faces      : Array of faces to process  
 *   face_count : Number of faces
 * 
 * ALGORITHM:
 *   For each face:
 *   - Initialize z_min with very large value (9999.0)
 *   - Initialize display_flag as true (displayable)
 *   - For each vertex of the face:
 *     * Check if vertex is behind camera (zo <= 0)
 *     * If ANY vertex is behind camera, set display_flag = false
 *     * Update z_min with minimum zo value found (closest vertex)
 *   - Store both z_min and display_flag in the face structure
 * 
 * CULLING LOGIC:
 *   - If ANY vertex has zo <= 0, the entire face is marked as non-displayable
 *   - This prevents rendering artifacts from perspective projection errors
 *   - Improves performance by eliminating faces early in the pipeline
 * 
 * NOTES:
 *   - Must be called AFTER transformToObserver() or processModelFast()
 *   - Uses zo coordinates (observer system depth)
 *   - Lower z_min value means face is closer to camera (should draw first in painter's algorithm)
 *   - display_flag = 1 means visible, 0 means hidden (behind camera)
 */
void calculateFaceDepths(Model3D* model, Face3D* faces, int face_count) {
    int i, j;
    VertexArrays3D* vtx = &model->vertices;
    FaceArrays3D* face_arrays = &model->faces;
    
    for (i = 0; i < face_count; i++) {
        Fixed32 z_min = FLOAT_TO_FIXED(9999.0);  // Initialize to very large value (closest)
        Fixed32 z_max = FLOAT_TO_FIXED(-9999.0); // Initialize to very small value (farthest)
        int display_flag = 1;
        Fixed32 sum = 0;
        int n = face_arrays->vertex_count[i];
        int minx = 9999, maxx = -9999, miny = 9999, maxy = -9999;
        
        // Access indices from the packed buffer using the offset
        int offset = face_arrays->vertex_indices_ptr[i];
        for (j = 0; j < n; j++) {
            int vertex_idx = face_arrays->vertex_indices_buffer[offset + j] - 1;
            if (vertex_idx >= 0) {
                Fixed32 zo = vtx->zo[vertex_idx];
                if (zo < 0) display_flag = 0; // strictly behind camera
                if (zo < z_min) z_min = zo;  // Find minimum (closest)
                if (zo > z_max) z_max = zo;  // Find maximum (farthest)
                sum += zo;
                int x2d = vtx->x2d[vertex_idx];
                int y2d = vtx->y2d[vertex_idx];
                if (x2d < minx) minx = x2d;
                if (x2d > maxx) maxx = x2d;
                if (y2d < miny) miny = y2d;
                if (y2d > maxy) maxy = y2d;
            }
        }
        // Compute plane coefficients (a,b,c,d) using only the first 3 vertices (observer space)
        // Formules (Fixed32 arithmetic, NO normalization):
        // a := y1 * (z2 - z3) + y2 * (z3 - z1) + y3 * (z1 - z2);
        // b := -x1 * (z2 - z3) + x2 * (z1 - z3) - x3 * (z1 - z2);
        // c := x1 * (y2 - y3) - x2 * (y1 - y3) + x3 * (y1 - y2);
        // d := -x1 * (y2 * z3 - y3 * z2) + x2 * (y1 * z3 - y3 * z1) - x3 * (y1 * z2 - y2 * z1);
        Fixed32 a = 0, b = 0, c = 0, d = 0;
        if (!display_flag || n < 3) {
            // face is behind camera or degenerate: zero coefficients
            face_arrays->plane_a[i] = 0;
            face_arrays->plane_b[i] = 0;
            face_arrays->plane_c[i] = 0;
            face_arrays->plane_d[i] = 0;
        } else {
            int idx0 = face_arrays->vertex_indices_buffer[offset] - 1;
            int idx1 = face_arrays->vertex_indices_buffer[offset + 1] - 1;
            int idx2 = face_arrays->vertex_indices_buffer[offset + 2] - 1;
            if (idx0 < 0 || idx1 < 0 || idx2 < 0) {
                face_arrays->plane_a[i] = 0;
                face_arrays->plane_b[i] = 0;
                face_arrays->plane_c[i] = 0;
                face_arrays->plane_d[i] = 0;
            } else {
                Fixed32 x1 = vtx->xo[idx0], y1 = vtx->yo[idx0], z1 = vtx->zo[idx0];
                Fixed32 x2 = vtx->xo[idx1], y2 = vtx->yo[idx1], z2 = vtx->zo[idx1];
                Fixed32 x3 = vtx->xo[idx2], y3 = vtx->yo[idx2], z3 = vtx->zo[idx2];

                // a = y1*(z2-z3) + y2*(z3-z1) + y3*(z1-z2)
                a = FIXED_ADD(FIXED_ADD(FIXED_MUL_64(y1, FIXED_SUB(z2, z3)), FIXED_MUL_64(y2, FIXED_SUB(z3, z1))), FIXED_MUL_64(y3, FIXED_SUB(z1, z2)));

                // b = -x1*(z2-z3) + x2*(z1-z3) - x3*(z1-z2)
                b = FIXED_SUB(FIXED_ADD(FIXED_NEG(FIXED_MUL_64(x1, FIXED_SUB(z2, z3))), FIXED_MUL_64(x2, FIXED_SUB(z1, z3))), FIXED_MUL_64(x3, FIXED_SUB(z1, z2)));

                // c = x1*(y2-y3) - x2*(y1-y3) + x3*(y1-y2)
                c = FIXED_ADD(FIXED_SUB(FIXED_MUL_64(x1, FIXED_SUB(y2, y3)), FIXED_MUL_64(x2, FIXED_SUB(y1, y3))), FIXED_MUL_64(x3, FIXED_SUB(y1, y2)));

                // d = -x1*(y2*z3 - y3*z2) + x2*(y1*z3 - y3*z1) - x3*(y1*z2 - y2*z1)
                Fixed32 t1 = FIXED_SUB(FIXED_MUL_64(y2, z3), FIXED_MUL_64(y3, z2));
                Fixed32 t2 = FIXED_SUB(FIXED_MUL_64(y1, z3), FIXED_MUL_64(y3, z1));
                Fixed32 t3 = FIXED_SUB(FIXED_MUL_64(y1, z2), FIXED_MUL_64(y2, z1));
                d = FIXED_ADD(FIXED_ADD(FIXED_NEG(FIXED_MUL_64(x1, t1)), FIXED_MUL_64(x2, t2)), FIXED_NEG(FIXED_MUL_64(x3, t3)));

                face_arrays->plane_a[i] = a;
                face_arrays->plane_b[i] = b;
                face_arrays->plane_c[i] = c;
                face_arrays->plane_d[i] = d;
            }
        }

        // Diagnostic: report suspiciously negative z_max values to help debug
        if (!PERFORMANCE_MODE) {
            float zmax_f = FIXED_TO_FLOAT(z_max);
            if (zmax_f < -100.0f) {
                printf("[DEBUG] Face %d: z_min=%.2f z_max=%.2f display_flag=%d n=%d\n", i, FIXED_TO_FLOAT(z_min), zmax_f, display_flag, n);
                // Print per-vertex observer-space zo values
                printf("[DEBUG]  vertex zo: ");
                for (int jj = 0; jj < n; ++jj) {
                    int vidx = face_arrays->vertex_indices_buffer[offset + jj] - 1;
                    if (vidx >= 0) printf("(%d: %.2f) ", vidx, FIXED_TO_FLOAT(vtx->zo[vidx]));
                }
                printf("\n");
            }
        }

        face_arrays->z_min[i] = z_min;  // Store minimum depth for this face (closest)
        face_arrays->z_max[i] = z_max;  // Store maximum depth for this face (farthest)
        face_arrays->display_flag[i] = display_flag;
        if (n > 0) {
            face_arrays->z_mean[i] = sum / n;
            face_arrays->minx[i] = minx;
            face_arrays->maxx[i] = maxx;
            face_arrays->miny[i] = miny;
            face_arrays->maxy[i] = maxy;
        } else {
            face_arrays->z_mean[i] = 0;
            face_arrays->z_min[i] = 0;
            face_arrays->z_max[i] = 0;
            face_arrays->minx[i] = 0;
            face_arrays->maxx[i] = 0;
            face_arrays->miny[i] = 0;
            face_arrays->maxy[i] = 0;
        }
    }
}






// Dump face plane coefficients and depth stats to CSV
// Columns: face,a,b,c,d,z_min,z_mean,z_max,vertex_indices
void dumpFaceEquationsCSV(Model3D* model, const char* csv_filename) {
    if (model == NULL || csv_filename == NULL) return;
    FILE* f = fopen(csv_filename, "w");
    if (!f) {
        printf("Error: cannot open '%s' for writing\n", csv_filename);
        return;
    }
    FaceArrays3D* faces = &model->faces;
    int face_count = faces->face_count;
    fprintf(f, "face,a,b,c,d,z_min,z_mean,z_max,vertex_indices\n");
    for (int i = 0; i < face_count; ++i) {
        float a = FIXED_TO_FLOAT(faces->plane_a[i]);
        float b = FIXED_TO_FLOAT(faces->plane_b[i]);
        float c = FIXED_TO_FLOAT(faces->plane_c[i]);
        float d = FIXED_TO_FLOAT(faces->plane_d[i]);
        float zmin = FIXED_TO_FLOAT(faces->z_min[i]);
        float zmean = FIXED_TO_FLOAT(faces->z_mean[i]);
        float zmax = FIXED_TO_FLOAT(faces->z_max[i]);
        fprintf(f, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,\"", i, a, b, c, d, zmin, zmean, zmax);
        int offset = faces->vertex_indices_ptr[i];
        int n = faces->vertex_count[i];
        for (int j = 0; j < n; ++j) {
            if (j) fputc(' ', f);
            int vid = faces->vertex_indices_buffer[offset + j]; // keep OBJ 1-based index
            fprintf(f, "%d", vid);
        }
        fprintf(f, "\"\n");
    }
    fclose(f);
    printf("Wrote face equations to %s (%d faces)\n", csv_filename, face_count);
}

// DEPRECATED: computeDistanceToFit
// ---------------------------------
// This vertex-based function computed an observer distance by scanning all
// vertices (O(n)). It was originally used to auto-fit the model into view.
//
// It has been superseded by the bounding-sphere approach:
//  - computeModelBoundingSphere() computes a centroid+radius once at model load
//  - computeDistanceFromBoundingSphere() estimates distance in O(1) from that sphere
//
// The function is retained here as a fallback and for historical/diagnostic
// purposes. Prefer the bounding-sphere helpers for production paths.

// Compute an observation distance (Fixed32) that fits the model within the view.
// Uses the model bounding box and projection scale to estimate a conservative distance
// so the model fits with the requested margin (0..1).
// Optimized: avoid per-vertex float conversions by computing min/max in Fixed32 then
// convert to float only once per axis (fewer expensive float ops for large models).
Fixed32 computeDistanceToFit(VertexArrays3D* vtx, float margin) {
    if (vtx == NULL || vtx->vertex_count <= 0) return FLOAT_TO_FIXED(30.0f);

    int n = vtx->vertex_count;
    // Initialize mins/maxs with first vertex to avoid large sentinels
    Fixed32 minx_f = vtx->x[0], maxx_f = vtx->x[0];
    Fixed32 miny_f = vtx->y[0], maxy_f = vtx->y[0];
    Fixed32 minz_f = vtx->z[0], maxz_f = vtx->z[0];

    for (int i = 1; i < n; ++i) {
        Fixed32 xi = vtx->x[i];
        Fixed32 yi = vtx->y[i];
        Fixed32 zi = vtx->z[i];
        if (xi < minx_f) minx_f = xi;
        if (xi > maxx_f) maxx_f = xi;
        if (yi < miny_f) miny_f = yi;
        if (yi > maxy_f) maxy_f = yi;
        if (zi < minz_f) minz_f = zi;
        if (zi > maxz_f) maxz_f = zi;
    }

    // Convert deltas to float only once
    float width = FIXED_TO_FLOAT(FIXED_SUB(maxx_f, minx_f));
    float height = FIXED_TO_FLOAT(FIXED_SUB(maxy_f, miny_f));
    float depth = FIXED_TO_FLOAT(FIXED_SUB(maxz_f, minz_f));
    float half_w = width * 0.5f;
    float half_h = height * 0.5f;
    float half_depth = depth * 0.5f;

    const float scale = 100.0f; // matches projection scale
    float screenHalfX = (float)CENTRE_X * margin;
    float screenHalfY = (float)CENTRE_Y * margin;

    float d_for_w = (screenHalfX > 0.0f) ? (half_w * scale) / screenHalfX : 1.0f;
    float d_for_h = (screenHalfY > 0.0f) ? (half_h * scale) / screenHalfY : 1.0f;
    float d = fmaxf(d_for_w, d_for_h);

    if (!isfinite(d) || d < 1.0f) d = 1.0f;

    // add half depth so the center offset is respected and a small safety margin
    d += half_depth + 1.0f;
    return FLOAT_TO_FIXED(d);
}


// ==============================================================
// Auto-scaling helpers
// Non-destructive: subtract center, apply scale, store metadata
// Provide a revert function to restore original coordinates
// ==============================================================
void autoScaleModel(Model3D* model, float target_max_dim, float min_scale, float max_scale, int center_flag) {
    if (model == NULL) return;
    VertexArrays3D* vtx = &model->vertices;
    int n = vtx->vertex_count;
    if (n <= 0) return;

    // If the model was already auto-scaled, revert first to avoid cumulative scaling
    if (model->auto_scaled) {
        revertAutoScaleModel(model);
    }

    Fixed32 minx = vtx->x[0], maxx = vtx->x[0];
    Fixed32 miny = vtx->y[0], maxy = vtx->y[0];
    Fixed32 minz = vtx->z[0], maxz = vtx->z[0];
    for (int i = 1; i < n; ++i) {
        Fixed32 xi = vtx->x[i];
        Fixed32 yi = vtx->y[i];
        Fixed32 zi = vtx->z[i];
        if (xi < minx) minx = xi; if (xi > maxx) maxx = xi;
        if (yi < miny) miny = yi; if (yi > maxy) maxy = yi;
        if (zi < minz) minz = zi; if (zi > maxz) maxz = zi;
    }

    float width = FIXED_TO_FLOAT(FIXED_SUB(maxx, minx));
    float height = FIXED_TO_FLOAT(FIXED_SUB(maxy, miny));
    float depth = FIXED_TO_FLOAT(FIXED_SUB(maxz, minz));
    float maxdim = fmaxf(width, fmaxf(height, depth));
    if (maxdim <= 0.0f) return;

    float scale_f = target_max_dim / maxdim;
    if (scale_f < min_scale) scale_f = min_scale;
    if (scale_f > max_scale) scale_f = max_scale;
    Fixed32 scale_fixed = FLOAT_TO_FIXED(scale_f);

    Fixed32 center_x = FIXED_DIV(FIXED_ADD(minx, maxx), INT_TO_FIXED(2));
    Fixed32 center_y = FIXED_DIV(FIXED_ADD(miny, maxy), INT_TO_FIXED(2));
    Fixed32 center_z = FIXED_DIV(FIXED_ADD(minz, maxz), INT_TO_FIXED(2));

    // Backup original coordinates for exact revert and ensure non-cumulative scaling
    backupModelCoords(model);

    // Apply: subtract center (if requested) then scale, using original coordinates as base
    for (int i = 0; i < n; ++i) {
        Fixed32 ox = model->orig_x ? model->orig_x[i] : vtx->x[i];
        Fixed32 oy = model->orig_y ? model->orig_y[i] : vtx->y[i];
        Fixed32 oz = model->orig_z ? model->orig_z[i] : vtx->z[i];
        if (center_flag) {
            vtx->x[i] = FIXED_MUL_64(FIXED_SUB(ox, center_x), scale_fixed);
            vtx->y[i] = FIXED_MUL_64(FIXED_SUB(oy, center_y), scale_fixed);
            vtx->z[i] = FIXED_MUL_64(FIXED_SUB(oz, center_z), scale_fixed);
        } else {
            vtx->x[i] = FIXED_MUL_64(ox, scale_fixed);
            vtx->y[i] = FIXED_MUL_64(oy, scale_fixed);
            vtx->z[i] = FIXED_MUL_64(oz, scale_fixed);
        }
    }

    model->auto_scale = scale_fixed;
    model->auto_center_x = center_x;
    model->auto_center_y = center_y;
    model->auto_center_z = center_z;
    model->auto_scaled = 1;
    model->auto_centered = center_flag;
}

void revertAutoScaleModel(Model3D* model) {
    if (model == NULL) return;
    VertexArrays3D* vtx = &model->vertices;
    int n = vtx->vertex_count;
    if (n <= 0) { model->auto_scaled = 0; model->auto_centered = 0; model->auto_scale = FIXED_ONE; return; }

    // If we have an exact backup, restore it for perfect revert
    if (model->orig_x && model->orig_y && model->orig_z) {
        for (int i = 0; i < n; ++i) {
            vtx->x[i] = model->orig_x[i];
            vtx->y[i] = model->orig_y[i];
            vtx->z[i] = model->orig_z[i];
        }
        freeBackupModelCoords(model);
    } else {
        // Fallback: inverse the applied scale+center using stored params
        Fixed32 scale = model->auto_scale;
        if (scale != 0) {
            for (int i = 0; i < n; ++i) {
                vtx->x[i] = FIXED_ADD(FIXED_DIV_64(vtx->x[i], scale), model->auto_center_x);
                vtx->y[i] = FIXED_ADD(FIXED_DIV_64(vtx->y[i], scale), model->auto_center_y);
                vtx->z[i] = FIXED_ADD(FIXED_DIV_64(vtx->z[i], scale), model->auto_center_z);
            }
        }
    }

    /* Ensure meta flags and scale are cleared after revert so a future auto-fit starts from
       a clean state and does not attempt to double-revert or apply inverse scaling. */
    model->auto_scaled = 0;
    model->auto_centered = 0;
    model->auto_scale = FIXED_ONE;

    /* Recompute the bounding sphere from the (restored) model coordinates so that
       distance estimations and subsequent auto-fit operations are based on current data. */
    computeModelBoundingSphere(model);

    // free radius buffer when reverting (optional but keeps memory tidy)
    if (model->radius_buf) { free(model->radius_buf); model->radius_buf = NULL; model->radius_buf_capacity = 0; }
}

void backupModelCoords(Model3D* model) {
    if (model == NULL) return;
    VertexArrays3D* vtx = &model->vertices;
    int n = vtx->vertex_count;
    if (n <= 0) return;

    // If previous backup existed, free it first
    if (model->orig_x) { free(model->orig_x); model->orig_x = NULL; }
    if (model->orig_y) { free(model->orig_y); model->orig_y = NULL; }
    if (model->orig_z) { free(model->orig_z); model->orig_z = NULL; }

    model->orig_x = (Fixed32*)malloc(n * sizeof(Fixed32));
    model->orig_y = (Fixed32*)malloc(n * sizeof(Fixed32));
    model->orig_z = (Fixed32*)malloc(n * sizeof(Fixed32));
    if (!model->orig_x || !model->orig_y || !model->orig_z) {
        // If allocation failed, ensure a consistent state and bail out
        freeBackupModelCoords(model);
        return;
    }

    for (int i = 0; i < n; ++i) {
        model->orig_x[i] = vtx->x[i];
        model->orig_y[i] = vtx->y[i];
        model->orig_z[i] = vtx->z[i];
    }
}

void freeBackupModelCoords(Model3D* model) {
    if (model == NULL) return;
    if (model->orig_x) { free(model->orig_x); model->orig_x = NULL; }
    if (model->orig_y) { free(model->orig_y); model->orig_y = NULL; }
    if (model->orig_z) { free(model->orig_z); model->orig_z = NULL; }
    if (model->radius_buf) { free(model->radius_buf); model->radius_buf = NULL; model->radius_buf_capacity = 0; }
}

// Fit model to view using sphere-based metric with percentile trimming
void fitModelToView(Model3D* model, ObserverParams* params, float target_max_dim, float margin, float percentile, int center_flag) {
    if (model == NULL || params == NULL) return;
    VertexArrays3D* vtx = &model->vertices;
    int n = vtx->vertex_count;
    if (n <= 0) return;

    // Compute centroid using original coordinates if available
    // Vertex sampling: sample up to max_samples vertices uniformly (no face traversal)
    const int max_samples = 4096;
    int step = 1;
    int sampleCount = n;
    float minR = 0.0f, maxR = 0.0f; // declared here for scope to be used later
    if (n > max_samples) {
        step = (n + max_samples - 1) / max_samples; // ceil division
        if (step < 1) step = 1;
        sampleCount = (n + step - 1) / step;
    }

    // Compute centroid on sampled vertices (fast, no allocations)
    float cx = 0.0f, cy = 0.0f, cz = 0.0f;
    int sc = 0;
    for (int vi = 0; vi < n; vi += step) {
        Fixed32 xi = model->orig_x ? model->orig_x[vi] : vtx->x[vi];
        Fixed32 yi = model->orig_y ? model->orig_y[vi] : vtx->y[vi];
        Fixed32 zi = model->orig_z ? model->orig_z[vi] : vtx->z[vi];
        cx += FIXED_TO_FLOAT(xi);
        cy += FIXED_TO_FLOAT(yi);
        cz += FIXED_TO_FLOAT(zi);
        sc++;
    }

    if (sc < 3) {
        // Fallback: process the full dataset (rare)
        if (model->coord_buf_capacity < n) {
            if (model->coord_buf) free(model->coord_buf);
            model->coord_buf = (float*)malloc((size_t)n * 3 * sizeof(float));
            if (!model->coord_buf) { model->coord_buf_capacity = 0; return; }
            model->coord_buf_capacity = n;
        }
        float *coords = model->coord_buf;
        cx = cy = cz = 0.0f;
        for (int i = 0; i < n; ++i) {
            Fixed32 xi = model->orig_x ? model->orig_x[i] : vtx->x[i];
            Fixed32 yi = model->orig_y ? model->orig_y[i] : vtx->y[i];
            Fixed32 zi = model->orig_z ? model->orig_z[i] : vtx->z[i];
            float xf = FIXED_TO_FLOAT(xi);
            float yf = FIXED_TO_FLOAT(yi);
            float zf = FIXED_TO_FLOAT(zi);
            coords[i*3 + 0] = xf;
            coords[i*3 + 1] = yf;
            coords[i*3 + 2] = zf;
            cx += xf; cy += yf; cz += zf;
        }
        cx /= (float)n; cy /= (float)n; cz /= (float)n;

        // Ensure radius buffer capacity
        if (model->radius_buf_capacity < n) {
            if (model->radius_buf) free(model->radius_buf);
            model->radius_buf = (float*)malloc((size_t)n * sizeof(float));
            if (!model->radius_buf) { model->radius_buf_capacity = 0; return; }
            model->radius_buf_capacity = n;
        }

        // Fill buffer with squared radii and compute min/max in one pass
        float dx0 = coords[0*3 + 0] - cx;
        float dy0 = coords[0*3 + 1] - cy;
        float dz0 = coords[0*3 + 2] - cz;
        float r20 = dx0*dx0 + dy0*dy0 + dz0*dz0;
        model->radius_buf[0] = r20;
        float minR = r20, maxR = r20;
        for (int i = 1; i < n; ++i) {
            float dx = coords[i*3 + 0] - cx;
            float dy = coords[i*3 + 1] - cy;
            float dz = coords[i*3 + 2] - cz;
            float r2 = dx*dx + dy*dy + dz*dz;
            model->radius_buf[i] = r2;
            if (r2 < minR) minR = r2;
            if (r2 > maxR) maxR = r2;
        }
        sampleCount = n; // use full dataset
    } else {
        // Use sampled vertices (sc samples)
        if (model->radius_buf_capacity < sc) {
            if (model->radius_buf) free(model->radius_buf);
            model->radius_buf = (float*)malloc((size_t)sc * sizeof(float));
            if (!model->radius_buf) { model->radius_buf_capacity = 0; return; }
            model->radius_buf_capacity = sc;
        }
        cx /= (float)sc; cy /= (float)sc; cz /= (float)sc;
        float minR = 0.0f, maxR = 0.0f;
        int ii = 0;
        for (int vi = 0; vi < n; vi += step) {
            Fixed32 xi = model->orig_x ? model->orig_x[vi] : vtx->x[vi];
            Fixed32 yi = model->orig_y ? model->orig_y[vi] : vtx->y[vi];
            Fixed32 zi = model->orig_z ? model->orig_z[vi] : vtx->z[vi];
            float dx = FIXED_TO_FLOAT(xi) - cx;
            float dy = FIXED_TO_FLOAT(yi) - cy;
            float dz = FIXED_TO_FLOAT(zi) - cz;
            float r2 = dx*dx + dy*dy + dz*dz;
            model->radius_buf[ii] = r2;
            if (ii == 0) { minR = maxR = r2; } else { if (r2 < minR) minR = r2; if (r2 > maxR) maxR = r2; }
            ii++;
        }
        sampleCount = sc;
    }

    // Select the percentile squared radius using optimized quickselect (median-of-three pivot)
    int idx = (int)(percentile * sampleCount) - 1;
    if (idx < 0) idx = 0; if (idx >= sampleCount) idx = sampleCount-1;

    // We already computed minR/maxR during radius filling earlier (minR/maxR available)
    if (minR == maxR) {
        float Rf = sqrtf(model->radius_buf[idx]);
        if (Rf <= 0.0f) return;
        double scale_f = (double)(target_max_dim / (2.0f * Rf));
        if (!isfinite(scale_f) || scale_f <= 0.0) return;
        if (scale_f < 0.01) scale_f = 0.01; if (scale_f > 1000.0) scale_f = 1000.0;
        Fixed32 scale_fixed = FLOAT_TO_FIXED((float)scale_f);

        // Backup original coords and apply scale+center
        backupModelCoords(model);

        Fixed32 center_x = FLOAT_TO_FIXED((float)cx);
        Fixed32 center_y = FLOAT_TO_FIXED((float)cy);
        Fixed32 center_z = FLOAT_TO_FIXED((float)cz);

        for (int i = 0; i < n; ++i) {
            Fixed32 ox = model->orig_x[i];
            Fixed32 oy = model->orig_y[i];
            Fixed32 oz = model->orig_z[i];
            if (center_flag) {
                vtx->x[i] = FIXED_MUL_64(FIXED_SUB(ox, center_x), scale_fixed);
                vtx->y[i] = FIXED_MUL_64(FIXED_SUB(oy, center_y), scale_fixed);
                vtx->z[i] = FIXED_MUL_64(FIXED_SUB(oz, center_z), scale_fixed);
            } else {
                vtx->x[i] = FIXED_MUL_64(ox, scale_fixed);
                vtx->y[i] = FIXED_MUL_64(oy, scale_fixed);
                vtx->z[i] = FIXED_MUL_64(oz, scale_fixed);
            }
        }

        model->auto_scale = scale_fixed;
        model->auto_center_x = center_x;
        model->auto_center_y = center_y;
        model->auto_center_z = center_z;
        model->auto_scaled = 1;
        model->auto_centered = center_flag;

        // Update bounding sphere parameters according to applied scaling and set distance from sphere (fast, O(1))
        {
            float scale_f32 = FIXED_TO_FLOAT(scale_fixed);
            if (model->bs_valid) {
                if (center_flag) {
                    model->bs_cx = (model->bs_cx - FIXED_TO_FLOAT(center_x)) * scale_f32;
                    model->bs_cy = (model->bs_cy - FIXED_TO_FLOAT(center_y)) * scale_f32;
                    model->bs_cz = (model->bs_cz - FIXED_TO_FLOAT(center_z)) * scale_f32;
                } else {
                    model->bs_cx = model->bs_cx * scale_f32;
                    model->bs_cy = model->bs_cy * scale_f32;
                    model->bs_cz = model->bs_cz * scale_f32;
                }
                model->bs_r *= scale_f32;
                model->bs_valid = 1;
                params->distance = computeDistanceFromBoundingSphere(model, margin);
            } else {
                // fallback (rare): compute from vertices (slower)
                params->distance = computeDistanceToFit(&model->vertices, margin);
            }
        }
        printf("[FIT] p%g radius=%.3f scale=%.4f applied. New distance: %.2f\n", percentile*100.0, (double)Rf, (double)scale_f, FIXED_TO_FLOAT(params->distance));
        return;
    }
    // Quickselect (iterative) with median-of-three pivot selection on the sample
    int left = 0, right = sampleCount - 1;
    while (left < right) {
        int mid = left + ((right - left) >> 1);
        // median-of-three: pick median of left, mid, right
        float a = model->radius_buf[left], b = model->radius_buf[mid], c = model->radius_buf[right];
        float pivot = b;
        if ((a <= b && b <= c) || (c <= b && b <= a)) pivot = b;
        else if ((b <= a && a <= c) || (c <= a && a <= b)) pivot = a;
        else pivot = c;

        int i = left, j = right;
        while (i <= j) {
            while (model->radius_buf[i] < pivot) i++;
            while (model->radius_buf[j] > pivot) j--;
            if (i <= j) {
                float tmp = model->radius_buf[i]; model->radius_buf[i] = model->radius_buf[j]; model->radius_buf[j] = tmp;
                i++; j--;
            }
        }
        if (idx <= j) right = j;
        else if (idx >= i) left = i;
        else break;
    }

    float R2f = model->radius_buf[idx];
    float Rf = sqrtf(R2f);
    if (Rf <= 0.0f) return;

    // Compute scale from sphere metric
    double scale_f = (double)(target_max_dim / (2.0f * Rf));
    if (!isfinite(scale_f) || scale_f <= 0.0) return;
    if (scale_f < 0.01) scale_f = 0.01; if (scale_f > 1000.0) scale_f = 1000.0;
    Fixed32 scale_fixed = FLOAT_TO_FIXED((float)scale_f);

    // Backup original coords and apply scale+center
    backupModelCoords(model);

    Fixed32 center_x = FLOAT_TO_FIXED((float)cx);
    Fixed32 center_y = FLOAT_TO_FIXED((float)cy);
    Fixed32 center_z = FLOAT_TO_FIXED((float)cz);

    for (int i = 0; i < n; ++i) {
        Fixed32 ox = model->orig_x[i];
        Fixed32 oy = model->orig_y[i];
        Fixed32 oz = model->orig_z[i];
        if (center_flag) {
            vtx->x[i] = FIXED_MUL_64(FIXED_SUB(ox, center_x), scale_fixed);
            vtx->y[i] = FIXED_MUL_64(FIXED_SUB(oy, center_y), scale_fixed);
            vtx->z[i] = FIXED_MUL_64(FIXED_SUB(oz, center_z), scale_fixed);
        } else {
            vtx->x[i] = FIXED_MUL_64(ox, scale_fixed);
            vtx->y[i] = FIXED_MUL_64(oy, scale_fixed);
            vtx->z[i] = FIXED_MUL_64(oz, scale_fixed);
        }
    }

    model->auto_scale = scale_fixed;
    model->auto_center_x = center_x;
    model->auto_center_y = center_y;
    model->auto_center_z = center_z;
    model->auto_scaled = 1;
    model->auto_centered = center_flag;

    // Update bounding sphere parameters according to applied scaling and set distance from sphere (fast, O(1))
    {
        float scale_f32 = FIXED_TO_FLOAT(scale_fixed);
        if (model->bs_valid) {
            if (center_flag) {
                model->bs_cx = (model->bs_cx - FIXED_TO_FLOAT(center_x)) * scale_f32;
                model->bs_cy = (model->bs_cy - FIXED_TO_FLOAT(center_y)) * scale_f32;
                model->bs_cz = (model->bs_cz - FIXED_TO_FLOAT(center_z)) * scale_f32;
            } else {
                model->bs_cx = model->bs_cx * scale_f32;
                model->bs_cy = model->bs_cy * scale_f32;
                model->bs_cz = model->bs_cz * scale_f32;
            }
            model->bs_r *= scale_f32;
            model->bs_valid = 1;
            params->distance = computeDistanceFromBoundingSphere(model, margin);
        } else {
            // fallback (rare): compute from vertices (slower)
            params->distance = computeDistanceToFit(&model->vertices, margin);
        }
    }
    printf("[FIT] p%g radius=%.3f scale=%.4f applied. New distance: %.2f\n", percentile*100.0, (double)Rf, (double)scale_f, FIXED_TO_FLOAT(params->distance));
}



// Helper macro to swap face indices in the sorted_face_indices array
// (We swap indices, not the faces themselves, to keep the buffer intact)
#define SWAP_FACE(faces, i, j) \
    do { \
        int temp_idx = faces->sorted_face_indices[i]; \
        faces->sorted_face_indices[i] = faces->sorted_face_indices[j]; \
        faces->sorted_face_indices[j] = temp_idx; \
    } while (0)


// Function to draw polygons with QuickDraw
void drawPolygons(Model3D* model, int* vertex_count, int face_count, int vertex_count_total) {
    int i, j;
    VertexArrays3D* vtx = &model->vertices;
    FaceArrays3D* faces = &model->faces;
    Handle polyHandle;
    DynamicPolygon *poly;
    int min_x, max_x, min_y, max_y;
    int valid_faces_drawn = 0;
    int invalid_faces_skipped = 0;
    int triangle_count = 0;
    int quad_count = 0;
    Pattern pat;
    
    // Precompute screen scale and screen bounds for culling
    int screenScale = mode / 320;
    int screenW = screenScale * (CENTRE_X * 2);
    int screenH = screenScale * (CENTRE_Y * 2);

    // Use global persistent handle to avoid repeated NewHandle/DisposeHandle
    // Each call allocates fresh if needed, but reuses same handle block
    if (globalPolyHandle == NULL) {
        int max_polySize = 2 + 8 + (4 * 4);  // Max for quad (4 vertices)
        globalPolyHandle = NewHandle((long)max_polySize, userid(), 0xC014, 0L);
        if (globalPolyHandle == NULL) {
            printf("Error: Unable to allocate global polygon handle\n");
            return;
        }
    }
    
    polyHandle = globalPolyHandle;
    
    // Make sure handle is unlocked before locking
    if (poly_handle_locked) {
        HUnlock(polyHandle);
        poly_handle_locked = 0;
    }
    HLock(polyHandle);
    poly_handle_locked = 1;

    SetPenMode(0);
    // printf("\nDrawing polygons on screen:\n");

    // Set fill pen once per frame (reduces state changes)
    SetSolidPenPat(14);

    // Use sorted_face_indices to draw in correct depth order
    // Draw ALL faces - painter's algorithm handles occlusion
    int start_face = 0;
    int max_faces_to_draw = face_count;
    for (i = start_face; i < start_face + max_faces_to_draw; i++) {
        int face_id = faces->sorted_face_indices[i];
        if (faces->display_flag[face_id] == 0) continue;
        if (faces->vertex_count[face_id] >= 3) {
            int offset = faces->vertex_indices_ptr[face_id];
            int vcount_face = faces->vertex_count[face_id];
            int *indices_base = &faces->vertex_indices_buffer[offset];

            // Quick validity pass: ensure all indices are valid to avoid undefined points
            int all_valid = 1;
            for (j = 0; j < vcount_face; ++j) {
                int vi = indices_base[j] - 1;
                if (vi < 0 || vi >= vtx->vertex_count) { all_valid = 0; break; }
            }
            if (!all_valid) { invalid_faces_skipped++; continue; }

            // Calculate polySize for this specific face
            int polySize = 2 + 8 + (vcount_face * 4);
            poly = (DynamicPolygon *)*polyHandle;
            poly->polySize = polySize;

            // Cache arrays
            int *x2d = vtx->x2d;
            int *y2d = vtx->y2d;

            // Initialize min/max with the first vertex to avoid sentinel checks
            int first_vi = indices_base[0] - 1;
            int x = x2d[first_vi];
            int y = y2d[first_vi];
            poly->polyPoints[0].h = screenScale * x;
            poly->polyPoints[0].v = y;
            min_x = max_x = x;
            min_y = max_y = y;

            // Fill remaining points using cached data
            for (j = 1; j < vcount_face; ++j) {
                int vi = indices_base[j] - 1;
                x = x2d[vi];
                y = y2d[vi];
                poly->polyPoints[j].h = screenScale * x;
                poly->polyPoints[j].v = y;
                if (x < min_x) min_x = x;
                if (x > max_x) max_x = x;
                if (y < min_y) min_y = y;
                if (y > max_y) max_y = y;
            }

            poly->polyBBox.h1 = min_x;
            poly->polyBBox.v1 = min_y;
            poly->polyBBox.h2 = max_x;
            poly->polyBBox.v2 = max_y;

            // Bounding-box culling (convert to screen pixels)
            int sc_min_x = screenScale * min_x;
            int sc_max_x = screenScale * max_x;
            int sc_min_y = screenScale * min_y;
            int sc_max_y = screenScale * max_y;
            if (sc_max_x < 0 || sc_min_x >= screenW || sc_max_y < 0 || sc_min_y >= screenH) {
                // Off-screen; skip drawing
            } else {
                if (framePolyOnly) {
                    // Frame-only rendering (no fill)
                    SetSolidPenPat(7);
                    FramePoly(polyHandle);
                    SetSolidPenPat(14); // keep fill pen as default for next faces
                } else {
                    GetPenPat(pat);
                    FillPoly(polyHandle, pat);
                    SetSolidPenPat(7);
                    FramePoly(polyHandle);
                    SetSolidPenPat(14); // restore fill pen
                }
                valid_faces_drawn++;
                if (vcount_face == 3) triangle_count++;
                else if (vcount_face == 4) quad_count++;
            }
        } else {
            invalid_faces_skipped++;
        }
    }
    // Print statistics after drawing
//     printf("Display statistics: %d valid faces drawn, %d invalid faces skipped\n", valid_faces_drawn, invalid_faces_skipped);
//     printf("Triangles: %d, Quads: %d\n", triangle_count, quad_count);
}

void DoColor() {
        Rect r;
        unsigned char pstr[4];  // Pascal string: [length][characters...]]

        SetRect (&r, 0, 1, mode / 320 *10, 11);
        for (int i = 0; i < 16; i++) {
            SetSolidPenPat(i);
            PaintRect(&r);

            if (i == 0) {
                SetSolidPenPat(15); // White frame for black background
                FrameRect(&r);
            }

            MoveTo(r.h1, r.v2+10);
            // Create a Pascal string to display the number
            if (i < 10) {
                pstr[0] = 1;           // Length: 1 character
                pstr[1] = '0' + i;     // Digit 0-9
            } else {
                pstr[0] = 2;           // Length: 2 characters
                pstr[1] = '0' + (i / 10);      // Tens (1 for 10-15)
                pstr[2] = '0' + (i % 10);      // Units (0-5 for 10-15)
            }
            DrawString(pstr);
            OffsetRect(&r, 20, 0);
        }
}

void DoText() {
        shroff();
        putchar((char) 12); // Clear screen    
}


// ==============================================================
// THIS IS THE MAIN PROGRAM
// ==============================================================
//
    int main() {
        Model3D* model;
        ObserverParams params;
        char filename[100];
        char input[50];
        int colorpalette = 0; // default color palette


    newmodel:
        printf("===================================\n");
        printf("       3D OBJ file viewer\n");
        printf("===================================\n\n");
        printf("A tribute to Robert DONY\n");
        printf("Author of \"Calcul des parties cachees\" (Masson, 1986)\n\n");

        // Creer le modele 3D
        model = createModel3D();
        if (model == NULL) {
            printf("Error: Unable to allocate memory for 3D model\n");
            printf("Press any key to quit...\n");
            keypress();
            return 1;
        }

        // Ask for filename (loop until a non-empty filename is entered and the model loads)
        while (1) {
            printf("Enter the filename to read (ENTER to exit): ");
            if (fgets(filename, sizeof(filename), stdin) != NULL) {
                size_t len = strlen(filename);
                if (len > 0 && filename[len-1] == '\n') {
                    filename[len-1] = '\0';
                }
            } else {
                // EOF or input error - exit gracefully
                printf("\nInput error or EOF. Exiting.\n");
                destroyModel3D(model);
                return 1;
            }

            if (filename[0] == '\0') {
                printf("No filename entered. Exiting.\n");
                destroyModel3D(model);
                return 0; // Exit program when user presses ENTER with empty filename
            }

            // Try to load the model; if it fails, inform the user, reset model state, and re-prompt
            if (loadModel3D(model, filename) < 0) {
                // printf("\nError loading file '%s'. Please try again.\n", filename);
                // Destroy and recreate model to ensure clean state for next attempt
                destroyModel3D(model);
                model = createModel3D();
                if (model == NULL) {
                    printf("Error: Unable to allocate memory for 3D model after failed load. Exiting.\n");
                    printf("Press any key to quit...\n");
                    keypress();
                    return 1;
                }
                continue;
            }
            // Successfully loaded
            break;
        }

        // Get observer parameters
        getObserverParams(&params, model);


    bigloop:
        // Process model with parameters - OPTIMIZED VERSION
            // Process model with parameters - OPTIMIZED VERSION
        printf("Processing model...\n");
        if (framePolyOnly) {
            // Wireframe mode: only project vertices and set simple face visibility—skip face sorting
            processModelWireframe(model, &params, filename);
        } else {
            processModelFast(model, &params, filename);
        }


    loopReDraw:
        {
            int key = 0;
            char input[50];

            if (model->faces.face_count > 0) {
                // Initialize QuickDraw
                startgraph(mode);
                // Draw 3D object
                drawPolygons(model, model->faces.vertex_count, model->faces.face_count, model->vertices.vertex_count);
                // display available colors
                if (colorpalette == 1) { 
                    DoColor(); 
                }

                // Wait for key press and get key code
        asm 
            {
            sep #0x20
        loop:
            lda >0xC000     // Read the keyboard status from memory address 0xC000
            bpl loop        // Wait until no key is pressed (= until bit 7 on)
            and #0x007f     // Clear the high bit
            sta >0xC010     // Clear the keypress by writing to 0xC010
            sta key         // Store the key code in variable 'key'
            rep #0x30
            }

        endgraph();        // Close QuickDraw
        }

        DoText();           // Show text screen

    #if ENABLE_DEBUG_SAVE
        sprintf(input, "You pressed key code: %d\n", key);
        printf("%s", input);
    #endif

        // Handle keyboard input with switch statement
        switch (key) {
            case 32:  // Space bar - display info and redraw
                printf("===================================\n");
                printf(" Model information and parameters\n");
                printf("===================================\n");
                printf("Model: %s\n", filename);
                printf("Vertices: %d, Faces: %d\n", model->vertices.vertex_count, model->faces.face_count);
                printf("Observer Parameters:\n");
                printf("    Distance: %.2f\n", FIXED_TO_FLOAT(params.distance));
                printf("    Horizontal Angle: %d deg\n", params.angle_h);
                printf("    Vertical Angle: %d deg\n", params.angle_v);
                printf("    Screen Rotation Angle: %d deg\n", params.angle_w);
                if (model->auto_scaled) {
                    printf("    Auto-scale: ON (factor %.4f, centered: %s)\n", FIXED_TO_FLOAT(model->auto_scale), model->auto_centered ? "yes" : "no");
                } else {
                    printf("    Auto-scale: OFF\n");
                }
                printf("===================================\n");
                printf("\n");
                printf("Press any key to continue...\n");
                keypress();
                goto loopReDraw;

            case 82:  // 'R' - revert auto-scale
            case 114: // 'r'
                if (model != NULL && model->auto_scaled) {
                    revertAutoScaleModel(model);
                    printf("Auto-scale reverted.\n");
                } else {
                    printf("No auto-scale to revert.\n");
                }
                goto bigloop;

            case 43:  // '+' - ensure auto-fit then increase distance by 10%
            case 61:  // '=' also acts as '+' on some keyboards
                if (model != NULL) {
                    if (!model->auto_scaled) {
                        float default_target = 200.0f;
                        fitModelToView(model, &params, default_target, 0.92f, 0.99f, 1);
                        printf("Auto-fit applied (factor %.4f). Press 'r' to revert or +/- to adjust distance.\n", FIXED_TO_FLOAT(model->auto_scale));
                    }
                    params.distance = params.distance + (params.distance / 10);
                    printf("Distance increased -> %.2f\n", FIXED_TO_FLOAT(params.distance));
                } else {
                    printf("No model loaded.\n");
                }
                goto bigloop;

            case 45:  // '-' - ensure auto-fit then decrease distance by 10%
                if (model != NULL) {
                    if (!model->auto_scaled) {
                        float default_target = 200.0f;
                        fitModelToView(model, &params, default_target, 0.92f, 0.99f, 1);
                        printf("Auto-fit applied (factor %.4f). Press 'r' to revert or +/- to adjust distance.\n", FIXED_TO_FLOAT(model->auto_scale));
                    }
                    params.distance = params.distance - (params.distance / 10);
                    printf("Distance decreased -> %.2f\n", FIXED_TO_FLOAT(params.distance));
                } else {
                    printf("No model loaded.\n");
                }
                goto bigloop;

            case 65:  // 'A' - decrease distance
            case 97:  // 'a'
                params.distance = params.distance - (params.distance / 10);
                printf("Distance decreased -> %.2f\n", FIXED_TO_FLOAT(params.distance));
                goto bigloop;

            case 90:  // 'Z' - increase distance  
            case 122: // 'z'
                params.distance = params.distance + (params.distance / 10);
                printf("Distance increased -> %.2f\n", FIXED_TO_FLOAT(params.distance));
                goto bigloop;

            case 21:  // Right arrow - increase horizontal angle
                params.angle_h = normalize_deg(params.angle_h + 10);
                goto bigloop;

            case 8:   // Left arrow - decrease horizontal angle
                params.angle_h = normalize_deg(params.angle_h - 10);
                goto bigloop;

            case 10:  // Down arrow - decrease vertical angle
                params.angle_v = normalize_deg(params.angle_v - 10);
                goto bigloop;

            case 11:  // Up arrow - increase vertical angle
                params.angle_v = normalize_deg(params.angle_v + 10);
                goto bigloop;

            case 87:  // 'W' - increase screen rotation angle
            case 119: // 'w'
                params.angle_w = normalize_deg(params.angle_w + 10);
                goto bigloop;

            case 88:  // 'X' - decrease screen rotation angle
            case 120: // 'x'
                params.angle_w = normalize_deg(params.angle_w - 10);
                goto bigloop;
        
            case 67:  // 'C' - toggle color palette display
            case 99:  // 'c'
                colorpalette ^= 1; // Toggle between 0 and 1
                goto loopReDraw;

case 70:  // 'F' - toggle fast/normal painter
case 102: // 'f'
    painterFastMode ^= 1;
    printf("Fast painter: %s\n", painterFastMode ? "ON (tests 1-3 only)" : "OFF (full tests)");
    if (model != NULL) {
        printf("Reprocessing model with current mode...\n");
        processModelFast(model, &params, filename);
    }
    goto loopReDraw;

case 80:  // 'P' - toggle frame-only polygon rendering (was 'F')
case 112: // 'p'
                framePolyOnly ^= 1;
                printf("Frame-only polygons: %s\n", framePolyOnly ? "ON" : "OFF");
                if (!framePolyOnly && model != NULL) {
                    // Switched back to filled polygons — re-run full processing to recompute depths & ordering
                    printf("Switching to filled mode: reprocessing model (sorting faces)...\n");
                    processModelFast(model, &params, filename);
                }
                goto loopReDraw;

            case 69:  // 'E' - dump face equations to equ.csv
            case 101: // 'e'
                if (model != NULL) {
                    dumpFaceEquationsCSV(model, "equ.csv");
                }
                goto loopReDraw;

            case 78:  // 'N' - load new model
            case 110: // 'n'
                destroyModel3D(model);
                goto newmodel;

            case 75:  // 'K' - edit angles/distance (no reload; ENTER may auto-fit)
            case 107: // 'k'
                getObserverParams(&params, model);
                printf("Observer parameters updated.\n");
                goto bigloop;
        
            // dispaly help
            case 72:  // 'H'
            case 104: // 'h'
                printf("===================================\n");
                printf("    HELP - Keyboard Controller\n");
                printf("===================================\n\n");
                printf("Space: Display model info\n");
                printf("A/Z: Increase/Decrease distance\n");
                printf("+/-: Apply auto-fit if none, then increase/decrease distance (use 'r' to revert auto-fit)\n");
                printf("K: Edit angles/distance (ENTER may trigger auto-fit)\n");
                printf("Arrow Left/Right: Decrease/Increase horizontal angle\n");
                printf("Arrow Up/Down: Increase/Decrease vertical angle\n");
                printf("W/X: Increase/Decrease screen rotation angle\n");
                printf("C: Toggle color palette display\n");
                printf("F: Toggle fast painter (default: ON — tests 1-3 only)\n");
                printf("P: Toggle frame-only polygons (default: OFF)\n");
                printf("E: Dump face equations to equ.csv (debug)\n");
                printf("N: Load new model\n");
                printf("H: Display this help message\n");
                printf("ESC: Quit program\n");
                printf("===================================\n");
                printf("\n");
                printf("Press any key to continue...\n");
                keypress();
                goto loopReDraw;

            case 27:  // ESC - quit
                goto end;
            
            default:  // All other keys - redraw
                goto loopReDraw;
        }
        }  // End of loopReDraw block

        end:
        // Cleanup and exit
        // Dispose of the global polygon handle if it was allocated
        if (globalPolyHandle != NULL) {
            if (poly_handle_locked) {
                HUnlock(globalPolyHandle);
            }
            DisposeHandle(globalPolyHandle);
            globalPolyHandle = NULL;
        }
        destroyModel3D(model);
        return 0;
    }