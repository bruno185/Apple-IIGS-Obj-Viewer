from pascal_sim import parse_obj, transform_and_project, compute_faces, to_fixed
import csv

OBJ='Windows\\3d objects\\q2.obj'
PAIR='Windows\\build_vs\\bin\\Release\\pair_values_v4.csv'
FACE_IDX=1

# Pascal fixed pipeline (use D = s_distance = 300)
verts, faces = parse_obj(OBJ)
# apply C loader Y/Z swap used earlier in compare_obs_pascal
verts_c = [(x, z, y) for (x,y,z) in verts]
H,V,W,D = 30,20,0,39.214455
xo,yo,zo,x2d,y2d = transform_and_project(verts_c, H, V, W, D)
face_info = compute_faces(verts_c, faces, xo, yo, zo, x2d, y2d)
# extract pascal fixed ints for face
fi = faces[FACE_IDX]
idx0,idx1,idx2 = fi[0], fi[1], fi[2]
x1f = xo[idx0]; y1f = yo[idx0]; z1f = zo[idx0]
x2f = xo[idx1]; y2f = yo[idx1]; z2f = zo[idx1]
x3f = xo[idx2]; y3f = yo[idx2]; z3f = zo[idx2]
# compute terms as pascal_sim does
term1 = (y1f * (z2f - z3f)) >> 16
term2 = (y2f * (z3f - z1f)) >> 16
term3 = (y3f * (z1f - z2f)) >> 16
a_f = term1 + term2 + term3
b_f = -((x1f * (z2f - z3f)) >> 16) + ((x2f * (z1f - z3f)) >> 16) - ((x3f * (z1f - z2f)) >> 16)
c_f = ((x1f * (y2f - y3f)) >> 16) - ((x2f * (y1f - y3f)) >> 16) + ((x3f * (y1f - y2f)) >> 16)
t1f = ((y2f * z3f) >> 16) - ((y3f * z2f) >> 16)
t2f = ((y1f * z3f) >> 16) - ((y3f * z1f) >> 16)
t3f = ((y1f * z2f) >> 16) - ((y2f * z1f) >> 16)
d_f = -((x1f * t1f) >> 16) + ((x2f * t2f) >> 16) - ((x3f * t3f) >> 16)
print('Pascal fixed ints:')
print(' x1f,y1f,z1f =', x1f,y1f,z1f)
print(' x2f,y2f,z2f =', x2f,y2f,z2f)
print(' x3f,y3f,z3f =', x3f,y3f,z3f)
print(' term1,term2,term3 =', term1,term2,term3)
print(' a_f,b_f,c_f,d_f =', a_f,b_f,c_f,d_f)
print(' a,b,c,d floats =', a_f/65536.0, b_f/65536.0, c_f/65536.0, d_f/65536.0)

# Read C logged fixed ints for same face from PVF
c_vals = {}
with open(PAIR,'r') as f:
    for line in f:
        if ',plane_fixed' in line or ',plane_raw' in line or ',plane_arr' in line:
            parts = line.strip().split(',')
            # Look for face index and name
            try:
                fidx = int(parts[4])
            except:
                continue
            if fidx != FACE_IDX: continue
            # parts like: instr,plane_fixed,0,0,0,1,0,plane,-1,term1, , -207211.000000000
            expr = parts[9]
            val = float(parts[-1])
            c_vals[expr] = val

print('\nC logged values (from PVF):')
for key in sorted(c_vals.keys()):
    print(key, c_vals[key])

# Compare key ints
print('\nCompare Pascal ints vs C a_f_int,b_f_int,c_f_int,d_f_int if present:')
print('pas a_f=', a_f, 'c a_f_int=', int(c_vals.get('a_f', 0)))
print('pas b_f=', b_f, 'c b_f_int=', int(c_vals.get('b_f', 0)))
print('pas c_f=', c_f, 'c c_f_int=', int(c_vals.get('c_f', 0)))
print('pas d_f=', d_f, 'c d_f_int=', int(c_vals.get('d_f', 0)))
