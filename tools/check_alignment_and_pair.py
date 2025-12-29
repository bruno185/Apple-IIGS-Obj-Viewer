from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float, test_pair
import csv

OBJ='Windows\\3d objects\\q2.obj'
PAIR_CSV='Windows\\build_vs\\bin\\Release\\pair_values_v4.csv'

# Pascal transform
verts, faces = parse_obj(OBJ)
H=30; V=20; W=0; D=39.214455
xo_p, yo_p, zo_p, x2d, y2d = transform_and_project(verts, H, V, W, D)
face_info = compute_faces(verts, faces, xo_p, yo_p, zo_p, x2d, y2d)

# Read C obs_dump values from pair_values_v4.csv
obs_c = {}
with open(PAIR_CSV, 'r', newline='') as f:
    for line in f:
        if not line.startswith('instr,obs_dump'): continue
        parts = line.strip().split(',')
        vid = int(parts[8]); field = parts[9]; val = float(parts[-1]); obs_c.setdefault(vid, {})[field] = val

# Compute centroid difference for face 1 vertices
face_idx = 1
idxs = face_info[face_idx]['idxs']
px = [to_float(xo_p[v]) for v in idxs]
py = [to_float(yo_p[v]) for v in idxs]
pz = [to_float(zo_p[v]) for v in idxs]
cx_pas = sum(px)/len(px); cy_pas = sum(py)/len(py); cz_pas = sum(pz)/len(pz)
cx_c = sum(obs_c[v]['xo'] for v in idxs)/len(idxs)
cy_c = sum(obs_c[v]['yo'] for v in idxs)/len(idxs)
cz_c = sum(obs_c[v]['zo'] for v in idxs)/len(idxs)
print('centroid_pas',cx_pas,cy_pas,cz_pas)
print('centroid_C',cx_c,cy_c,cz_c)

# Apply translation to Pascal coords so centroid aligns to C centroid
tx = cx_c - cx_pas; ty = cy_c - cy_pas; tz = cz_c - cz_pas
print('translation to apply to Pascal X,Y,Z =',tx,ty,tz)
# apply translation to xo_p,yo_p,zo_p
xo_p2 = {i:to_float(xo_p[i])+tx for i in range(len(xo_p))}
yo_p2 = {i:to_float(yo_p[i])+ty for i in range(len(yo_p))}
zo_p2 = {i:to_float(zo_p[i])+tz for i in range(len(zo_p))}

# recompute face planes from shifted Pascal coords (float arithmetic approximate)
def compute_face_plane_from_shifted(idxs):
    x1,y1,z1 = xo_p2[idxs[0]], yo_p2[idxs[0]], zo_p2[idxs[0]]
    x2,y2,z2 = xo_p2[idxs[1]], yo_p2[idxs[1]], zo_p2[idxs[1]]
    x3,y3,z3 = xo_p2[idxs[2]], yo_p2[idxs[2]], zo_p2[idxs[2]]
    a = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
    b = -x1*(z2 - z3) + x2*(z1 - z3) - x3*(z1 - z2)
    c = x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2)
    t1 = y2*z3 - y3*z2; t2 = y1*z3 - y3*z1; t3 = y1*z2 - y2*z1
    d = -x1*t1 + x2*t2 - x3*t3
    return a,b,c,d

print('\nOriginal Pascal face 1 plane (float):', to_float(face_info[1]['a']), to_float(face_info[1]['b']), to_float(face_info[1]['c']), to_float(face_info[1]['d']))
print('C plane from trace (face 1):')
# extract C plane from pair_values_v4.csv
with open(PAIR_CSV,'r') as f:
    for line in f:
        if ',planes,' in line and ',1,0,' in line:
            parts=line.strip().split(',')
            print(parts[6:10])
            break

print('\nShifted Pascal plane:', compute_face_plane_from_shifted(idxs))
# Evaluate pair (1,0) using shifted Pascal planes
print('\nRunning pascal test_pair on shifted coords:')
# build new face_info2 with shifted coords
face_info2 = []
for fi, idxs in enumerate(faces):
    n = len(idxs)
    zmin = min(zo_p2[i] for i in idxs); zmax=max(zo_p2[i] for i in idxs); zmean=sum(zo_p2[i] for i in idxs)//n
    a,b,c,d = compute_face_plane_from_shifted(idxs)
    face_info2.append({'a':a,'b':b,'c':c,'d':d,'zmin':zmin,'zmax':zmax,'zmean':zmean,'idxs':idxs})
res = test_pair(face_info2, 1, 0)
print('test_pair on shifted Pascal for (1,0):', res)
