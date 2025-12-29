from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float
import csv
import sys
import os

OBJ='Windows\\3d objects\\q2.obj'
PAIR_CSV='Windows\\build_vs\\bin\\Release\\pair_values_v4.csv'

# Attempt to read D from equ_windows.csv (meta line)
D = None
try:
    tmp = os.environ.get('TEMP')
    if tmp:
        fn = os.path.join(tmp, 'equ_windows.csv')
        if os.path.exists(fn):
            with open(fn, 'r') as f:
                first = f.readline().strip()
                if 'distance=' in first:
                    part = first.split('distance=')[-1]
                    try:
                        D = float(part)
                    except:
                        D = None
except Exception:
    D = None
if D is None:
    D = 39.214455

# Pascal transform
verts, faces = parse_obj(OBJ)
# Apply same Y/Z swap as C loader: GS3Dp stores x, y := OBJ z, z := OBJ y
verts_c = [(x, z, y) for (x,y,z) in verts]
H=30; V=20; W=0
xo_p, yo_p, zo_p, x2d, y2d = transform_and_project(verts_c, H, V, W, D)
face_info = compute_faces(verts_c, faces, xo_p, yo_p, zo_p, x2d, y2d)

# Read C obs_dump values from pair_values_v4.csv
obs_c = {} # vid -> {xo,yo,zo}
try:
    with open(PAIR_CSV, 'r', newline='') as f:
        for line in f:
            if not line.startswith('instr,obs_dump'): continue
            parts = line.strip().split(',')
            # instr,obs_dump,0,0,0,-1,0,obs,vid,field,,value
            if len(parts) < 12: continue
            vid = int(parts[8])
            field = parts[9]
            val = float(parts[-1])
            obs_c.setdefault(vid, {})
            # map 'xo' 'yo' 'zo'
            obs_c[vid][field] = val
except FileNotFoundError:
    print('C trace not found:', PAIR_CSV)
    sys.exit(1)

# Face index to inspect
face_idx = 1
face = face_info[face_idx]
if isinstance(face, dict):
    idxs = face['idxs']
else:
    # fallback if compute_faces returned a bare list
    idxs = face
print('Pascal face', face_idx, 'indices', idxs)

mismatches = []
print('\nComparing vertices:')
for vid in idxs:
    # Pascal fixed values -> float
    x_pas = to_float(xo_p[vid]); y_pas = to_float(yo_p[vid]); z_pas = to_float(zo_p[vid])
    c = obs_c.get(vid, {})
    xo_c = c.get('xo', None); yo_c = c.get('yo', None); zo_c = c.get('zo', None)
    print(f'vid {vid}: Pascal xo={x_pas:.12f} yo={y_pas:.12f} zo={z_pas:.12f} | C xo={xo_c} yo={yo_c} zo={zo_c}')
    if xo_c is None or yo_c is None or zo_c is None:
        mismatches.append((vid, 'missing'))
        continue
    dx = xo_c - x_pas; dy = yo_c - y_pas; dz = zo_c - z_pas
    if abs(dx) > 1e-6 or abs(dy) > 1e-6 or abs(dz) > 1e-6:
        mismatches.append((vid, dx, dy, dz))

# Also compare plane coefficients (Pascal fixed -> float) with C's planes row in pair_values_v4.csv
c_planes = None
try:
    with open(PAIR_CSV, 'r') as f:
        for line in f:
            if ',planes,' in line:
                parts = line.strip().split(',')
                # pass,i,j,f1,f2,planes,Ap,Bp,Cp,Dp,...
                if len(parts) >= 11:
                    passn,i,j,f1,f2,stage = parts[0:6]
                    Ap = float(parts[6]); Bp = float(parts[7]); Cp = float(parts[8]); Dp = float(parts[9])
                    c_planes = (Ap,Bp,Cp,Dp)
                    break
except Exception:
    c_planes = None

pf = face_info[face_idx]
pa = to_float(pf['a']); pb = to_float(pf['b']); pc = to_float(pf['c']); pd = to_float(pf['d'])
print('\nPascal face coeffs (float): a={:.12f} b={:.12f} c={:.12f} d={:.12f}'.format(pa,pb,pc,pd))
if c_planes:
    print('C planes (from trace): Ap={:.12f} Bp={:.12f} Cp={:.12f} Dp={:.12f}'.format(*c_planes))
    # compare diffs
    da = (c_planes[0]-pa); db = (c_planes[1]-pb); dc = (c_planes[2]-pc); dd = (c_planes[3]-pd)
    print('Diffs C-Pascal (a,b,c,d)=', da, db, dc, dd)

print('\nFindings:')
if not mismatches:
    print('All obs vertices match Pascal transform within tolerance.')
else:
    for m in mismatches:
        if m[1] == 'missing': print('VID', m[0], 'missing in C trace')
        else: print('VID', m[0], 'diffs dx,dy,dz =', m[1], m[2], m[3])

# exit non-zero if significant plane mismatch or vertex mismatches
plane_issue = False
if c_planes:
    if abs(c_planes[0]-pa) > 0.01 or abs(c_planes[1]-pb) > 0.01 or abs(c_planes[2]-pc) > 0.01 or abs(c_planes[3]-pd) > 0.01:
        plane_issue = True

sys.exit(0 if (not mismatches and not plane_issue) else 2)
