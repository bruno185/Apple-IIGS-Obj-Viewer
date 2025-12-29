from pascal_sim import parse_obj, transform_and_project, compute_faces
import csv, os
OBJ='Windows\\3d objects\\q2.obj'
PVF='Windows\\build_vs\\bin\\Release\\pair_values_v4.csv'
# Find D from equ_windows if present
D = None
tmp = os.environ.get('TEMP')
if tmp:
    fn = os.path.join(tmp, 'equ_windows.csv')
    if os.path.exists(fn):
        with open(fn,'r') as f:
            first = f.readline()
            if 'distance=' in first:
                try:
                    part = first.split('distance=')[-1]
                    D = float(part)
                except:
                    D=None
if D is None: D = 39.214455

verts, faces = parse_obj(OBJ)
# apply loader Y/Z swap used earlier
verts_c = [(x,z,y) for (x,y,z) in verts]
H,V,W = 30,20,0
xo,yo,zo,x2d,y2d = transform_and_project(verts_c, H, V, W, D)

# Pascal fixed ints per face
pascal_fixed = {}
for fi, idxs in enumerate(faces):
    if len(idxs) < 3: continue
    i0,i1,i2 = idxs[0], idxs[1], idxs[2]
    x1f,y1f,z1f = xo[i0], yo[i0], zo[i0]
    x2f,y2f,z2f = xo[i1], yo[i1], zo[i1]
    x3f,y3f,z3f = xo[i2], yo[i2], zo[i2]
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
    pascal_fixed[fi] = {'a_f':a_f,'b_f':b_f,'c_f':c_f,'d_f':d_f}

# Read PVF logged a_f_int etc
c_fixed = {}
with open(PVF,'r') as f:
    for line in f:
        if ',plane_raw,' in line or ',plane_fixed,' in line or ',plane_arr,' in line:
            parts = line.strip().split(',')
            # parts[4] is fidx in our format
            try:
                fidx = int(parts[4])
            except:
                continue
            key = parts[9]
            try:
                val = float(parts[-1])
            except:
                continue
            if fidx not in c_fixed: c_fixed[fidx] = {}
            c_fixed[fidx][key] = val

# Compare
for fi in range(len(faces)):
    pf = pascal_fixed.get(fi)
    cf = c_fixed.get(fi, {})
    if not pf:
        continue
    print('Face',fi)
    # compare int values
    for k in ['a_f','b_f','c_f','d_f']:
        pas = pf[k]
        ckey = k if k in cf else (k if k in cf else None)
        cval = int(cf.get(k, cf.get(k+'_int', 0)))
        equal = pas == cval
        print(' ',k, 'pas=',pas, 'c=',cval, 'match=', equal)
    # also compare float arrondi values
    a_f = pf['a_f']/65536.0; b_f = pf['b_f']/65536.0; c_f = pf['c_f']/65536.0; d_f = pf['d_f']/65536.0
    a_arr = cf.get('a_arr', None); b_arr = cf.get('b_arr', None); c_arr = cf.get('c_arr', None); d_arr = cf.get('d_arr', None)
    print('  floats pas=',a_f,b_f,c_f,d_f)
    print('  floats c =',a_arr,b_arr,c_arr,d_arr)
    print()