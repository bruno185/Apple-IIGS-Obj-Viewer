#!/usr/bin/env python3
import math
from pascal_sim import parse_obj

OBJ='Windows\\3d objects\\q2.obj'
H=30; V=20; W=0; D=39.214455

verts, faces = parse_obj(OBJ)

# C-like float transform
cos_h = math.cos(math.radians(H)); sin_h = math.sin(math.radians(H))
cos_v = math.cos(math.radians(V)); sin_v = math.sin(math.radians(V))
cos_h_cos_v = cos_h * cos_v
sin_h_cos_v = sin_h * cos_v
cos_h_sin_v = cos_h * sin_v
sin_h_sin_v = sin_h * sin_v

distance = D

obs = []
for (x,y,z) in verts:
    term1 = x * cos_h_cos_v
    term2 = y * sin_h_cos_v
    term3 = z * sin_v
    zo = -term1 - term2 - term3 + distance
    if zo > 0.0:
        xo = -x * sin_h + y * cos_h
        yo = -x * cos_h_sin_v - y * sin_h_sin_v + z * cos_v
    else:
        xo = 0.0; yo = 0.0
    obs.append((xo,yo,zo))

# compute Newell-based plane coeffs using obs coords
planes = []
for fidx, idxs in enumerate(faces):
    if len(idxs) < 3:
        planes.append((0,0,0,0)); continue
    i0, i1, i2 = idxs[0], idxs[1], idxs[2]
    x1,y1,z1 = obs[i0]
    x2,y2,z2 = obs[i1]
    x3,y3,z3 = obs[i2]
    a = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
    b = -x1*(z2 - z3) + x2*(z1 - z3) - x3*(z1 - z2)
    c = x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2)
    t1 = y2*z3 - y3*z2; t2 = y1*z3 - y3*z1; t3 = y1*z2 - y2*z1
    d = -x1 * t1 + x2 * t2 - x3 * t3
    planes.append((a,b,c,d))

for i,p in enumerate(planes):
    print('face',i,'simC a,b,c,d =',p)

print('\nNow showing pair_debug planes from C:')
with open('Windows\\build_vs\\bin\\Release\\pair_debug_v4.csv','r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split(',')
        if parts[0]=='0' and parts[1]=='1':
            print('C pair 0,1 planes from CSV', parts[4:12])
            break
