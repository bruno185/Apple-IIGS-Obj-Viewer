#!/usr/bin/env python3
import math

# Read obj
verts=[]
faces=[]
with open('Windows/3d objects/q2.obj','r') as f:
    for line in f:
        line=line.strip()
        if not line or line.startswith('#'): continue
        parts=line.split()
        if parts[0]=='v': verts.append((float(parts[1]),float(parts[2]),float(parts[3])))
        elif parts[0]=='f': idxs=[int(p.split('/')[0])-1 for p in parts[1:]]; faces.append(idxs)

# C compute_obs_vertices
H=30; V=20; D=39.214455
cos_h=math.cos(math.radians(H)); sin_h=math.sin(math.radians(H))
cos_v=math.cos(math.radians(V)); sin_v=math.sin(math.radians(V))
cos_h_cos_v = cos_h * cos_v
sin_h_cos_v = sin_h * cos_v
cos_h_sin_v = cos_h * sin_v
sin_h_sin_v = sin_h * sin_v

obs=[]
for (x,y,z) in verts:
    term1 = x * cos_h_cos_v
    term2 = y * sin_h_cos_v
    term3 = z * sin_v
    zo = -term1 - term2 - term3 + D
    if zo > 0:
        xo = -x * sin_h + y * cos_h
        yo = -x * cos_h_sin_v - y * sin_h_sin_v + z * cos_v
    else:
        xo=0.0; yo=0.0
    obs.append((xo,yo,zo))

# Compute planes using formulas
for fi,face in enumerate(faces):
    i0,i1,i2 = face[0],face[1],face[2]
    x1,y1,z1 = obs[i0]
    x2,y2,z2 = obs[i1]
    x3,y3,z3 = obs[i2]
    a = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
    b = -x1*(z2 - z3) + x2*(z1 - z3) - x3*(z1 - z2)
    c = x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2)
    t1 = y2*z3 - y3*z2; t2 = y1*z3 - y3*z1; t3 = y1*z2 - y2*z1
    d = -x1 * t1 + x2 * t2 - x3 * t3
    print(f"Face {fi}: a={a:.6f}, b={b:.6f}, c={c:.6f}, d={d:.6f}")

# Use pascal_sim compute_faces to compare
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float
xo,yo,zo,x2d,y2d = transform_and_project(verts,H,V,0,D)
faces_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)
for fi in range(len(faces_info)):
    a = to_float(faces_info[fi]['a']); b = to_float(faces_info[fi]['b']); c = to_float(faces_info[fi]['c']); d = to_float(faces_info[fi]['d'])
    print(f"Pascal Face {fi}: a={a:.6f}, b={b:.6f}, c={c:.6f}, d={d:.6f}")
