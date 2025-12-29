#!/usr/bin/env python3
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float

def arr(r):
    return 0.0 if abs(r) < 0.01 else r
OBJ='Windows\\3d objects\\q2.obj'
H=30;V=20;W=0;D=39.214455
verts,faces=parse_obj(OBJ)
xo,yo,zo,x2d,y2d = transform_and_project(verts,H,V,W,D)
fi = compute_faces(verts,faces,xo,yo,zo,x2d,y2d)
for idx in range(len(fi)):
    print('face',idx,'py arr a,b,c,d =', arr(to_float(fi[idx]['a'])), arr(to_float(fi[idx]['b'])), arr(to_float(fi[idx]['c'])), arr(to_float(fi[idx]['d'])))
