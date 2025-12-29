#!/usr/bin/env python3
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float
from pascal_sim import fixed_mul_64, to_fixed
from math import isclose
OBJ='Windows\\3d objects\\q2.obj'
D=39.214455
for factor in [1.0, 0.25, 4.0]:
    d = D*factor
    xo,yo,zo,x2d,y2d = transform_and_project(parse_obj(OBJ)[0],30,20,0,d)
    fi = compute_faces(parse_obj(OBJ)[0], parse_obj(OBJ)[1], xo,yo,zo,x2d,y2d)
    a0 = to_float(fi[0]['a']); a1 = to_float(fi[1]['a']);
    print(f'factor={factor}: face0.a={a0:.6f}, face1.a={a1:.6f}')
