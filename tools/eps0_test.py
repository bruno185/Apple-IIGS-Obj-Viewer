#!/usr/bin/env python3
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_fixed, to_float
import epsilon_sweep

path='cube.obj'
H=30; V=20; W=0; D=180.0
verts, faces = parse_obj(path)
xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
face_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)
eps_fixed = to_fixed(0.0)
print('eps=0 test:')
pa = epsilon_sweep.simulate_painter(face_info, xo, yo, zo, eps_fixed)
pi = epsilon_sweep.simulate_painter_insert(face_info, xo, yo, zo, eps_fixed)
pc = epsilon_sweep.full_pascal_order(face_info, xo, yo, zo, eps_fixed)
print(' painter_adj', pa)
print(' painter_insert', pi)
print(' pascal', pc)
for p in [(4,0),(4,1),(5,0),(5,1),(0,1)]:
    res = epsilon_sweep.test_pair_eps(face_info, xo, yo, zo, p[0], p[1], eps_fixed)
    print(f' pair {p}: {res}')
