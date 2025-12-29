from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float
OBJ='Windows\\3d objects\\q2.obj'
H=30; V=20; W=0; D=39.214455
verts, faces = parse_obj(OBJ)
xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
faces_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)
fi = 1
f = faces_info[fi]
print('face', fi, 'display', f['display'])
print('a_fixed', f['a'])
print('a_float', to_float(f['a']))
print('arr(a_float)', 0.0 if abs(to_float(f['a']))<0.01 else to_float(f['a']))
print('vertices idxs', f['idxs'])
for idx in f['idxs']:
    print(idx, 'xo,yo,zo (fixed)', xo[idx], yo[idx], zo[idx], '-> float', to_float(xo[idx]), to_float(yo[idx]), to_float(zo[idx]))
