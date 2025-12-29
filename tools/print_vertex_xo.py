from pascal_sim import parse_obj, transform_and_project
verts, faces = parse_obj('Windows\\3d objects\\q2.obj')
H=30; V=20; W=0; D=300.0
xo,yo,zo,x2d,y2d = transform_and_project(verts,H,V,W,D)
face = faces[1]
print('face 1 idxs',face)
for idx in face[:3]:
    print('idx', idx, 'xo,yo,zo fixed =', xo[idx], yo[idx], zo[idx], '-> floats', xo[idx]/65536.0, yo[idx]/65536.0, zo[idx]/65536.0)
