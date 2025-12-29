from pascal_sim import parse_obj, transform_and_project, compute_faces
path='Windows\\3d objects\\q2.obj'
H=30;V=20;W=0;D=300.0
verts,faces=parse_obj(path)
xo,yo,zo,x2d,y2d=transform_and_project(verts,H,V,W,D)
faces_info=compute_faces(verts,faces,xo,yo,zo,x2d,y2d)

def arr(r):
    return 0.0 if abs(r) < 0.01 else r

for i,fi in enumerate(faces_info):
    print(i, arr(fi['a']), arr(fi['b']), arr(fi['c']), arr(fi['d']))
