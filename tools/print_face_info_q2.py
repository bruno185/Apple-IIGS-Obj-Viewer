from pascal_sim import parse_obj, transform_and_project, compute_faces
verts, faces = parse_obj('Windows\\3d objects\\q2.obj')
H=30; V=20; W=0; D=300.0
xo,yo,zo,x2d,y2d = transform_and_project(verts,H,V,W,D)
face_info = compute_faces(verts,faces,xo,yo,zo,x2d,y2d)
for i,fi in enumerate(face_info):
    if i<3:
        print('Face',i,'idxs',fi['idxs'],'a,b,c,d=',fi['a'],fi['b'],fi['c'],fi['d'])
