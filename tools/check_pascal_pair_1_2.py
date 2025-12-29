from pascal_sim import parse_obj, transform_and_project, compute_faces, test_pair
path='Windows\\3d objects\\q2.obj'
H=30;V=20;W=0;D=39.214455
verts,faces=parse_obj(path)
xo,yo,zo,x2d,y2d=transform_and_project(verts,H,V,W,D)
import pascal_sim as ps
ps.xo=xo; ps.yo=yo; ps.zo=zo
faces_info=compute_faces(verts,faces,xo,yo,zo,x2d,y2d)
print('Pair 1 2 ->', test_pair(faces_info,1,2))
