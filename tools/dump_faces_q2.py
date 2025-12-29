from pascal_sim import parse_obj, transform_and_project, compute_faces
path='Windows\\3d objects\\q2.obj'
H=30;V=20;W=0;D=300.0
verts,faces=parse_obj(path)
print('faces count',len(faces))
for i,f in enumerate(faces): print(i,f)
xo,yo,zo,x2d,y2d=transform_and_project(verts,H,V,W,D)
faces_info=compute_faces(verts,faces,xo,yo,zo,x2d,y2d)
for i in range(len(faces_info)):
    fi = faces_info[i]
    print('face',i,'idxs',fi['idxs'],'a,b,c,d=',fi['a'],fi['b'],fi['c'],fi['d'])
    
print('\nFull faces_info dump:')
import pprint
pprint.pprint(faces_info)

