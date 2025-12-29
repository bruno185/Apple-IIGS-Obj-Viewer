from pascal_sim import parse_obj, transform_and_project, compute_faces, test_pair
path='Windows\\3d objects\\q2.obj'
H=30;V=20;W=0;D=39.214455
verts,faces=parse_obj(path)
xo,yo,zo,x2d,y2d=transform_and_project(verts,H,V,W,D)
# export xo/yo/zo into pascal_sim module globals so test_pair can access them
import pascal_sim as ps
ps.xo = xo; ps.yo = yo; ps.zo = zo
faces_info=compute_faces(verts,faces,xo,yo,zo,x2d,y2d)
for a,b in [(0,1),(0,2),(1,2)]:
    try:
        res = test_pair(faces_info,a,b)
        print('Pair',a,b,'->', res)
    except Exception as e:
        print('Pair',a,b,'-> ERROR', e)
        # dump some values to help debug
        f1=faces_info[a]; f2=faces_info[b]
        print(' f1.a,b,c,d =', f1['a'], f1['b'], f1['c'], f1['d'])
        print(' f2.a,b,c,d =', f2['a'], f2['b'], f2['c'], f2['d'])
        print(' xo sample:', xo[ f2['idxs'][0] ], xo[ f1['idxs'][0] ])
        raise
