from pascal_sim import parse_obj, transform_and_project, compute_faces
path='Windows\\3d objects\\q2.obj'
H=30;V=20;W=0;D=300.0
verts,faces=parse_obj(path)
xo,yo,zo,x2d,y2d=transform_and_project(verts,H,V,W,D)
faces_info=compute_faces(verts,faces,xo,yo,zo,x2d,y2d)

def arrondi(r):
    return 0.0 if abs(r) < 0.01 else r

import math
xobs=arrondi(D * math.cos(math.radians(H)) * math.cos(math.radians(V)))
yobs=arrondi(D * math.sin(math.radians(H)) * math.cos(math.radians(V)))
zobs=arrondi(D * math.sin(math.radians(V)))

def eval_ab(a,b):
    Ap=arrondi(faces_info[a]['a']); Bp=arrondi(faces_info[a]['b']); Cp=arrondi(faces_info[a]['c']); Dp=arrondi(faces_info[a]['d'])
    pos=neg=0
    for idx in faces_info[b]['idxs']:
        tempo=arrondi(Ap*verts[idx][0] + Bp*verts[idx][1] + Cp*verts[idx][2] + Dp)
        if tempo>0: pos+=1
        elif tempo<0: neg+=1
    if pos>0 and neg>0: return 0
    cotepoint = 1 if pos>0 else -1
    position = 1 if (Ap*xobs + Bp*yobs + Cp*zobs + Dp) > 0.0 else -1
    return 1 if position == cotepoint else -1

pairs=[(0,1),(0,2),(1,2)]
for a,b in pairs:
    try:
        res_ab = eval_ab(a,b)
    except Exception as e:
        res_ab = f'ERROR: {e}'
    try:
        res_ba = eval_ab(b,a)
    except Exception as e:
        res_ba = f'ERROR: {e}'
    print('pair',a,b,'devant(a,b)=',res_ab,'devant(b,a)=',res_ba)

