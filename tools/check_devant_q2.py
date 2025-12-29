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

print('observer a=',xobs,yobs,zobs)
for a in range(len(faces)):
    for b in range(a+1,len(faces)):
        Ap=arrondi(faces_info[a]['a']); Bp=arrondi(faces_info[a]['b']); Cp=arrondi(faces_info[a]['c']); Dp=arrondi(faces_info[a]['d'])
        pos=neg=0
        for idx in faces_info[b]['idxs']:
            tempo=arrondi(Ap*verts[idx][0] + Bp*verts[idx][1] + Cp*verts[idx][2] + Dp)
            if tempo>0: pos+=1
            elif tempo<0: neg+=1
        res_ab = 0
        if not (pos>0 and neg>0):
            cotepoint = 1 if pos>0 else -1
            position = 1 if (Ap*xobs + Bp*yobs + Cp*zobs + Dp) > 0.0 else -1
            res_ab = 1 if position == cotepoint else -1
        # symmetric
        Ap2=arrondi(faces_info[b]['a']); Bp2=arrondi(faces_info[b]['b']); Cp2=arrondi(faces_info[b]['c']); Dp2=arrondi(faces_info[b]['d'])
        pos2=neg2=0
        for idx in faces_info[a]['idxs']:
            tempo=arrondi(Ap2*verts[idx][0] + Bp2*verts[idx][1] + Cp2*verts[idx][2] + Dp2)
            if tempo>0: pos2+=1
            elif tempo<0: neg2+=1
        res_ba = 0
        if not (pos2>0 and neg2>0):
            cotepoint2 = 1 if pos2>0 else -1
            position2 = 1 if (Ap2*xobs + Bp2*yobs + Cp2*zobs + Dp2) > 0.0 else -1
            res_ba = 1 if position2 == cotepoint2 else -1
        print('pair',a,b,'Devant(a,b)=',res_ab,'Devant(b,a)=',res_ba)
