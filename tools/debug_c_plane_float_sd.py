import math
s_angle_h=30; s_angle_v=20; s_distance=39.214455
cos_h = math.cos(math.radians(s_angle_h)); sin_h=math.sin(math.radians(s_angle_h))
cos_v = math.cos(math.radians(s_angle_v)); sin_v=math.sin(math.radians(s_angle_v))
cos_h_cos_v = cos_h * cos_v
sin_h_cos_v = sin_h * cos_v
cos_h_sin_v = cos_h * sin_v
sin_h_sin_v = sin_h * sin_v
distance = s_distance * 4.0
verts=[]
with open('Windows\\3d objects\\q2.obj','r') as f:
    for line in f:
        p=line.strip().split()
        if len(p)>0 and p[0]=='v':
            x=float(p[1]); y=float(p[2]); z=float(p[3])
            verts.append((x,y,z))
# build g_obsv
g_obsv=[]
for (x,y,z) in verts:
    term1 = x * cos_h_cos_v
    term2 = y * sin_h_cos_v
    term3 = z * sin_v
    zo = -term1 - term2 - term3 + distance
    if zo > 0.0:
        xo = -x * sin_h + y * cos_h
        yo = -x * cos_h_sin_v - y * sin_h_sin_v + z * cos_v
    else:
        xo=0.0; yo=0.0
    g_obsv.append((xo,yo,zo))
# face indices from Pascal: [1,2,8,4]
idxs=[1,2,8,4]
Aidx, Bidx, Cidx = idxs[0], idxs[1], idxs[2]
A=g_obsv[Aidx]; B=g_obsv[Bidx]; C=g_obsv[Cidx]
x1,y1,z1=A; x2,y2,z2=B; x3,y3,z3=C
print('obs verts (float):',Aidx,A,Bidx,B,Cidx,C)
a = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
b = -x1*(z2 - z3) + x2*(z1 - z3) - x3*(z1 - z2)
c = x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2)
t1 = y2*z3 - y3*z2; t2 = y1*z3 - y3*z1; t3 = y1*z2 - y2*z1
d = -x1 * t1 + x2 * t2 - x3 * t3
print('raw a,b,c,d =',a,b,c,d)
def arr(r): return 0.0 if abs(r) < 0.01 else r
print('arr a,b,c,d =',arr(a),arr(b),arr(c),arr(d))
