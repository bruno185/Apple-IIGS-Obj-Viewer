import math
# small ad-hoc check for local occlusion on q1.obj
obj='q1.obj'
angle_h=240.0; angle_v=0.0; distance=300.0

# parse
verts=[]; faces=[]
with open(obj,'r',encoding='utf-8') as f:
    for line in f:
        line=line.strip()
        if not line or line.startswith('#'): continue
        if line.startswith('v '):
            p=line.split(); x=float(p[1]); y=float(p[2]); z=float(p[3]); verts.append((x,z,y))
        elif line.startswith('f '):
            parts=line.split(); idx=[]
            for p in parts[1:]: v=p.split('/')[0]; idx.append(int(v)-1)
            faces.append(idx)
# transform

def deg2rad(d): return d*math.pi/180.0

def transform(v):
    x,y,z=v
    ch=math.cos(deg2rad(angle_h)); sh=math.sin(deg2rad(angle_h))
    cv=math.cos(deg2rad(angle_v)); sv=math.sin(deg2rad(angle_v))
    ch_cv=ch*cv; sh_cv=sh*cv; ch_sv=ch*sv; sh_sv=sh*sv
    term1 = x*ch_cv
    term2 = y*sh_cv
    term3 = z*sv
    zo = -term1 - term2 - term3 + distance
    if zo>0:
        xo = -x*sh + y*ch
        yo = - x*ch_sv - y*sh_sv + z*cv
    else:
        xo=yo=0.0
    return xo,yo,zo

v_obs=[None]*len(verts)
for i,v in enumerate(verts): v_obs[i]={'xo':transform(v)[0],'yo':transform(v)[1],'zo':transform(v)[2]}

# compare face 0 and 1
f0=faces[0]; f1=faces[1]

# compute projected bboxes
def proj_pts(face):
    pts=[]
    for vi in face:
        v=v_obs[vi]
        if v['zo']<=0: continue
        pts.append((v['xo']/v['zo'], v['yo']/v['zo'], v['zo']))
    return pts
p0=proj_pts(f0); p1=proj_pts(f1)
print('p0',p0)
print('p1',p1)
if not p0 or not p1:
    print('behind camera')
else:
    x0min=min([p[0] for p in p0]); x0max=max([p[0] for p in p0]); y0min=min([p[1] for p in p0]); y0max=max([p[1] for p in p0])
    x1min=min([p[0] for p in p1]); x1max=max([p[0] for p in p1]); y1min=min([p[1] for p in p1]); y1max=max([p[1] for p in p1])
    oxmin=max(x0min,x1min); oxmax=min(x0max,x1max); oymin=max(y0min,y1min); oymax=min(y0max,y1max)
    print('overlap bbox',oxmin,oxmax,oymin,oymax)
    if oxmax<=oxmin or oymax<=oymin:
        print('no overlap')
    else:
        cx=(oxmin+oxmax)/2; cy=(oymin+oymax)/2
        def point_in_poly(px,py,poly):
            crossings=0
            n=len(poly)
            for i in range(n):
                x0,y0=poly[i][0],poly[i][1]; x1,y1=poly[(i+1)%n][0],poly[(i+1)%n][1]
                cond=((y0<=py) and (y1>py)) or ((y0>py) and (y1<=py))
                if cond:
                    xint = x0 + (py - y0) * (x1 - x0) / (y1 - y0)
                    if xint > px: crossings += 1
            return (crossings & 1)
        in0=point_in_poly(cx,cy,p0); in1=point_in_poly(cx,cy,p1)
        print('centroid in0,in1',in0,in1)
        # grid sample
        dx=oxmax-oxmin; dy=oymax-oymin
        skeep=sswap=sund=svalid=0; sumd0=sumd1=0.0
        for ix in (-1,0,1):
            for iy in (-1,0,1):
                px = cx + ix * (dx * 0.25)
                py = cy + iy * (dy * 0.25)
                if not point_in_poly(px,py,p0): continue
                if not point_in_poly(px,py,p1): continue
                # depth at point in triangle fan
                def depth_at(poly,px,py):
                    n=len(poly)
                    for k in range(1,n-1):
                        Ax,Ay,Az = poly[0]
                        Bx,By,Bz = poly[k]
                        Cx,Cy,Cz = poly[k+1]
                        o1=(Bx-Ax)*(py-Ay) - (By-Ay)*(px-Ax)
                        o2=(Cx-Bx)*(py-By) - (Cy-By)*(px-Bx)
                        o3=(Ax-Cx)*(py-Cy) - (Ay-Cy)*(px-Cx)
                        if not ((o1>=0 and o2>=0 and o3>=0) or (o1<=0 and o2<=0 and o3<=0)): continue
                        denom = (By - Cy)*(Ax - Cx) + (Cx - Bx)*(Ay - Cy)
                        if abs(denom) < 1e-9: return None
                        wA = ((By - Cy)*(px - Cx) + (Cx - Bx)*(py - Cy)) / denom
                        wB = ((Cy - Ay)*(px - Cx) + (Ax - Cx)*(py - Cy)) / denom
                        wC = 1.0 - wA - wB
                        if wA < -1e-5 or wB < -1e-5 or wC < -1e-5: return None
                        return wA*Az + wB*Bz + wC*Cz
                    return None
                d0 = depth_at(p0,px,py); d1 = depth_at(p1,px,py)
                print('sample',px,py,'d0',d0,'d1',d1)
                if d0 is None or d1 is None: continue
                svalid += 1; sumd0 += d0; sumd1 += d1
                epsd=1e-4
                if d0 < d1 - epsd: skeep += 1
                elif d1 < d0 - epsd: sswap += 1
                else: sund += 1
        print('samples',svalid,'keep',skeep,'swap',sswap,'und',sund)
