from pascal_sim import parse_obj

verts, faces = parse_obj('cube.obj')

# compute plane using first 3 vertex coords (float)
def plane_from_triple(p1,p2,p3):
    x1,y1,z1=p1; x2,y2,z2=p2; x3,y3,z3=p3
    a = y1*(z2-z3) + y2*(z3-z1) + y3*(z1-z2)
    b = -x1*(z2-z3) + x2*(z1-z3) - x3*(z1-z2)
    c = x1*(y2-y3) - x2*(y1-y3) + x3*(y1-y2)
    d = -x1*(y2*z3 - y3*z2) + x2*(y1*z3 - y3*z1) - x3*(y1*z2 - y2*z1)
    return a,b,c,d

for i,idxs in enumerate(faces):
    p1=verts[idxs[0]]; p2=verts[idxs[1]]; p3=verts[idxs[2]]
    a,b,c,d = plane_from_triple(p1,p2,p3)
    print(i, a,b,c,d)
