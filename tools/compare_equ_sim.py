import csv
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float
# parameters
path='cube.obj'
H=40; V=20; W=0; D=45.0
verts, faces = parse_obj(path)
xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
face_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)
# load equ.csv
expected = {}
with open('equ.csv','r', newline='') as f:
    r=csv.DictReader(f)
    for row in r:
        fi=int(row['face'])
        expected[fi]={
            'a':float(row['a']), 'b':float(row['b']), 'c':float(row['c']), 'd':float(row['d']),
            'z_min':float(row['z_min']), 'z_mean':float(row['z_mean']), 'z_max':float(row['z_max']),
            'idxs': [int(x) for x in row['vertex_indices'].strip().split()]
        }

# compare
print('face, eq_a, sim_a, eq_b, sim_b, eq_c, sim_c, eq_d, sim_d, zmean_eq, zmean_sim')
for i in range(len(face_info)):
    sim=face_info[i]
    a_sim = to_float(sim['a'])
    b_sim = to_float(sim['b'])
    c_sim = to_float(sim['c'])
    d_sim = to_float(sim['d'])
    zmean_sim = to_float(sim['zmean'])
    eq = expected.get(i)
    print(i, f"{eq['a']:.6f}", f"{a_sim:.6f}", f"{eq['b']:.6f}", f"{b_sim:.6f}", f"{eq['c']:.6f}", f"{c_sim:.6f}", f"{eq['d']:.6f}", f"{d_sim:.6f}", f"{eq['z_mean']:.6f}", f"{zmean_sim:.6f}")