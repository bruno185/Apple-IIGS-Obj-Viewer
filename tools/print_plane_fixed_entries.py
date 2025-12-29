PAIR='Windows\\build_vs\\bin\\Release\\pair_values_v4.csv'
with open(PAIR,'r') as f:
    for line in f:
        if 'plane_fixed' in line or 'plane_raw' in line or 'plane_arr' in line:
            print(line.strip())
