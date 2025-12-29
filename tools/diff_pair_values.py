#!/usr/bin/env python3
import csv
from collections import defaultdict

cfile = 'Windows/build_vs/bin/Release/pair_debug_v4.csv'
pfile = 'Windows/build_vs/bin/Release/pair_values_pascal.csv'

def parse_c(fname):
    data = {}
    with open(fname,'r') as f:
        for line in f:
            if line.startswith('#DETAIL'):
                parts = line.strip().split(',')
                # format: #DETAIL,pair,f1,f2,vertex_pq,vid,tempo_obs=...,tempo_world=...
                if len(parts) < 7: continue
                tag = parts[4]
                f1 = int(parts[2]); f2 = int(parts[3])
                if tag.startswith('vertex_'):
                    vid = int(parts[5])
                    # extract tempo_obs and tempo_world
                    obspart = parts[6]
                    worldpart = parts[7] if len(parts)>7 else ''
                    # obspart like 'tempo_obs=-28.856323'
                    to = float(obspart.split('=')[1])
                    tw = float(worldpart.split('=')[1])
                    data[(f1,f2,tag,vid)] = (to, tw)
    return data


def parse_pascal(fname):
    data = {}
    with open(fname,'r') as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            stage = r['stage']
            f1 = int(r['f1']); f2 = int(r['f2'])
            if stage.startswith('vertex_'):
                vid = int(r['idx'])
                tempo_obs = float(r['tempo_obs'])
                tempo_world = float(r['tempo_world'])
                data[(f1,f2,stage,vid)] = (tempo_obs, tempo_world)
    return data

cmap = parse_c(cfile)
pmap = parse_pascal(pfile)

# Find first mismatch
for key in sorted(set(list(cmap.keys())+list(pmap.keys()))):
    cvals = cmap.get(key)
    pvals = pmap.get(key)
    if cvals != pvals:
        print('Mismatch at', key)
        print(' C:', cvals)
        print(' P:', pvals)
        break
else:
    print('No mismatch found')
