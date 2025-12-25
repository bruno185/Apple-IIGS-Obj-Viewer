import glob
import re
from pathlib import Path

root = Path('.').resolve()
pairdiag = root / 'PAIRDIAG.CSV'
py_files = sorted(root.glob('PAIRPY_*.CSV'))

def extract_py_info(path):
    info = {'file': str(path), 'tri_vote_first3': None, 'tri_vote_newell': None, 'pair_should_swap_obs_first3': None, 'pair_should_swap_obs_newell': None, 'pair_should_swap_model_first3': None, 'pair_should_swap_model_newell': None}
    with path.open('r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line=line.strip()
            if line.startswith('tri_vote,first3'):
                m = re.search(r'decision=(\-?\d+)', line)
                if m: info['tri_vote_first3'] = int(m.group(1))
            if line.startswith('tri_vote,newell'):
                m = re.search(r'decision=(\-?\d+)', line)
                if m: info['tri_vote_newell'] = int(m.group(1))
            if line.startswith('pair_should_swap_obs_first3'):
                parts=line.split(',')
                info['pair_should_swap_obs_first3'] = int(parts[1]) if len(parts)>1 else None
            if line.startswith('pair_should_swap_obs_newell'):
                parts=line.split(',')
                info['pair_should_swap_obs_newell'] = int(parts[1]) if len(parts)>1 else None
            if line.startswith('pair_should_swap_model_first3'):
                parts=line.split(',')
                info['pair_should_swap_model_first3'] = int(parts[1]) if len(parts)>1 else None
            if line.startswith('pair_should_swap_model_newell'):
                parts=line.split(',')
                info['pair_should_swap_model_newell'] = int(parts[1]) if len(parts)>1 else None
    return info


def extract_c_info(path):
    info = {'file': str(path), 'pair_should_swap': None, 'tri_vote': None, 'tri_local': [], 'swap_events': [], 'swap_skipped': [], 'passes': []}
    if not path.exists():
        return info
    with path.open('r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            s=line.strip()
            if s.startswith('pair_should_swap'):
                parts=s.split(',')
                try:
                    info['pair_should_swap'] = int(parts[1])
                except:
                    pass
            if s.startswith('TRI_VOTE'):
                m = re.search(r'TRI_VOTE,?(\d+),(\d+),(\d+)', s)
                if m:
                    info['tri_vote'] = {'keep': int(m.group(1)), 'swap': int(m.group(2)), 'und': int(m.group(3))}
                else:
                    # older CSV style
                    parts = re.split('[,;=]', s)
            if s.startswith('TRI_LOCAL_OCCLUSION'):
                info['tri_local'].append(s)
            if s.startswith('SWAP_EVENT'):
                info['swap_events'].append(s)
            if s.startswith('SWAP_SKIPPED'):
                info['swap_skipped'].append(s)
            if s.startswith('PASS'):
                info['passes'].append(s)
    return info

if __name__ == '__main__':
    py_infos = [extract_py_info(p) for p in py_files]
    c_info = extract_c_info(pairdiag)

    print('\n=== Comparison Report ===\n')
    print('C log: {}'.format(c_info['file']))
    print('  pair_should_swap (final comparator): {}'.format(c_info['pair_should_swap']))
    if c_info['tri_vote']:
        print('  triangle vote: keep={}, swap={}, und={}'.format(c_info['tri_vote']['keep'], c_info['tri_vote']['swap'], c_info['tri_vote']['und']))
    print('  TRI_LOCAL_OCCLUSION entries: {}'.format(len(c_info['tri_local'])) )
    for i,l in enumerate(c_info['tri_local'][:5]):
        print('    [{}] {}'.format(i, l))
    print('  SWAP_EVENT count: {}'.format(len(c_info['swap_events'])))
    print('  SWAP_SKIPPED count: {}'.format(len(c_info['swap_skipped'])))
    print('  PASS entries: {}'.format(len(c_info['passes'])))

    for p in py_infos:
        print('\nPython log: {}'.format(p['file']))
        print('  tri_vote_first3: {}'.format(p['tri_vote_first3']))
        print('  tri_vote_newell: {}'.format(p['tri_vote_newell']))
        print('  pair_should_swap_obs_first3: {}'.format(p['pair_should_swap_obs_first3']))
        print('  pair_should_swap_obs_newell: {}'.format(p['pair_should_swap_obs_newell']))
        print('  pair_should_swap_model_first3: {}'.format(p['pair_should_swap_model_first3']))
        print('  pair_should_swap_model_newell: {}'.format(p['pair_should_swap_model_newell']))

    # Simple comparison note
    print('\nNotes:')
    for p in py_infos:
        # if both tri votes indicate swap but C's pair_should_swap indicates keep or swap
        if (p['tri_vote_first3'] is not None and p['tri_vote_newell'] is not None):
            dec_py = p['tri_vote_first3']
            dec_pyN = p['tri_vote_newell']
            print(' - {}: Python votes (first3={}, newell={})'.format(Path(p['file']).name, dec_py, dec_pyN))
            if c_info['pair_should_swap'] is not None:
                if dec_py == 1 and c_info['pair_should_swap'] == -1:
                    print('   => MISMATCH: Python prefers SWAP but C final comparator returned KEEP')
                elif dec_py == 1 and c_info['pair_should_swap'] == 1:
                    print('   => AGREE: both prefer SWAP')
                elif dec_py == -1 and c_info['pair_should_swap'] == -1:
                    print('   => AGREE: both prefer KEEP')
                else:
                    print('   => Possible minor disagreement (see details above)')

    print('\nReport generated by tools/compare_pair_logs.py')
