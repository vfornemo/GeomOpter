# Assert energy with reference file

import sys
import os
import regex as re

for dir in os.listdir('.'):
    if os.path.isdir(dir):
        mol_name = dir
        ref_file = os.path.join(dir, mol_name+'.out3')
        if not os.path.exists(ref_file):
            print(f'{mol_name} Failed, ref_file not found')
            continue
        bench_file = os.path.join(dir, mol_name+'.log')
        # print('ref_file:', ref_file)
        regex = 'Final energy at mimimum:\s+(-?\d+\.\d+)\s+kcal\/mol'
        with open(ref_file, 'r') as f:
            ref_energy = float(re.search(regex, f.read()).group(1))
        
        regex2 = 'Final energy at mimimum:\s+(-?\d+\.\d+)'
        regex3 = 'Time Usage:\s+(\d+\s+ms)'
        
        with open(bench_file, 'r') as f:
            bench_energy = float(re.search(regex2, f.read()).group(1))

        with open(bench_file, 'r') as f:
            time_usage = re.search(regex3, f.read()).group(1)
        
        if abs(ref_energy - bench_energy) < 1e-6:
            print(f'{mol_name} Passed, time usage: {time_usage}, ref_energy: {ref_energy}, bench_energy: {bench_energy}')
        else:
            print(f'{mol_name} Failed, time usage: {time_usage}, ref_energy: {ref_energy}, bench_energy: {bench_energy}')