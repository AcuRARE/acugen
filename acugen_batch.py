import os, sys

directory = os.path.abspath(sys.argv[1])
samples = os.listdir(directory)

a = sorted(samples)
e = enumerate(a)
pair_list = []

for i, f in e:
    if 'R1' in a[i] and 'R2' in a[i+1]:
        if a[i][0:a[i].index('R1')] + a[i][a[i].index('R1')+2]  == a[i+1][0:a[i+1].index('R2')] + a[i+1][a[i+1].index('R2')+2]:
            pair_list.append([os.path.join(directory,a[i]), os.path.join(directory,a[i+1])])

if len(pair_list) == 0:
    raise Exception('No pairs found in directory ' + directory + '.', directory)

import acugen_main

for k in pair_list:
    if sys.argv[2]:
        acugen_main.main(*k, sys.argv[2])
    else:
        acugen_main.main(*k)