# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# accepts as input the output of vcftools --missing-indv 
import sys

file = sys.argv[1]
file_out=('.').join(file.split('.')[0:-1])

fout = open(file_out + '.bad_samples.txt', 'w')

with open(file) as fin:
    for line in fin:
        new = line.split()
        if new[0] =='INDV':
            continue
        elif float(new[-1]) > 0.20:
                fout.write(new[0]+'\n')
fout.close()
                
                

