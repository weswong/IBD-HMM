#!/usr/bin/env python 

'''simple script that takes the output from vcf-to-tab and converts it into a format the Steve Schaffner's HMM can run.
Should be run as zcat file.vcf.gz | vcf-to-tab | thisscript.py output_file name'''
import sys

output = sys.argv[1]
fout = open(output, 'w')

chrom_dict = {'#CHROM':'chrom',
              'M76611':'15',
              'PFC10_API_IRAB':'16',
              'Pf3D7_01_v3':'1',
              'Pf3D7_02_v3':'2',
              'Pf3D7_03_v3':'3',
              'Pf3D7_04_v3':'4',
              'Pf3D7_05_v3':'5',
              'Pf3D7_06_v3':'6',
              'Pf3D7_07_v3':'7',
              'Pf3D7_08_v3':'8',
              'Pf3D7_09_v3':'9',
              'Pf3D7_10_v3':'10',
              'Pf3D7_11_v3':'11',
              'Pf3D7_12_v3':'12',
              'Pf3D7_13_v3':'13',
              'Pf3D7_14_v3':'14'}              


for line in sys.stdin:
	new = line.split()
        if new[0] =='#CHROM':
            new[0], new[1],new[2] = 'chrom','pos','ref_allele'
        else:
            new[0] = chrom_dict[new[0]]
            base_array  = []
            for element in new[3:]:
                base = element.split('/')[0]
                base_array.append(base)
            new = new[0:3] + base_array
	fout.write(('\t').join(new) + '\n')
fout.close()

