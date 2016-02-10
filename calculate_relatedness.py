
import sys
import matplotlib.pyplot as plt
import numpy as np

chr_lengths = {1:643292,
               2:947102,
               3:1060087,
               4:1204112,
               5:1343552,
               6:1418244,
               7:1501717,
               8:1419563,
               9:1541723,
               10: 1687655,
               11:2038337,
               12:2271478,
               13:2895605,
               14:3291871,
               15:5967,
               16:100}

genome_size = 0.
for key in chr_lengths:
    genome_size +=chr_lengths[key]
        
        
class chromosome:
    def __init__(self):
        self.ibd = []
        self.dbd= []

def prepare_data(file):
    fin = open(file)
    data = [x.strip().split() for x in fin.readlines()]
    
    samples_dict = {}
    states_dict = {}
    samp1, samp2 = data[0][0], data[0][1]
    final_line = data[-1]
    
    for line in data:
        if line[0] != samp1 or line[1] != samp2 or line == final_line:
            samples_dict[samp1 + '_' + samp2] = states_dict
            samp1, samp2 = line[0], line[1]
            states_dict = {}

        chr = int(line[2])
        if chr == 15 or chr == 16:
            continue
        elif line[-1] == 'None':
            continue
        else:
            chr_positions = line[-2].split(',')[0:-1]
            states =  list(line[-1])
            if chr not in states_dict.keys():
                states_dict[chr] = [chr_positions, states]
            else:
                print 'error'
                

    for comparison in samples_dict:
        states_dict = samples_dict[comparison]
        for key in states_dict:
            flip_points = []
            chr_positions = states_dict[key][0]
            states = states_dict[key][1]
        
            if states[0] == '0':
                flip_points.append(('IBD', 0))
            else:
                flip_points.append(('DBD', 0))
            
            for i, state in enumerate(states):
                if i != 0:
                    previous_state = states[i-1]
                else:
                    previous_state = states[0]
                
                if state == previous_state:
                    continue
                
                else:
                    if state == '0':
                        flip_points.append(('IBD', int(chr_positions[i])))
                    else:
                        flip_points.append(('DBD', int(chr_positions[i])))
            states_dict[key].append(flip_points)
    
    graph_samples_dict = {}
    graph_dict= {}
    for comparison in samples_dict:
        states_dict= samples_dict[comparison]
        for key in states_dict:
            chr = chromosome()
            for i, element in enumerate(states_dict[key][2]):
                if i != len(states_dict[key][2]) - 1:
                    if element[0] == 'IBD':
                        ibd_range = (element[1], states_dict[key][2][i+1][1] - element[1])
                        chr.ibd.append(ibd_range)
                    else:
                        dbd_range = (element[1], states_dict[key][2][i+1][1] - element[1])
                        chr.dbd.append(dbd_range)
                elif i == len(states_dict[key][2]) -1:
                    if element[0] == 'IBD':
                        ibd_range =(element[1],chr_lengths[key] - element[1])
                        chr.ibd.append(ibd_range)
                    else:
                        dbd_range =(element[1], chr_lengths[key] - element[1])
                        chr.dbd.append(dbd_range)
            graph_dict[key] = chr.__dict__
        graph_samples_dict[comparison] = graph_dict
        graph_dict = {}
    return graph_samples_dict

def calculate_coef_relatedness(file):
    out_file = ('_').join(file.split('_')[0:-2]) +'_relatedness.txt'
    fout = open(out_file, 'w')
    relatedness_dict = prepare_data(file)
    for comparison in relatedness_dict:
        coef_relatedness = 0
        for key in relatedness_dict[comparison]:
            for segment in relatedness_dict[comparison][key]['ibd']:
                coef_relatedness += segment[1] /genome_size
        fout.write(comparison + '\t' + str(coef_relatedness) + '\n')
    fout.close()


if __name__ == '__main__':
    file = sys.argv[1]
    calculate_coef_relatedness(file)