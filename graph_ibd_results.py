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
    
    for line in data:
        if line[0] != samp1 or line[1] != samp2:
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

def calculate_coef_relatedness(file, comparison):
    relatedness_dict = prepare_data(file)
    coef_relatedness = 0
    for key in relatedness_dict[comparison]:
        for segment in relatedness_dict[comparison][key]['ibd']:
            coef_relatedness += segment[1] /genome_size
    return coef_relatedness


def graph_broken_barh(key,graph_dict, comparison):
    plt.broken_barh(graph_dict[key]['ibd'] , (key-0.25, 0.5),label='IBD', facecolors='orange')
    plt.broken_barh(graph_dict[key]['dbd'] , (key-0.25, .5),label='DBD', facecolors = 'blue')
    
    
    plt.legend(loc='right')
    plt.ylim(0,15)
    plt.yticks(np.arange(15))
    plt.xlabel('Position')
    plt.ylabel('Chromosome')
    plt.xscale('linear')
    plt.title(comparison)
    plt.tick_params(axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off',         # ticks along the top edge are off
                    labelbottom='off')


def generate_graph(file):
    plt.figure()
    sample = file.split('_')[0]
    graph_dict = prepare_data(file)
    for comparison in graph_dict:
        for key in graph_dict[comparison]:
            graph_broken_barh(key, graph_dict[comparison], comparison)
        coef_relatedness = calculate_coef_relatedness(file, comparison)
        plt.annotate('prop_IBD={number}'.format(number = coef_relatedness),
                 xy=(0,0),
                 textcoords='figure fraction',
                 xytext=(0.5, 0.2))
        plt.savefig(comparison + '.png')
        plt.close()

if __name__ == '__main__':
    file = sys.argv[1]
    generate_graph(file)
    
    



