# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from pandas import DataFrame, read_csv
import operator
import json
import sys
# <codecell>

def parse_file(file):
    fin = open(file)
    data=[x.strip().split() for x in fin.readlines()]
    return data

# <codecell>

def find_major_minor(id, freq):
    data = zip(id, freq)
    data.sort(key = lambda x: x[1], reverse = True)
    major, minor = data[0][0], data[1][0]
    raw_major_freq, raw_minor_freq = float(data[0][1]), float(data[1][1])
    major_freq = raw_major_freq / (raw_major_freq + raw_minor_freq)
    minor_freq = raw_minor_freq / (raw_major_freq + raw_minor_freq)
    return major, major_freq, minor, minor_freq

# <codecell>

class Position:
    def __init__(self, chromosome, position, major, major_freq, minor, minor_freq):
        self.chromosome = chromosome
        self.position = position
        self.major = major
        self.major_freq = major_freq
        self.minor = minor
        self.minor_freq = minor_freq
    
    @classmethod
    def interpret_line(cls, line):
        if line[0] == 'M76611':
            chr = 15
        elif line[0] == 'PFC10_API_IRAB':
            chr = 16
        else:
            chr = int(line[0].split('_')[1])
        
        pos = line[1]
        
        allele_id = []
        allele_freq = []
        for number in range(4, len(line)):
            allele_id.append(line[number].split(':')[0])
            allele_freq.append(line[number].split(':')[1])
            
        
        major, major_freq, minor, minor_freq = find_major_minor(allele_id, allele_freq)
        #print major, major_freq
        #print minor, minor_freq
        return Position(chr, pos,major, major_freq, minor, minor_freq)
                    

# <codecell>

class Simulation:
    def __init__(self, input_file):
        self.file = input_file
        self.file_data = parse_file(input_file)
        
    def run(self):
        snp_dict = {}
        for line in self.file_data[1:]:
            position = Position.interpret_line(line)
            snp_dict[str(position.chromosome) + ':' + str(position.position)] = position.__dict__
        json.dump( snp_dict, open('snp_dict.json', 'w'))
            

# <codecell>

if __name__ == '__main__':
    freq_file = sys.argv[1]
    s = Simulation(freq_file)
    s.run()

# <codecell>


