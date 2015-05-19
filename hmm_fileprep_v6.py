
import sys
import numpy as np
import itertools
import math
from numpy import unique
from collections import Counter
import cPickle as pickle
import json

def allele_count(allele_list):
    counts = Counter(allele_list)
    counts.pop('.', 0)
    count = len(counts.keys())
    return count

def allele_freq(line):
    missing_count, a1_count, a2_count = 0,0,0
    counts = Counter(line)
    # removes and returns the value of '.'
    missing_count = counts.pop('.', 0)
    #only bi-allelic snps allowed
    if len(counts.keys()) == 2:
        minor = sorted(counts.items(), key=lambda x: x[1])[0][0]
        major = sorted(counts.items(), key=lambda x: x[1])[1][0]
        
        total = float(sorted(counts.items(), key=lambda x: x[1])[0][1] + sorted(counts.items(), key=lambda x: x[1])[1][1])
        
        minor_prop = sorted(counts.items(), key=lambda x: x[1])[0][1] / total
        major_prop = sorted(counts.items(), key=lambda x: x[1])[1][1] / total
        
    return major, major_prop, minor, minor_prop

def assign_allele_freq(cleaned_file):
    allele_freq_list = []
    with open(cleaned_file) as cleaned_fin:
            for raw_line in cleaned_fin:
                line = np.array(raw_line.strip().split())
                if line[0] =='chrom':
                    header = line
                else:
                    major, major_prop, minor, minor_prop = allele_freq(line[3:])
                    allele_freq_list.append([major, major_prop, minor, minor_prop])
    return allele_freq_list
# <codecell>

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def data_prep(input_file, bad_samples_file):
    min_snpD = 10
    tri_allele= 0
    
    output_file = ('').join(input_file.split('.')[0])
    fout = open(output_file + 'cleaned.txt', 'w')
    with open(input_file) as fin:
        pos = {}
        for raw_line in fin:
            raw_line=raw_line.split()
            if raw_line[0] =='chrom':
                bad_samples= []
                with open(bad_samples_file) as fin_bad:
                    for bad_id in fin_bad:
                        bad_id = bad_id.strip()
                        print 'removing {bad_id}'.format(bad_id=bad_id)
                        bad_samples.append(raw_line.index(bad_id))
                header = [i for j, i in enumerate(raw_line) if j not in bad_samples]
                # -3 because first three columns are not samples
                combo_indices = [x for x in itertools.combinations([number for number in range(3,len(header[3:]) + 2)],2)]
                
                print 'nsamples present: {nsamples}'.format(nsamples=str(len(raw_line) -3))
                print 'nsamples used: {nused}'.format(nused=str(len(raw_line) -3))
                print 'nsamples removed: {nremoved}'.format(nremoved=str(len(bad_samples)))
                print 'ncomparisons: {ncomparisons}'.format(ncomparisons=len(combo_indices))
                fout.write(('\t').join(header) + '\n')
                
            else:
                line = [i for j, i in enumerate(raw_line) if j not in bad_samples]
                # start of a new chromosome
                if line[0] not in pos.keys():
                    if allele_count(line[3:]) == 2:
                        pos[line[0]] = line[1]
                        fout.write(('\t').join(line) + '\n')
                    else:
                        tri_allele+=1
                        print 'reject triallele chr:{chr} pos:{pos}'.format(chr=line[0], pos=line[1])
                else:
                    if np.absolute(int(pos[line[0]]) - int(line[1])) < min_snpD:
                        print 'reject chr:{chr} pos:{pos} (too close to previous SNP)'.format(chr=line[0], pos=line[1])
                        continue
                    else:
                        if allele_count(line[3:]) == 2:
                            pos[line[0]] = line[1]
                            fout.write(('\t').join(line) + '\n')
                        else:
                            tri_allele+=1
                            print 'reject triallele chr:{chr} pos:{pos}'.format(chr=line[0], pos=line[1])
    print 'tri_alleles: {tri_allele}'.format(tri_allele=tri_allele)
    fout.close()
    return combo_indices

class Psi:
    '''an object that stores the most likely path of the viterbi algorithm. It also records the probabilities of each path'''
    def __init__(self):
        # last element of each path consists of a tuple: (delta_tminus IBD, delta_tminus_DBD)
        # refers to the path that always assumes the next step is IBD
        self.path1= []
        # refers to the path that always assumes that the next step is DBD
        self.path2 = []
        self.snp_list = []

        # end probabilities that tells us which path is the correct one
        self.delta1 = None
        self.delta2 = None

class Comparator:
    '''an object that stores the number of sites where the minor alele is present in at least one sample AND:
        1) it is the same (concordant)
        2) it is different (discordant)
        OR
        3) number of sites that are missing
        4) Total number of SNPs used overall. This does not correct for number of missing sites'''
    max_discordance = 0.85
    snpcount = 0 #we set the snpcount at runtime
    def __init__(self,sample1, sample2, concordant=0, discordant=0, missing=0):
        self.concordant_sites = float(concordant)
        self.discordant_sites= float(discordant)
        self.missing = float(missing)
        self.sample1, self.sample2 = sample1, sample2
        
    def frac_miss_calc(self):
        try:
            self.frac_miss = float(self.missing) / float(Comparator.snpcount)
        except:
            self.frac_miss = None
    
    def discordance_calc(self):
        try:
            self.discordance = self.discordant_sites / (self.discordant_sites + self.concordant_sites)
        except:
            self.discordance = None
        
    def viable_check(self):
        # determines whether or not we actually use this comparison
        self.discordance_calc()
        if self.discordance > Comparator.max_discordance:
            self.use = False
        else:
            self.use = True
    
    def update(self, c, d, m):
        #updates these values by adding it to the existing one
        self.concordant_sites += c
        self.discordant_sites += d
        self.missing += m



class SNP:
    '''SNP object....kind of self explanatory'''
    def __init__(self, chr, position, major, major_freq, minor, minor_freq):
        self.chr, self.position, self.major, self.major_freq, self.minor, self.minor_freq = chr, position, major, major_freq, minor, minor_freq
        self.sample_list = {}
        self.combo_list = {}
        
    def sample_list_maker(self, sampleID, SNPidentity):
        self.sample_list[sampleID] = SNPidentity

#class Sequence:
#    '''object for each sample.'''
#    def __init__(self, name):
#        self.name = name
#        self.snp_dict = {}
        
def good_combosfinder(cleaned_file):
    '''filters file to get the good combos'''
    #allele_freq_list = assign_allele_freq(cleaned_file)
    Comparator.snpcount = float(file_len(cleaned_file) -1)
    combo_dict = {}
    index = 0
    snp_dict = {}
    fout = open('snp_properties.txt', 'wb')
    fout.write(('\t').join(['name', 'chr', 'position', 'major', 'major_prop', 'minor', 'minor_prop']) + '\n')
    
    with open(cleaned_file) as cleaned_fin:
        line = np.array(cleaned_fin.readline().strip().split())
        header = line
        
        for icombo in combo_indices:
            comparison = str(header[icombo[0]]) + ':' + str(header[icombo[1]])
            # initiate the Comparator object
            combo_dict[comparison] = Comparator(header[icombo[0]],header[icombo[1]])
            
        for raw_line in cleaned_fin:
            line = np.array(raw_line.strip().split())
            chr = line[0]
            position = line[1]
            snp_name = chr + '_' + position
            major, major_prop, minor, minor_prop = allele_freq(line[3:]) 
            fout.write(('\t').join([snp_name, chr, position, major, str(major_prop), minor, str(minor_prop)]) + '\n')

            #snp_dict[snp_name] = SNP(chr, position, major, major_prop, minor, minor_prop)
                    
            #for i, element in enumerate(line[3:]):
            #    snp_dict[snp_name].sample_list_maker(header[i+3], line[i+3])
                                   
            for icombo in combo_indices:
                comparison = str(header[icombo[0]]) + ':' + str(header[icombo[1]])
                same_minor, diff_minor, missing = 0,0,0
                isamp1= icombo[0]
                isamp2= icombo[1]

                # 0 = chrom, 1 = position, 2 = reference allele
                #print major, major_prop, minor, minor_prop
                if line[isamp1] != '.' and line[isamp2] != '.':
                    if line[isamp1] != major or line[isamp2] != major:
                        if line[isamp1] == line[isamp2]:
                            combo_dict[comparison].concordant_sites +=1
                        else:
                            combo_dict[comparison].discordant_sites +=1
                else:
                    combo_dict[comparison].missing += 1
                    
            index += 1
    fout.close()
    good_combos = {}
    for key in combo_dict:
        combo_dict[key].frac_miss_calc()
        combo_dict[key].viable_check()
        if combo_dict[key].use == True:
            good_combos[key] = combo_dict[key]
    
    fout_goodcombos = open('good_combos.txt', 'wb')
    for key in good_combos:
        fout_goodcombos.write(key + '\n')
    
    fout_goodcombos.close()
    
    '''with open('combo_dict.pkl', 'wb') as fp:
        pickle.dump(combo_dict, fp)
            
    with open('snp_dict.pkl', 'wb') as fp:
        pickle.dump(snp_dict, fp)'''
              
                

# <codecell>

if __name__ == "__main__":
    input_file = sys.argv[1]
    bad_samples_file = sys.argv[2]

    cleaned_file = ('').join(input_file.split('.')[0]) + 'cleaned.txt'
    combo_indices = data_prep(input_file, bad_samples_file)
    good_combosfinder(cleaned_file)
    
    
        
        

# <codecell>


