import sys
import numpy as np
import itertools
from numpy import unique
from collections import Counter
from pandas import DataFrame, read_csv
import scipy.spatial

def allele_count(df_row):
    counts = Counter(df_row[3:])
    counts.pop('.', 0)
    count = len(counts.keys())
    return count
    
def major_find(df_row):
    line = df_row[3:]
    counts = Counter(line)
    # removes and returns the value of '.'
    missing_count = counts.pop('.', 0)
    major = sorted(counts.items(), key=lambda x: x[1])[1][0]
    return major

def minor_find(df_row):
    line = df_row[3:]
    counts = Counter(line)
    # removes and returns the value of '.'
    missing_count = counts.pop('.', 0)
    minor = sorted(counts.items(), key=lambda x: x[1])[0][0]
    return minor

def minor_prop_find(df_row):
    line = df_row[3:]
    counts = Counter(line)
    # removes and returns the value of '.'
    missing_count = counts.pop('.', 0)
    total = float(sorted(counts.items(), key=lambda x: x[1])[0][1] + sorted(counts.items(), key=lambda x: x[1])[1][1])
    minor_prop = sorted(counts.items(), key=lambda x: x[1])[0][1] / total
    return minor_prop

def major_prop_find(df_row):
    line = df_row[3:]
    counts = Counter(line)
    # removes and returns the value of '.'
    missing_count = counts.pop('.', 0)
    total = float(sorted(counts.items(), key=lambda x: x[1])[0][1] + sorted(counts.items(), key=lambda x: x[1])[1][1])
    major_prop = sorted(counts.items(), key=lambda x: x[1])[1][1] / total
    return major_prop

#@profile
def data_prep(input_file, bad_samples_file):
    min_snpD = 10
    tri_allele= 0
    
    output_file = ('').join(input_file.split('.')[0]) + 'cleaned.txt'
    
    bad_samples = [sample.strip() for sample in open(bad_samples_file)]                                              
    df = DataFrame(read_csv(input_file, sep = '\t'))
    #remove bad samples
    df.drop(bad_samples, inplace = True, axis =1)
    #remove non-biallelic alleles
    df.drop(df[df.apply(allele_count, axis = 1) != 2].index, inplace = True)
    
    #remove SNPs that are too close to one another
    df['diff'] = df.groupby('chrom')['pos'].diff()
    df.fillna('first', inplace = True)
    #df.to_csv('test_df.txt', sep = '\t')
    # BUG NOTE MUST FIX THE DAISY CHAIN PROBLEM
    df = df.query('diff > 10 or diff == "first"')
    df.drop('diff', axis = 1, inplace = True)
    
    #calculate the major and minor allele
    major = df.apply(major_find, axis =1 )
    minor = df.apply(minor_find, axis =1 )
    major_prop = df.apply(major_prop_find, axis =1 )
    minor_prop = df.apply(minor_prop_find, axis = 1)
    
    #inserting this stuff into dataframe for future use
    df.insert(3, 'minor_prop', minor_prop)
    df.insert(3, 'minor', minor)
    df.insert(3, 'major_prop', major_prop)
    df.insert(3, 'major', major)
    
    df.to_csv(output_file, sep = '\t', index= False)
    return df

#@profile
def good_combos_find(df):
    fout = open('good_combos.txt', 'w')
    max_discordance = 0.85
    fout_discordance = open('discordance.txt', 'w')
    #minor_df = DataFrame()
    
    super_missing_df = DataFrame()
    super_minor_df = DataFrame()
    
    for sample in df.keys()[7:]:
        #exclude missing sites
        super_missing_df[sample] = df[sample] == '.'
        super_minor_df[sample] = df[sample] == df['minor']

    def viablility_check(sample1, sample2):
        missing_bool= np.logical_or(super_missing_df[sample1],super_missing_df[sample2])
        missing_sites = missing_bool[missing_bool == True]
        missing_df = DataFrame(missing_sites)
        #sans missing sites
        minor_df = super_minor_df.drop(missing_df.index)
        #jaccard distance = (C_TF + C_FT) / (C_TT + C_TF + C_FT)
        discordance = scipy.spatial.distance.jaccard(minor_df[sample1], minor_df[sample2])
        
        fout_discordance.write(('\t').join([sample1, sample2, str(len(minor_df)), str(discordance)]) + '\n')
        if discordance < max_discordance:
            fout.write(('\t').join([sample1, sample2, str(discordance)]) + '\n')
        
        
    samples = [sample for sample in df.keys()[7:]]
    for combo in itertools.combinations(samples, 2):
        viablility_check(combo[0], combo[1])
    fout.close()
    
if __name__ == "__main__":
    input_file = sys.argv[1]
    bad_samples_file = sys.argv[2]
    df = data_prep(input_file, bad_samples_file)
    good_combos_find(df)
