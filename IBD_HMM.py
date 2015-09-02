#@profile
def data_prep(data_df):
    eps = 0.001
    k_rec = 2.0
    rec_rate = 5.8e-7
    #discordance condition
    data_df.insert(7,'d_IBD', log(2*eps * (1-eps)))
    d_DBD = (2*data_df['minor_prop']*data_df['major_prop'])*((1-eps)*(1-eps) + eps*eps) + (data_df['minor_prop']**2 + data_df['major_prop']**2) * 2 * eps * (1-eps)
    data_df.insert(8,'d_DBD', log(d_DBD))
    
    #concordance condition
    data_df.insert(9,'c_IBD', log((1-eps)*(1-eps) + eps*eps))
    c_DBD = data_df['minor_prop']**2*(1-eps)**2 + (data_df['major_prop'])**2*eps**2 + 2*data_df['minor_prop']*(data_df['major_prop'])*eps*(1-eps)
    data_df.insert(10,'c_DBD',log(c_DBD))
    
    diff = data_df.groupby('chrom')['pos'].diff()
    data_df.insert(11, 'diff', diff)
    
    p_trans = k_rec*rec_rate*data_df['diff']
    p_notrans = 1 - p_trans
    data_df.insert(12, 'p_trans', log(p_trans))
    data_df.insert(13, 'p_notrans', log(p_notrans))
    
    
    super_minor_df = DataFrame()
    for sample in data_df.keys()[7:]:
        super_minor_df[sample] = data_df[sample] == data_df['minor']
        
    super_missing_df = DataFrame()
    for sample in data_df.keys()[14:]:
        super_missing_df[sample] = data_df[sample] == '.'
    
    return data_df, super_missing_df, super_minor_df
#logical operator
#True = missing
#False = not missing
#Finds all spots that have some level of missing data, then flips it to serve as an index

# <codecell>
#@profile
def pairwise_df(good_df):
    fout = open(base+'_ibdhmm_results.txt', 'w')
    for sample1, sample2 in zip(good_df[0], good_df[1]):
        print sample1, sample2
        missing_bool= np.logical_or(super_missing_df[sample1],super_missing_df[sample2])
        frac_missing = sum(missing_bool) / float(len(data_df))
        missing_sites = missing_bool[missing_bool == True]
        missing_df = DataFrame(missing_sites)
        missing_df['b_IBD'], missing_df['b_DBD'] =frac_missing, frac_missing
        missing_df['reason'] = 'missing'
            
        minor_df = super_minor_df.drop(missing_df.index)
        
        #Major Major
        major_bool = np.logical_not(np.logical_or(minor_df[sample1], minor_df[sample2]))
        major_sites = (major_bool[major_bool == True])
        major_df = DataFrame(major_sites)
        major_df['b_IBD'], major_df['b_DBD'] = log(0.5), log(0.5)
        major_df['reason'] = 'major'
        
        # discordance
        discord_bool = np.logical_xor(minor_df[sample1], minor_df[sample2])
        discord_sites = discord_bool[discord_bool == True]
        discord_df = DataFrame(discord_sites)
        discord_df['b_IBD'] = data_df['d_IBD'].loc[discord_df.index]
        discord_df['b_DBD'] = data_df['d_DBD'].loc[discord_df.index]
        discord_df['reason'] = 'discordance'
        
        #concordance
        concord_bool = np.logical_and(minor_df[sample1], minor_df[sample2])
        concord_sites = concord_bool[concord_bool == True]
        concord_df = DataFrame(concord_sites)
        concord_df['b_IBD'] = data_df['c_IBD'].loc[concord_df.index]
        concord_df['b_DBD'] = data_df['c_DBD'].loc[concord_df.index]
        concord_df['reason'] = 'concordance'
        
        concatenized_df = (pd.concat([missing_df, major_df, discord_df, concord_df],axis=0, join = 'outer')).sort()
        
        concatenized_df.insert(0,'chrom', data_df['chrom'])
        concatenized_df.insert(1, 'pos', data_df['pos'])
        concatenized_df['p_trans'] = data_df['p_trans']
        concatenized_df['p_notrans'] = data_df['p_notrans']
        
        viterbi_dict = viterbi(concatenized_df)
        
        for chromosome in range(1,17): 
            print len(viterbi_dict[chromosome][1]), len(viterbi_dict[chromosome][0])
	    fout.write(('\t').join([sample1, sample2, str(chromosome)]) +'\t' + (',').join([str(x) for x in viterbi_dict[chromosome][1]]) + '\t' + viterbi_dict[chromosome][0] + '\n')

# <codecell>
#@profile
def viterbi(df):
    viterbi_dict = {}
    for number in range(1,17):
        concatenized_df = df[df['chrom'] == number]
       	psi_0 = []
        psi_1 = []
        start_indices = concatenized_df.loc[concatenized_df['p_trans'].isnull()].index
        for index, b_IBD, b_DBD, p_trans, p_notrans in zip(concatenized_df.index, concatenized_df['b_IBD'],
                                                           concatenized_df['b_DBD'], concatenized_df['p_trans'],
                                                           concatenized_df['p_notrans']):
            
            if index in start_indices:
                delta_IBD = b_IBD + log(0.5)
                delta_DBD = b_DBD + log(0.5)
            else:
                delta_IBD_tminus, delta_DBD_tminus = delta_IBD, delta_DBD
                
                delta_IBD = max(delta_IBD_tminus + p_notrans + b_IBD, delta_DBD_tminus + p_trans + b_IBD)
                #IBD_tminus_index = np.argmax(np.array([delta_IBD_tminus + p_notrans, delta_DBD_tminus + p_trans]))
                IBD_tminus_index = [delta_IBD_tminus + p_notrans, delta_DBD_tminus + p_trans].index(max([delta_IBD_tminus + p_notrans, delta_DBD_tminus + p_trans]))
                psi_0.append(IBD_tminus_index)
                    
                delta_DBD = max(delta_IBD_tminus + p_trans + b_DBD, delta_DBD_tminus + p_notrans + b_DBD)
                #DBD_tminus_index = np.argmax(np.array([delta_IBD_tminus + p_trans, delta_DBD_tminus + p_notrans]))
                DBD_tminus_index = [delta_IBD_tminus + p_trans, delta_DBD_tminus + p_notrans].index(max([delta_IBD_tminus + p_trans, delta_DBD_tminus + p_notrans]))
                psi_1.append(DBD_tminus_index)
        
        #reverse lookup to resolve path
        psi = zip(reversed(psi_0), reversed(psi_1))
        path = ''
        
        if delta_IBD > delta_DBD:
            istate = 0
            for element in psi:
                path+=str(element[istate])
                istate = element[istate]
        elif delta_IBD < delta_DBD:
            istate = 1
            for element in psi:
                path+=str(element[istate])
                istate = element[istate]
        else:
            path = 'None'
            
            
        viterbi_dict[number] = (path, concatenized_df['pos'])
    return viterbi_dict
                                 
                             

# <codecell>

if __name__ == "__main__":
    from pandas import DataFrame, read_csv
    import pandas as pd
    import sys
    from numpy import log
    import numpy as np
    
    good_combos_file = sys.argv[1]
    cleaned_file = sys.argv[2]

    base = good_combos_file.split('_')[0]
    good_df = DataFrame(read_csv(good_combos_file, sep = '\t', header = None))
    data_df = DataFrame(read_csv(cleaned_file, sep = '\t'))
    
    data_df, super_missing_df, super_minor_df = data_prep(data_df)
    
    pairwise_df(good_df)

# <codecell>


