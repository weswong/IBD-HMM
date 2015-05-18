from hmm_fileprep_v5 import *
import cPickle as pickle

 with open('combo_dict.pkl', 'wb') as fp:
    combo_dict = pickle.load(combo_dict, fp)
            
with open('snp_dict.pkl', 'wb') as fp:
    snp_dict = pickle.load(snp_dict, fp)
            
class Sequence:
    snps = 
    
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

def conditional_probability(comparison, snp1):
    # conditional probabilities are calculated per position
    eps = 0.001
    sample1, sample2 = comparison.sample1, comparison.sample2
    major_prop, minor_prop = snp1.major_freq, snp1.minor_freq
    
    # determining conditional probabilities (c_ refers to conditional probability)
    if snp1.sample_list[sample1] == '.' or snp1.sample_list[sample2] == '.':
        # if missing data
        b_IBD, b_DBD = comparison.frac_miss, comparison.frac_miss
    elif snp1.sample_list[sample1] == snp1.sample_list[sample2] and snp1.sample_list[sample1] == snp1.major:
        # if only the major allele is present, we set to 0.5 to avoid bias
        b_IBD, b_DBD = 0.5, 0.5
    elif snp1.sample_list[sample1] == snp1.sample_list[sample2] and snp1.sample_list[sample1] == snp1.minor:        
        # concordance condition and minor allele is is present
        b_IBD = (1-eps)*(1-eps) + eps*eps
        b_DBD = minor_prop**2*(1-eps)**2 + (major_prop)**2*eps**2 + 2*minor_prop*(major_prop)*eps*(1-eps)
    else:
        # discordance condition (which essentially means that minor allele is present)
        b_IBD = 2*eps * (1-eps)
        b_DBD = (2*minor_prop*major_prop)*((1-eps)*(1-eps) + eps*eps) + (minor_prop**2 + major_prop**2) * 2 * eps * (1-eps)
        
        # b_state is shorthand for P(observed | state) I don't know how to make a variable look more like that atm
    return b_IBD, b_DBD

def viterbi(comparison, snp1, snp2=None, psi=None):
    '''accepts SNP and Comparator objects and starts the viterbi algorithm
    if no snp2 is provided, we assume that this is the start of the chain
    Most of the abbrebiations are based off of Rabiner and Juan 1986 'An Introduction to HMM' '''    
    k_rec = 2.0
    rec_rate = 5.8e-7
    sample1, sample2 = comparison.sample1, comparison.sample2
    
    # if the start of the markov chain
    if not snp2:
        psi=Psi()
        b_IBD, b_DBD = conditional_probability(comparison, snp1)
        # delta is the joint probability of P(state) * P(observed | state)
        # pi_IBD = pi_DBD = 0.5 = initial state distribution
        delta_IBD = b_IBD * 0.5
        delta_DBD = b_DBD * 0.5
        
        psi.delta1 = delta_IBD
        psi.delta2 = delta_DBD
        
        return psi
    
    else:
        b_IBD, b_DBD = conditional_probability(comparison, snp2)
        p_trans = k_rec*rec_rate*(float(snp1.position) - float(snp2.position))
        p_notrans = 1-p_trans
        psi.snp_list.append(snp1.chr + '.' + snp1.position)
        
        #path 1 first : path 1 assumes that the state of the next point is always IBD
        delta_IBD_tminus, delta_DBD_tminus = psi.delta1, psi.delta2    
        
        delta_IBD = max(delta_IBD_tminus * p_notrans * b_IBD, delta_DBD_tminus * p_trans * b_IBD)        
        IBD_tminus_index = [delta_IBD_tminus * p_notrans * b_IBD, delta_DBD_tminus * p_trans * b_IBD].index(delta_IBD)
        if IBD_tminus_index == 0:
            psi.path1.append('IBD')
        else:
            psi.path1.append('DBD')
        
        delta_DBD = max(delta_IBD_tminus * p_trans * b_DBD, delta_DBD_tminus * p_notrans * b_DBD)
        DBD_tminus_index = [delta_IBD_tminus * p_trans * b_DBD, delta_DBD_tminus * p_notrans * b_DBD].index(delta_DBD)
        if DBD_tminus_index == 0:
            psi.path2.append('IBD')
        else:
            psi.path2.append('DBD')
        
        # update delta values
        psi.delta1, psi.delta2 = delta_IBD, delta_DBD
        
        return psi
        
        
def IBD_HMM(comparison, snp_dict):
    IBD_dict = {}
    sorting_tmp= [tuple(key.split('_')) for key in snp_dict.keys()]
    sorted_snps = [element[0] + '_' + element[1] for element in sorted(sorting_tmp, key=lambda x: (x[0], int(x[1])))]
    
    for first, second in itertools.izip(sorted_snps, sorted_snps[1:]):
        print first, second
        snp1, snp2 = snp_dict[first], snp_dict[second]
        if snp1.chr not in IBD_dict.keys() and snp2.chr == snp1.chr:
            IBD_dict[snp1.chr] = viterbi(comparison, snp1)
            IBD_dict[snp1.chr] = viterbi(comparison, snp1, snp2, IBD_dict[snp1.chr])
        elif snp1.chr in IBD_dict.keys() and snp2.chr == snp1.chr: 
            IBD_dict[snp1.chr] = viterbi(comparison, snp1, snp2, IBD_dict[snp1.chr])
    return IBD_dict