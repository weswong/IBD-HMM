
from hmm_fileprep_v6 import *

class Psi:
    '''an object that stores the most likely path of the viterbi algorithm. It also records the probabilities of each path
    0 =DBD , 1 = IBD'''
    def __init__(self):
        # last element of each path consists of a tuple: (delta_tminus IBD, delta_tminus_DBD)
        # refers to the path that always assumes the next step is IBD
        self.path1= ''
        # refers to the path that always assumes that the next step is DBD
        self.path2 = ''
        self.snp_list = []

        # end probabilities that tells us which path is the correct one
        self.delta1 = None
        self.delta2 = None
    
    def resolve_path(self):
        if self.delta1 > self.delta2:[
            self.optimal = self.path1
        elif self.delta1 < self.delta2:
            self.optimal = self.path2
        else:
            self.optimal = None
            
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
    return log(b_IBD), log(b_DBD)

    
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
        delta_IBD = b_IBD + log(0.5)
        delta_DBD = b_DBD + log(0.5)
        
        psi.delta1 = delta_IBD
        psi.delta2 = delta_DBD
        
        return psi
    else:
        b_IBD, b_DBD = conditional_probability(comparison, snp2)
        p_tmp_trans = k_rec*rec_rate*(float(snp2.position) - float(snp1.position))
        p_tmp_notrans = 1-p_tmp_trans
        
        p_trans, p_notrans = log(p_tmp_trans), log(p_tmp_notrans)
        psi.snp_list.append(snp1.chr + '.' + snp1.position)
        
        #path 1 first : path 1 assumes that the state of the next point is always IBD
        delta_IBD_tminus, delta_DBD_tminus = psi.delta1, psi.delta2
        log_fout.write(('\t').join([str(snp1.chr) + '_' + snp1.position, str(snp2.chr) + '_' + snp2.position,str(delta_IBD_tminus), str(delta_DBD_tminus), str(p_tmp_trans), str(p_trans), str(p_tmp_notrans), str(p_notrans)]) + '\n')        
        
        delta_IBD = max(delta_IBD_tminus + p_notrans + b_IBD, delta_DBD_tminus + p_trans + b_IBD)        
        IBD_tminus_index = [delta_IBD_tminus + p_notrans + b_IBD, delta_DBD_tminus + p_trans + b_IBD].index(delta_IBD)
        if IBD_tminus_index == 0:
            psi.path1 += '1'
        else:
            psi.path1 += '0'
        
        delta_DBD = max(delta_IBD_tminus + p_trans + b_DBD, delta_DBD_tminus + p_notrans + b_DBD)
        DBD_tminus_index = [delta_IBD_tminus + p_trans + b_DBD, delta_DBD_tminus + p_notrans + b_DBD].index(delta_DBD)
        if DBD_tminus_index == 0:
            psi.path2 += '1'
        else:
            psi.path2 += '0'
        
        # update delta values
        psi.delta1, psi.delta2 = delta_IBD, delta_DBD
        
        return psi
        
def snp_file_parse(snp_file, cleaned_file):
    fin = open(snp_file)
    snp_dict = {}
    for line in fin.readlines()[1:]:
        data =line.strip().split()
        snp_dict[data[0]] = SNP(data[1], data[2], data[3], float(data[4]), data[5], float(data[6]))
    fin.close()
    
    with open(cleaned_file) as cleaned_fin:
        line = np.array(cleaned_fin.readline().strip().split())
        header = line
        for raw_line in cleaned_fin:
            line = np.array(raw_line.strip().split())
            chr = line[0]
            position = line[1]
            snp_name = chr + '_' + position
            for i, element in enumerate(line[3:]):
                snp_dict[snp_name].sample_list_maker(header[i+3], line[i+3])
    return snp_dict

def good_combos_parse(file):
    comparison_list = []
    fin = open(file)
    for line in fin.readlines():
        sample1 = line.strip().split()[0].split(':')[0]
        sample2 = line.strip().split()[0].split(':')[1]
        frac_missing = float(line.strip().split()[1])
        comparison = Comparator(sample1, sample2)
        comparison.frac_miss = frac_missing
        comparison_list.append(comparison)
    return comparison_list

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
            
    for key in IBD_dict:
        IBD_dict[key].resolve_path()
    return IBD_dict


good_combos_file = 'good_combos.txt'
snp_file = 'snp_properties.txt'
cleaned_file = 'Mastercleaned.txt'
comparisons_list = good_combos_parse(good_combos_file)
snp_dict = snp_file_parse(snp_file,cleaned_file)
    
log_fout = open('viterbi_log.txt', 'w')
log_fout.write(('\t').join(['sampl1', 'sampl2','delta_IBD_tminus', 'delta_DBD_tminus', 'p_tmp_trans', 'p_trans', 'p_tmp_notrans', 'p_notrans'])+ '\n')
    
test = IBD_HMM(comparisons_list[2], snp_dict)

log_fout.close()