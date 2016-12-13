# IBD-HMM
HMM that uses the viterbi algorithm to infer relatedness between genome
Questions/Concerns should be directed to wesleywong@fas.harvard.edu


* This version of the IBD-HMM was used in the paper titled "Genetic relatedness analysis reveals the cotransmission of genetically related Plasmodium falciparum parasites in Thies, Senegal" but is no longer being actively maintained. For the most recent version, please refer to: https://github.com/glipsnort/hmmIBD


Requires:
python 2.7.0
pandas 0.16.1
itertools
numpy
scipy

Optional:
vcftools 0.1.12b

Input Scripts:

* IBD_HMM.py -- standard IBD HMM, uses population allele frequencies to infer IBD
* IBD_HMM_pedigree_known.py -- non-standard IBD HMM. Used identically to IBD_HMM.py but requires that the VCF file be prefiltered to contain only lab-cross offspring with a known parent. Furthermore, VCF file must be prefiltered to contain positions that are variant amongst the two parental strains.
* calculate_relatedness.py -- calculates the relatedness of parasites using the final HMM output
* freq_parse.py -- script that calculates the allele frequencies from the input file to the HMM.
* graph_ibd_pedigree.py -- script that graphs the IBD segments from the HMM output and maps the IBD segments of that particular polygenomic infection two monogenomic samples
* graph_ibd_results.py -- script that graphs the IBD segments from the HMM output
* hmm_fileprep.py -- script that creates the HMM input files from output of scriptVCF_IBD.py. Input file should only accepts a list of monogenomic infections.
* hmm_fileprep_coi.py -- identical to hmm_fileprep.py but allows user to pass in allele frequencies generated through freq_parse.py
* scriptVCF_IBD.py -- simple script that takes the output from vcf-to-tab and converts it into a format that the HMM can run. Should be run as zcat file.vcf.gz | vcf-to-tab | thisscript.py output_file name
* scriptVCF_badsamples.py -- takes the output of the vcftools --missing-indv function to identify bad samples within the dataset

Standard Pipeline:

Assuming you start off with a VCF file that is properly filtered/contains the information you need:

1) convert VCF file to a tab-delimited file, which is piped into scriptVCF_IBD.py

zcat {vcf.gz} | vcf-to-tab | scriptVCF_IBD.py {ibd_file_name}

2) create a file that contains all the bad samples (>20% uncallable sites)
vcftools --gzvcf {.vcf.gz} --indiv-missing --out {badsamples_file}

python scriptVCF_badsamples.py {badsamples_file}

3) Run the fileprep program for the IBD HMM to prepare the samples
python hmm_fileprep sequence_file bad_samples
This script has 3 outputs: a good combos file, a discordance file, and a the input for the HMM (referred to as a cleaned_file)

Discordance file format:
Tab-delimited file where 
column1 = sample1
column2 = sample 2
column 3 = the number of sites (this version only examines sites where the minor allele is present) examined for this particular pair
column 4 - the discordance

Good combos file:
Tab-delimited file where:
column1 = sample1
column 2 = sample2

HMM input file (cleaned_file):
Cleaned up data by removing sites that are too close from each other and calculates  major and minor allele frequencies.py

To specify allele frequencies, use the hmm_fileprep_coi.py script:
hmm_fileprep_coi file is the same as hmm_fileprep except it contains an optional parameter to pass in the frequencies as a json file. 
hmm_fileprep_coi.py sequence_file bad_samples {optional freq file json}

These frequencies are created using the freq_parse.py script (default titled snp_dict.json) and is a python dictionary where:
* {Chromosome:Position : {'major': identity of the major allele,
                          'major_freq': frequency of the minor allele,
                          'minor': identity of the minor allele,
                          'minor_freq': frequency of the minor allele}
Chromosome and position are provided as integers (example: If chromosome 1, position 2000 = 2000, Chromosome:position should be provided as "1:2000"

4)Run the IBD HMM 
python IBD_HMM.py good_combos_file cleaned_file

4) Run Graphing scripts (graph_ibd_pedigree.py, graph_ibd_results.py)
Input is the output of the hmm file.
