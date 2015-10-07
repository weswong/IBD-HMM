# IBD-HMM
HMM that uses the viterbi algorithm to infer relatedness between genome
Questions/Concerns should be directed to wesleywong@fas.harvard.edu


Requires:
python 2.7.0
pandas 0.16.1
itertools
numpy
scipy

Optional:
vcftools 0.1.12b

Assuming you start off with a VCF file that is properly filtered/contains the information you need:

1) convert VCF file to a tab-delimited file, which is piped into scriptVCF_IBD.py

zcat {vcf.gz} | vcf-to-tab | scriptVCF_IBD.py {ibd_file_name}

2) create a file that contains all the bad samples (>20% uncallable sites)
vcftools --gzvcf {.vcf.gz} --indiv-missing --out {badsamples_file}

python scriptVCF_badsamples.py {badsamples_file}

3) Run the fileprep program for the IBD HMM
python hmm_fileprep sequence_file bad_samples
* Two outputs: a good combos file and a discordance file
For the discordance file, column1 is sample1, column2 is sample 2, column 3 is the number of sites examined for this particular pair, column 4 is the discordance

*Note: hmm_fileprep_coi.py sequence_file bad_samples {optional freq file json}
hmm_fileprep_coi file is the same as hmm_fileprep except it contains an optional parameter to pass
in the frequencies as a json file. These frequencies are created using the freq_parse.py file

4)Run the IBD HMM
python IBD_HMM.py good_combos_file cleaned_file
