# IBD-HMM
HMM that uses the viterbi algorithm to infer relatedness between genome
Original version written by Steve Schaffner
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
* hmm_fileprep_coi.py -- identical to hmm_fileprep.py but should only contain the two randomly constructed pseudohaplotypes for the polygenomic infection
* scriptVCF_IBD.py simple script that takes the output from vcf-to-tab and converts it into a format that the HMM can run. Should be run as zcat file.vcf.gz | vcf-to-tab | thisscript.py output_file name

Standard Pipeline:

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
