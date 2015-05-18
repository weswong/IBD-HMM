#!/usr/bin/env perl
use warnings;
use strict;

my $file = "/Users/wesley/Wirth_Lab_Results/IBD_HMM/Miotto_analysis/Miotto.maf01.mm80.vcf.gz";

my $com = "zcat $file | ./vcf_genotypes.pl";
my $resp = `$com`;
print $resp;
