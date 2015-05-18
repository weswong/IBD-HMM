#!/usr/bin/env perl
use warnings;
use strict;

my $keep_barcode = 1;
my $mito = 0;
my $min_call = 0.80;
$| = 1;
my @header = ();
my @samples = ();

my @bc_chr;
my @bc_pos;
if ($keep_barcode == 1) {
    open(IN, "bcode_snp_list_1.txt") || die "could not open whatever . . .\n";
    while (my $line = <IN>) {
	chomp($line);
	my @pieces = split /\t/, $line;
	push @bc_chr, $pieces[0];
	push @bc_pos, $pieces[1];
    }
}

my $file = "data/seq.txt";
if ($mito == 1) {$file = "data/seq_org.txt";}
open(SEQ, ">$file") || die "could not open output\n";
open(DEP, ">data/depth.txt") || die "could not open depth file\n";

my $icode = 0;
while(<STDIN>) {
    if ($_ =~ m/^#CHROM/) {
	@header = split(/\s+/, $_);
	@samples = @header[9..(scalar(@header)-1)];
	print SEQ "chrom\tpos\tref_allele\t", join("\t", @samples);
	print SEQ "\n";
	print DEP "chrom\tpos\tdepth\t", join("\t", @samples);
	print DEP "\n";
	next;
    }
    if ($_ =~ m/^\#/) {next;}

    my @tokens = split /\s+/, $_;
    my $chr = $tokens[0];
    my $match = ($chr =~ m/^Pf3D7_(\d+)_v3/);
    if ($mito == 1) {
	if ($match == 1) {next;}
	if ($chr =~ m/M76611/) {$chr = 15;}
	elsif ($chr =~ m/API/) {$chr = 16;}
    }
    else {
	if ($match != 1) {next;}
	$chr = $1;
	$chr =~ s/^0+//;
    }
    
    my $depth = -1;
    if ($tokens[7] =~ m/DP=(\d+)\;/) {
#	$depth = $1;
    }
    
    $tokens[7] =~ m/AC=(\d+)\;/ or next; 
    my $AC = $1;
    my $ref_all = $tokens[3];
    my @alt_all = split /,/, $tokens[4];
    if ($ref_all !~ m/^[ACGT]$/i) {next;}
    if (length($ref_all) > 1) {next;}
    my $kill = 0;
    foreach my $alt_all (@alt_all) {
	if ($alt_all !~ m/^[ACGT]$/i) {$kill = 1;}
    }
    if ($kill == 1) {next;}
    my $format = $tokens[8];
    my @format = split ":", $format;
    my @genotypes = @tokens[9..(scalar(@header)-1)];
    my $pos = $tokens[1];
    my $outline = "$chr\t$pos\t$ref_all";
    my $depline = $outline;
    my $ncall = 0;
    my $nassay = 0;
    for (my $i = 0; $i < @samples; $i++) {
	my $nfield = 0;
	my @genotype = split ":", $genotypes[$i];
	my $g; 
	for (my $j = 0; $j < @format; $j++) {
	    if ($format[$j] eq "GT") {
		$g = $genotype[$j];
		$nfield++;
#		if ($nfield >= 1) {last;}
	    }
	    elsif ($format[$j] eq "DP") {
		$depth = $genotype[$j];
		if (!defined $depth) {$depth = 0;}
	    }
#	    if ($format[$j] eq "GT") {
#		$g = $genotype[$j];
#		$nfield++;
#		if ($nfield >= 2) {last;}
#	    }
	}
	my $all;
	if ($g =~/0/) {
	    $all = $ref_all;
	    $ncall++;
	}
	elsif ($g =~ m/1/) {
	    $all = $alt_all[0];
	    $ncall++;
	}
	elsif ($g =~ m/2/) {
	    $all = $alt_all[1];
	    $ncall++;
	}
	elsif ($g =~ m/3/) {
	    $all = $alt_all[2];
	    $ncall++;
	}
	else {$all = ".";}
	$nassay++;
	if (length($all) > 1) {print "   $all\n";}
	$outline .= "\t$all";
	$depline .= "\t$depth";
    }
    $outline .= "\n";
    $depline .= "\n";
    if ($ncall / $nassay >= $min_call) {
	print SEQ $outline;
	print DEP $depline;
    }
    else {
	if ($keep_barcode == 1) {
	    while ($icode < 24 && $bc_chr[$icode] < $chr) {
		$icode++;
	    }
	    while ($icode < 24 && $bc_chr[$icode] == $chr && $bc_pos[$icode] < $pos) {
		$icode++;
	    }
	    if ($icode < 24 && $bc_chr[$icode] == $chr && $bc_pos[$icode] == $pos) {
		print SEQ $outline;
	    }
	}
    }
}
