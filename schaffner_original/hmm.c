// Calculate pairwise genoype discordance, conditional on at least one having minor allele 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "translate.h"

int main(int argc, char **argv) {
  // mostly variable declaration
  // states: 0=IBD, 1=DBD
  //  const double mean_trans = 2;  // mean transitions per chromosome
  const char target_pop[32] = "Sen";
  double k_rec;    // recombination coefficient (= effective number of meioses, sort of)
  const double rec_rate = 5.8e-7;  // 5.8e-5 cM/bp, or 17kb/cM
  //eps?
  const double eps = .001;
  const int linesize = 5000;
  // min_disc_call = no discordance, or complete concordance
  const int min_disc_call = 0;
  // if I want to save compute time and ignore things that are very similar
  const double min_discord = -0.01;
  // why .85 for max discord? Based on eye to ignore sequences that are too different
  const double max_discord = 0.85;
  // assumes that we won't have > 150,000 SNPs, reasonable, for falciparum
  int max_snp = 150000;
  const int nchrom = 14;
  // each snp must have 10 snps between
  const int min_snp_sep = 10;
  // what is max_bad?
  const int max_bad = 15;
  // * are pointers that refer to that address of that variable
  char *newLine, *token, *running, **sample, *head, **geno, *majall=NULL, gi, gj;
  // char create an array, file, with 64 elements
  char file[64], *refall=NULL, **bad_samp=NULL;
  // create integer variables for these guys
  int itoken, nsample=0, isamp, chr, sum, iall, maxall, all, js, snp_ind;
  double discord, *Maf=NULL;   
  double hom, *est_sharing=NULL, *fmiss=NULL, *phi[2], pinit[2], b[2], a[2][2], ptrans;
  double maxval, max_phi=0, max_phiL;
  FILE *inf=NULL, *outf=NULL, *chrf=NULL, *pf=NULL;
  int *diff=NULL, *same_min=NULL, jsamp, allcount[5], *use_sample, *traj;
  int nsample_use, nsnp, ipair, npair, isnp, chrlen, *pos, *psi[2], max;
  int *nmiss_bypair=NULL, totall, *start_chr=NULL, *end_chr=NULL, is, maxlen;
  int **use_pair=NULL, nall, killit, nuse_pair=0;
  int printit=0, ntri=0, ibad, nbad;

  //  if (argc != 2) {
  //    fprintf(stderr, "Usage: hmm <recomb coeff>\n");
  //    exit(0);
  //  }
  //  k_rec = strtod(argv[1], NULL);
  k_rec = 2.0;

  bad_samp = malloc(max_bad * sizeof(char*));
  newLine = malloc((linesize+1) * sizeof(char));
  assert(newLine != NULL);
  head = malloc((linesize+1) * sizeof(char));
  assert(head != NULL);
  Maf = calloc(max_snp, sizeof(double));
  majall = malloc(max_snp * sizeof(char));
  refall = malloc(max_snp * sizeof(char));
  pos = malloc(max_snp * sizeof(int));
  start_chr = malloc((nchrom+1) * sizeof(int));
  end_chr = calloc(nchrom+1, sizeof(int));
  // for loop in C format = for (variable initialization; condition; variable update)
  // For (starting when chr =1, if chr <=nchrom; chr +=1
  // to initialize array in C, set some number 10000...
  for (chr = 1; chr <= nchrom; chr++) {start_chr[chr] = 10000000;}
  for (isamp = 0; isamp < max_bad; isamp++) {
    bad_samp[isamp] = malloc(32 * sizeof(char));
  }
  // bad_samples.txt I'm assuming there is a file that lists the files that we belive are terrible
  inf = fopen("data/bad_samples.txt", "r");
  assert(inf != NULL);
  nbad = 0;
  while (fgets(newLine, linesize, inf) != NULL) {
    size_t ln = strlen(newLine) - 1;
    if (newLine[ln] == '\n') {newLine[ln] = '\0';}
    strcpy(bad_samp[nbad], newLine);
    nbad++;
    assert(nbad < max_bad);
  }  
  fclose(inf);
  
  inf = fopen("data/seq.txt", "r");
  assert(inf != NULL);
  sprintf(file, "results/hmm_%s.txt", target_pop);
  outf = fopen(file, "w");
  assert(outf != NULL);
  chrf = fopen("results/chr_ends.txt", "w");
  assert(chrf != NULL);
  fprintf(outf, "pair\tchr\tstart\tend\tdifferent\n");
  fprintf(chrf, "chrom\tstart\tend\n");
  sprintf(file, "results/hmm_probs_%s.txt", target_pop);
  pf = fopen(file, "w");
  assert(pf != NULL);
  fprintf(pf, "pair\tlog_p\n");

  fgets(newLine, linesize, inf); // header
  newLine[strlen(newLine)-1] = 0;  // trailing newline
  if (strlen(newLine) > linesize-3) {fprintf(stderr, "You need a bigger boat\n"); exit(0);}
  strcpy(head, newLine);
  for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
    if (itoken > 1) {nsample++;}
  }

  sample = malloc(nsample * sizeof(char*));
  geno = malloc(nsample * sizeof(char*));
  use_sample = malloc(nsample * sizeof(int));
  use_pair = malloc(nsample * sizeof(int*));
  //isamp = sample index
  for (isamp = 0; isamp < nsample; isamp++) {
    geno[isamp] = malloc(max_snp * sizeof(char));
    assert(geno[isamp] != NULL);
    sample[isamp] = malloc(33 * sizeof(char));
    assert(sample[isamp] != NULL);
    use_pair[isamp] = calloc(nsample, sizeof(int));
    assert(use_pair[isamp] != NULL);
}
  
  isamp = nsample_use = 0;
  for (running = head, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
    if (itoken > 1) {
      strcpy(sample[isamp], token);
      use_sample[isamp] = 1;
      if (strncmp(token, "Ecoli", 5) == 0) {use_sample[isamp] = 0;}
      else if (strncmp(token, "ref_a", 5) == 0) {use_sample[isamp] = 0;}
      else if (strncmp(token, "Preich", 6) == 0) {use_sample[isamp] = 0;}
      else if (strncmp(token, target_pop, strlen(target_pop)) != 0) {use_sample[isamp] = 0;}
      else {
	for (ibad = 0; ibad < nbad; ibad++) {
	  if (strcmp(token, bad_samp[ibad]) == 0) {
	    use_sample[isamp] = 0;
	    fprintf(stderr, "killing %s\n", token);
	    break;
	  }
	}
      }
      if (use_sample[isamp] == 1) {
	nsample_use++;
      }
      isamp++;
    }
  }
  npair = nsample_use * (nsample_use-1) / 2;
  fprintf(stdout, "nsample: %d\t used: %d expected pairs: %d\n", nsample, nsample_use, 
	  npair);

  nmiss_bypair = calloc(npair, sizeof(int));
  same_min = malloc(npair * sizeof(int));
  diff = malloc(npair * sizeof(int));
  //estimated sharing -> not actually used
  est_sharing = malloc(npair * sizeof(double));
  fmiss = malloc(npair * sizeof(double));

  nsnp = 0;
  while (fgets(newLine, linesize, inf) != NULL) {
    totall = killit = 0;
    for (iall = 0; iall < 5; iall++) {allcount[iall] = 0;}
    for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
      if (itoken == 0) {chr = strtol(token, NULL, 10);}
      else if (itoken == 1) {
	pos[nsnp] = strtol(token, NULL, 10);
	if (abs(pos[nsnp] - pos[nsnp-1]) < min_snp_sep) {
	  killit = 1;
	  break;
	}
      }
      else if (itoken == 2) {
	refall[nsnp] = token[0];
      }
      else {
	geno[itoken-2][nsnp] = token[0];
	if (use_sample[itoken-2] == 1) {
	  all = base2index(token[0]);
	  allcount[all]++;
	  if (all > 0) {totall++;}
	}
      }
    }
    // process this snp
    //  find major allele and calculate MAF (major allele freq!)
    if (killit == 1) {continue;}
    maxall = 0;
    nall = 0;
    for (iall = 1; iall <= 4; iall++) {
      if (allcount[iall] > 0) {
	nall++;
      }
      if (allcount[iall] > maxall) {
	maxall = allcount[iall];
	majall[nsnp] = index2base(iall);
      }
    }
    if (maxall == totall) {continue;}
    if (nall > 2) {ntri++;}
    Maf[nsnp] = (double) maxall / totall;
    ipair = 0;
    for (isamp = 0; isamp < nsample; isamp++) {
      if (use_sample[isamp] == 0) {continue;}
      for (jsamp = isamp+1; jsamp < nsample; jsamp++) {
	if (use_sample[jsamp] == 0) {continue;}
	if (geno[isamp][nsnp] != '.' && geno[jsamp][nsnp] != '.') {
	  if (geno[isamp][nsnp] == geno[jsamp][nsnp]) {
	    if (geno[isamp][nsnp] != majall[nsnp]) {same_min[ipair]++;}
	  }
	  else {diff[ipair]++;}
	}
	else {
	  nmiss_bypair[ipair]++;
	}
	ipair++;
      }
    }
    if (nsnp < start_chr[chr]) {
      start_chr[chr] = nsnp;
    }
    if (nsnp > end_chr[chr]) {
      end_chr[chr] = nsnp;
    }
    nsnp++;
    assert(nsnp < max_snp);
  }

  fprintf(stdout, "nsnp: %d\ttriallelic: %d\n", nsnp, ntri);
  ipair = 0;
  for (isamp = 0; isamp < nsample; isamp++) {
    if (use_sample[isamp] == 0) {continue;}
    for (jsamp = isamp+1; jsamp < nsample; jsamp++) {
      if (use_sample[jsamp] == 0) {continue;}
      discord = -1; 
      sum = diff[ipair] + same_min[ipair];
      assert(sum > 0);
      discord = (double) diff[ipair] / sum;
      est_sharing[ipair] = -0.5961 * discord * discord - 0.511 * discord + 0.9898;
      fmiss[ipair] = (double) nmiss_bypair[ipair] / nsnp;
      if (discord <=  max_discord && sum >= min_disc_call && discord >= min_discord) {
	use_pair[isamp][jsamp] = 1;
	nuse_pair++;
      }
      ipair++;
    }
  }  
  fprintf(stderr, "pairs used: %d\n", nuse_pair);

  maxlen = 0;
  for (chr = 1; chr <= nchrom; chr++) {
    fprintf(chrf, "%d\t%d\t%d\n", chr, pos[start_chr[chr]], pos[end_chr[chr]]);
    if (end_chr[chr] - start_chr[chr] + 1 > maxlen) {maxlen = end_chr[chr] - start_chr[chr] + 1;}
  }
  for (is = 0; is < 2; is++) {
    psi[is] = malloc(maxlen * sizeof(int));
    phi[is] = malloc(maxlen * sizeof(double));
  }
  traj = malloc(maxlen * sizeof(int));

  ipair = 0;
  for (isamp = 0; isamp < nsample; isamp++) {
    if (use_sample[isamp] == 0) {continue;}
    for (jsamp = isamp+1; jsamp < nsample; jsamp++) {
      if (use_sample[jsamp] == 0) {continue;}
      if (use_pair[isamp][jsamp] == 1) {
	max_phi = 0;
	printit = 0;
	pinit[0] = 0.5;  // flat prior
	pinit[1] = 0.5;
	for (chr = 1; chr <= nchrom; chr++) {
	  chrlen = pos[end_chr[chr]] - pos[start_chr[chr]];
	  for (isnp = start_chr[chr]; isnp <= end_chr[chr]; isnp++) {
	    //	    if (majall[isnp] != refall[isnp]) {continue;} 
	    snp_ind = isnp - start_chr[chr];
	    gi = geno[isamp][isnp];
	    gj = geno[jsamp][isnp];
	    hom = Maf[isnp]*Maf[isnp] + (1-Maf[isnp]) * (1-Maf[isnp]);
	    if (gi == '.' || gj == '.') {
	      // O = n (missing data)
	      b[0] = b[1] = fmiss[ipair];
	    }
	    else if (gi == gj && gi == majall[isnp]) {
	      //	    else if (gi == gj && gi == refall[isnp]) {
	      // O = aa  major
	      b[0] = (1-eps)*(1-eps) + eps * eps;
	      b[1] = Maf[isnp] * Maf[isnp] * (1-eps)*(1-eps) + 
		(1-Maf[isnp])*(1 - Maf[isnp]) * eps * eps +
		(1 - hom) * eps * (1-eps);
	      b[0] = b[1] = 0.5;  // ignoring homo major allele calls in case of ref bias
	    }
	    else if (gi == gj && gi != majall[isnp]) {
	      //	    else if (gi == gj && gi != refall[isnp]) {
	      // O = AA  minor
	      if (printit == 1) {fprintf(stderr, "    same: chr %d  pos %d\n", chr, pos[isnp]);}
	      b[0] = (1-eps)*(1-eps) + eps * eps;
	      b[1] = (1-Maf[isnp])*(1 - Maf[isnp]) * (1-eps)*(1-eps) + 
		Maf[isnp] * Maf[isnp] * eps * eps +
		(1 - hom) * eps * (1-eps);	
	    }
	    else {
	      // O = aA 
	      if (printit == 1) {fprintf(stderr, "diff: chr %d  pos %d\n", chr, pos[isnp]);}
	      b[0] = 2 * eps * (1 - eps);
	      b[1] = (1 - hom) * ((1-eps)*(1-eps) + eps * eps) + 
		hom * 2 * eps * (1 - eps);
	    }
	    if (isnp == start_chr[chr]) {
	      psi[0][snp_ind] = psi[1][snp_ind] = 0;
	      for (is = 0; is < 2; is++) {phi[is][snp_ind] = log(pinit[is]) + log(b[is]);}
	    }
	    else {
	      //	      ptrans = mean_trans * (pos[isnp] - pos[isnp-1]) / chrlen;
	      ptrans = k_rec * rec_rate * (pos[isnp] - pos[isnp-1]);
	      a[1][0] = ptrans;
	      a[0][1] = ptrans;
	      a[0][0] = 1 - a[0][1];
	      a[1][1] = 1 - a[1][0];
	      for (js = 0; js < 2; js++) {
		maxval = -10000000;
		for (is = 0; is < 2; is++) {
		  if (phi[is][snp_ind-1] + log(a[is][js]) > maxval ) {
		    maxval = phi[is][snp_ind-1] + log(a[is][js]);
		    psi[js][snp_ind] = is;
		  }
		  phi[js][snp_ind] = maxval + log(b[js]);
		}
	      }
	    }   // end if initializing/continuing
	  }  // end snp loop
	  max_phiL = phi[1][snp_ind];
	  if (phi[0][snp_ind] > phi[1][snp_ind]) {max_phiL = phi[0][snp_ind];}
	  max = (phi[1][snp_ind] > phi[0][snp_ind]) ? 1 : 0;
	  traj[snp_ind] = max;
	  max_phi += max_phiL;
	  
	  // inverse loop 
	  for (isnp = snp_ind-1; isnp >= 0; isnp--) {
	    traj[isnp] = psi[max][isnp+1];
	    max = traj[isnp];
	  }
	  // start
	  fprintf(outf, "%s_%s\t%d\t%d", sample[isamp], sample[jsamp], chr, pos[0+start_chr[chr]]);
	  for (isnp = 1; isnp < end_chr[chr] - start_chr[chr] + 1; isnp++) {
	    if (traj[isnp] != traj[isnp-1]) {
	      // end one and start another
	      fprintf(outf, "\t%d\t%d\n", pos[isnp - 1 + start_chr[chr]], traj[isnp-1]);
	      fprintf(outf, "%s_%s\t%d\t%d", sample[isamp], sample[jsamp], chr, pos[isnp + start_chr[chr]]);
	    }
	  }
	  isnp = end_chr[chr] - start_chr[chr];
	  fprintf(outf, "\t%d\t%d\n", pos[end_chr[chr]], traj[isnp]);
	}  // end chrom loop
	fprintf(pf, "%s_%s\t%.5e\n", sample[isamp], sample[jsamp], max_phi);
      }   // end if use pair
      ipair++;
    }
  }
  
  return(0);
}
