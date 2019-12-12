# Pipeline for unmatched variant calling, filtering, and tree building used across projects for sequencing normal tissues

This works for both SNVs (CaVEMan) and indels (Pindel).

## Variant calling and initial filtering

SNVs are called by CaVEMan against an in-silico generated matched normal based on the reference genome (GRCh37). 
Indels can be called by Pindel using the same unmatched normal. 
For SNVs, it is recommended to filter on PASS, as well as ASMD>=140 and CLPM==0. For indels, PASS and a quality cutoff of >300.

## Germline filtering

To filter out germline variants, I fit a binomial distribution to the combined read counts of all normal samples from one patient per SNV/indel site, with the total depth as the number of trials, and the total number of reads supporting the variant as number of successes.  
Germline and somatic variants were differentiated based on a one-sided exact binomial test. 
For this test, the null hypothesis is that the number of reads supporting the variants across copy number normal samples is drawn from a binomial with p=0.5 (p=0.95 for copy number equal to one), and the alternative hypothesis drawn from a distribution with p<0.5 (or p<0.95). 
Resulting p-values were corrected for multiple testing with the Benjamini-Hochberg method and a cut-off was set at q < 10-5 to minimize false positives as on average, roughly 40,000 variants were subjected to this statistical test. 
Variants for which the null hypothesis could be rejected were classified as somatic, otherwise as germline. 
When parental genomes were sequenced, de novo germline mutations were taken to be those variants classified as germline in the child, but absent in both parents. 
In addition, it's very advisable to look at the average depth per variant across the normal samples. 
The resulting distribution should roughly look normal; any outliers should be filtered out. 
These might come from repeat regions in the germline and wouldn't be picked up by the exact binomial. 
Be careful to split autosomal and XY-chromosomal variants in male samples.

The code for this exact binomial test can be found in germline_exact_binom.R. 

## Filtering out noise using beta-binomial overdispersion estimation

This works best when individual samples are clonal, i.e. have a VAF distributed around 0.5. 
It will likely work with poly- and oligoclonal samples as well, but the cutoff for rho might need to be adjusted.

To filter out artefact, we fit a beta-binomial distribution to the number of variant supporting reads and total number of reads across samples from the same patient. 
For every somatic SNV, we determined the maximum likelihood overdispersion parameter (ρ) in a grid-based way (ranging the value of ρ from 10-6 to 10-0.05). 
A low overdispersion captures artefactual variants as they appear seemingly randomly across samples and can be modelled as drawn from a binomial distribution. 
In contrast, true somatic variants will be present at a clonal level in some, but not all crypt genomes, and are thus best represented by a beta-binomial with a high overdispersion. 
To distinguish artefacts from true variants, we used ρ=0.1 as a threshold, below which variants were considered artefacts. 
The code for this filtering approach is an adaptation of the Shearwater variant caller (1) and can be found in beta_binom_filter.R.

1. Gerstung M. et al. Subclonal variant calling with multiple samples and prior knowledge. Bioinformatics, 30(9), 1198-1204 (2014).

## Estimate/filter subclonal contribution using a binomial mixture model 




