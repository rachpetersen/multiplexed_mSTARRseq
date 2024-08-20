# multiplexed_mSTARRseq

This project contains code for analyses and figures related to:

Applying multiplexed mSTARR-seq to understand methylation-dependent genetic effects on regulatory element function in diverse 
genomes
Rachel M. Petersen, Christopher M. Vockley, Amanda J. Lea

Specifically, the scripts provided here are:

1) mSTARR_trim_map_counts.sh
Code to: trim raw fastq files and map them to the human genome; split the human genome into non-overlapping 400 bp windows; 
calculate coverage of each 400 bp window with DNA and RNA reads; filter windows for those with ≥ 1 DNA read in at least half
of the replicates in both conditions, and ≥ 1 RNA read in half of the replicates in either condition

2) mSTARR_DE_mashr.Rmd
Code to: normalize counts and run differential expression analysis to test for regulatory activity, ran for data generated from the methylated and unmethylated replicates separately; run mashr to test for differences in regulatory activity between the methylated and unmethylated conditions

3) Comparison_to_prior_studies.Rmd
Code to: reanalyze count data generated from two prior mSTARR-seq studies (Lea et al. 2018 and Johnston et al. 2024) using the same analysis pipeline as the present study; perform a correlation analysis on beta values generated for windows tested in prior studies versus the present study; perform a fisher's exact test to test whether windows identified as regulatory in prior studies were more likely to be identified as regulatory in the current study

4) mSTARR_ASE_mashr.Rmd
Code to: run beta binomial model to test for allele specific expression, ran for data generated from the methylated and unmethylated replicates separately; run mashr to test of differences in allele specific expression between the methylated and unmethylated conditions

5) chromHMM_annotations.Rmd
Code to: calculate the number of regulatory windows and ASE sites that overlap with 15 chromHMM genome annotations; perform fisher's exact test for enrichment of each genomic annotation within regulatory windows and ASE sites

6) SNP_matching.Rmd
Code to: generate background SNP sets that are matched with sites in a test set based on MAF, LD score, and gene density. Based on code and data from Mostafavi et al. 2024 which can be found here: 
    data: https://zenodo.org/records/6618073 
    code: https://github.com/hakha-most/gwas_eqtl

7) GWAS_EWAS_overlap.Rmd
Code to: Intersect MD genetic effect sites with GWAS and EWAS hits

8) Multiplex_subsampling.Rmd
Code to: count the number of SNPs located in regulatory windows when subsampling to different combininations of individuals; count the number of MD genetic effect sites that would be present in the assay when subsampling to different combinations of individuals

Data required as input for each script, as well as the output from the scripts, are hosted on Zenodo (xx). The raw data (FASTQ files) are available on NCBI SRA at accession     PRJNA1137064. 

Note that specific paths in the scripts will not work on your computer and you will need to change them accordingly.

Please contact me at rpetersen42@gmail.com with any questions. 

