# heritability

The main purpose of this project is to study heritability, considering the following datasets:

* 1000 genomes phase 3 - genotype data
* Geuvadis - gene expressions (phenotypes of study) and other features from individuals (sex, population and laboratory whose data were collected)

The structure of this project:

* 00.Descriptive - Programs created to analyse data before studying heritability:
    * 01.Filters - Programs to analyse how filters impact size of genes' set size:
      * 01.1.generateListSnpsPosition.sh - Solely created to produce SNPs' location on genome;
      * 01.2.GenerateHardyMaf.sh - Given different samples - all genotyped individuals, only individuals whose phenotypes are known (labeled as 'Geuvadis') and the subset of europeans from 'Geuvadis' - calculate HW equilibrium and MAF for each SNP;
      * 01.3.hardyComparisons.r - Analyse how filters based on latter metrics impact the size of SNPs' set and generate list of SNPs' based on thresholds to the next program;
      * 02.vif.sh - Based on lists of SNPs' generated on 01.3.hardyComparisons.r, generate new sets SNPs' for different VIF thresholds;
      * 03.comparisonsAllFilters.r  - Analyse how based on thresholds for MAF, HWE and VIF, the resulting SNPs' sets are influenced (in volume);
  *  02.Phenotype - Analysis of phenotypes and its relation with individuals features.
* 01.Pipeline - Pipeline built to study heritability:
  *  01.DataPrepGeneral.sh - Given lists of SNPs' previously generated in 00.Descriptive/01.Filters/, generate environment to produce .bed files and calculate GRM
  *  exec.py - Execute pipeline
* 99.Old - Old codes, used in previous studies


Main document:
https://www.overleaf.com/read/hcqdxcjzxzmr
