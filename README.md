# array_CNV_analysis
Author: Daniel P. Howrigan  
Start Date: January 2016

Welcome! This repository is a collection of R scripts used to analyze array-based CNV. The primary themes are getting CNV descriptive statistics, outlier detection, CNV burden, and individual locus association. I've added a small example dataset that mimics a number of proprties and issues related to large-scale CNV analysis. If all works well, the scripts provided should be able to run on the example dataset, and should give the user enough background to pursue further, more refined analyses, on their own data. 

Example Dataset
============

The example dataset compromises 2000 individuals from 4 separate datasets. I have anonymized the IDs to prevent any re-identification with the selected samples. The data is a formatted version of CNV calls that allow for CNV analysis in PLINK, and the phenotypes (.phe) have a selected set of principle components that were calculated from GWAS. For data management purposes outside of this example, it is critically important that the sample identifiers for array CNV data and their resepective principal component scores from array GWAS data have been properly matched. 

For the principal components, I've included 5 PCs that showed association to small CNV (< 100 kb) burden in previous analysis, however the relevant PCs used will likely vary from dataset to dataset. Of note, there may be other relevant covariates for each study, and this example is not necessarily a guide of which covariates are, or are not, important in any given analysis framework.

As these files closely follow the format of PLINK CNV files, it is best to look at the PLINK documentation for full details:
pngu.mgh.harvard.edu/~purcell/plink/cnv.shtml

**CNV_data.cnv** - CNV calls  
FID = Family ID (unique identifier across all datasets)  
IID = Individual ID (unique identifier within each dataset)  
CHR = Chromosome  
BP1 = CNV start position  
BP2 = CNV end position  
TYPE = 1 is deletion or copy number loss, 3 is duplication or copy number gain  
SCORE = Can be any value of interest for CNV - Here it is the number of genes overlapping the CNV provided by the original CNV calls  
SITES = number of genotyped SNPs in the CNV  
  
**CNV_data.fam** - pedigree file  
FID = Family ID  
IID = Individual ID  
PID = Paternal ID (set to zero)  
IID = Maternal ID (set to zero)  
SEX = male=1/female=2  
AFF = Affection status (unaffected=1,affected=2)  
  
**CNV_data.phe** - Phenotype file  
FID  
IID  
AFF = Affection status (control=1,case=2)  
SEX (male=1,female=2)  
ancestry = eur (European)  
CNV_platform = Genotyping chip used to call CNVs  
C1 to C8 = Principal Components 1,2,3,4, and 8 from GWAS  

I've also included a few additional files that are often used. A candidate CNV list, which list the genomic positions of selected regions for exlcusion or inclusion (**hg18_implicated_CNV.txt**), and a gene list (**hg18_refGene_plink.txt**) to map CNVs to specific genes.

