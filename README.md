# array_CNV_analysis
Author: Daniel P. Howrigan

Start Date: January 2016

Welcome! This repository is a collection of R scripts used to analyze array-based CNV. The primary themes are getting CNV descriptive statistics, outlier detection, CNV burden, and individual locus association. I've added a small example dataset that mimics a number of proprties and issues related to large-scale CNV analysis. If all works well, the scripts provided should be able to run on the example dataset, and should give the user enough background to pursue further, more refined analyses, on their own data. 

Example Dataset
============

The example dataset compromises 2000 individuals from 4 separate datasets. I have anonymized the IDs to prevent any re-identification with previous samples. The data is a formatted version of CNV calls that work for CNV analysis in PLINK, and the phenotypes (.phe) have a selected set of principle components that were calculated from GWAS. For data management purposes outside of this example, it is critically important that the sample identifiers for array CNV data and their resepective principal component scores from array GWAS data have been properly matched.  
