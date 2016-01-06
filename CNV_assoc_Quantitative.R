#! /usr/bin/Rscript

## CNV_assoc_Quantitative.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: January 2016

## Example usage:
## Rscript CNV_assoc_Quantitative.R [working directory] [job commands] [input] [phe] [perm] [output] [save]

## working directory: /path/to/working/directory
## job commands: Job submission commands
## input: PLINK .cfile filename (CNV_data.cnv is entered as 'CNV_data')
## phe: phenotype file
## perm: number of permutations
## output: output filename
## save: whether or not to save permuted stats

# wdir <- '/path/to/working/directory'
# job <- "bsub -P PGC_CNV -q hour -R 'rusage[mem=4]' -o CNV_assoc_Quantitative.out"
# input <- 'CNV_data'
# phe <- 'CNV_data_residual'
# perm <- 1000
# output <- 'CNV_data'
# save <- 'save'

## WARNING: saving permuted stats from a large number of permutations will create an extrememly large file (SNP x perm matrix)

## ----- REQUIREMENTS:

## PLINK and R installed
## PLINK .cfile: .cnv / .fam / .map file (need to do --cnv-make-map in PLINK if .map file is not provided)
## PLINK .phe file (can have different name than cfile)

## read in command line arguments
args <- commandArgs(TRUE)
wdir <- args[1]
job <- args[2]
input <- args[3]
phe <- args[4]
perm <- args[5]
output <- args[6]
save <- args[7]

## change to working directory
setwd(wdir)

## create temporary subdirectory within working directory (will not overwrite existing directory)
system('mkdir -p assoc')

type <- c('','--cnv-del','--cnv-dup')
type_name <- c('cnv','del','dup')


## Get summary statistics and association
for (a in 1:length(type)) {

  if (save=='no-save') { system(paste(job," plink --noweb --cfile ",input," ",type[a]," --mperm ",perm," --pheno ",phe,".phe --pheno-name AFF_residual --out assoc/",output,"_",type_name[a],sep="")) }
  if (save=='save') { system(paste(job," plink --noweb --cfile ",input," ",type[a]," --mperm ",perm," --pheno ",phe,".phe --pheno-name AFF_residual --mperm-save-all --out assoc/",output,"_",type_name[a],sep="")) }

} #END of a LOOP


