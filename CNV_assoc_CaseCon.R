#! /usr/bin/Rscript

## CNV_assoc_CaseCon.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: January 2016

## Example usage:
## Rscript CNV_assoc_CaseCon.R [working directory] [job commands] [input] [perm] [output] [save]

## working directory: /path/to/working/directory
## job commands: Job submission commands
## input: PLINK .cfile filename (CNV_data.cnv is entered as 'CNV_data')
## perm: number of permutations
## output: output filename
## save: whether or not to save permuted stats

# wdir <- '/path/to/working/directory'
# job <- "bsub -P PGC_CNV -q hour -R 'rusage[mem=4]' -o CNV_assoc_CaseCon.out"
# input <- 'CNV_data'
# perm <- 1000
# output <- 'CNV_data'
# save <- 'no-save'

## WARNING: saving permuted stats from a large number of permutations will create an extrememly large file (SNP x perm matrix)


## ----- REQUIREMENTS:

## PLINK and R installed
## PLINK .cfile: .cnv / .fam / .map file (need to do --cnv-make-map in PLINK if .map file is not provided)

## read in command line arguments
args <- commandArgs(TRUE)
wdir <- args[1]
job <- args[2]
input <- args[3]
perm <- args[4]
output <- args[5]
save <- args[6]

## change to working directory
setwd(wdir)

## create temporary subdirectory within working directory (will not overwrite existing directory)
system('mkdir -p assoc')

type <- c('','--cnv-del','--cnv-dup')
type_name <- c('cnv','del','dup')


## Get summary statistics and association
for (a in 1:length(type)) {

    if (save == 'no-save') { system(paste(job," plink --noweb --cfile ",input," ",type[a]," --mperm ",perm," --cnv-test-2sided --out assoc/",output,"_",type_name[a],sep="")) }
    if (save == 'save')  { system(paste(job," plink --noweb --cfile ",input," ",type[a]," --mperm ",perm," --cnv-test-2sided --mperm-save-all --out assoc/",output,"_",type_name[a],sep="")) }

} #END of a LOOP


