

## CNV_burden_platform.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: January 2016


## ----- LSF example usage:

## bsub -P CNV -q priority -R 'rusage[mem=8]' -o CNV_burden_platform.out Rscript --verbose CNV_burden_platform.R [working directory] [input] [gene file] [output]

## working directory: /path/to/working/directory
## input: PLINK .cfile filename (CNV_data.cnv is entered as 'CNV_data')
## gene_file: CHR/START/STOP with gene coordinates
## output: .burden output filename

# wdir <- '/path/to/working/directory'
# input <- 'CNV_data'
# gene_file <- 'hg18_refGene_plink.txt'
# output <- 'CNV_burden_platform'


## ----- REQUIREMENTS:

## PLINK and R installed
## PLINK .cfile: .cnv / .fam / .map file (need to do --cnv-make-map in PLINK if not provided)
## .phe file with platform/dataset information
## .cfile and .phe file must have the same name
## gene list file (currently using hg18_refGene_plink.txt)


## ----- DESCRIPTION:

## Automated script to run CNV burden tests on a variety of datasets, CNV types, and CNV property filters

## 2 main parts:

## 1) Applying dataset and CNV filters

## data_set - Filter on specific datasets
## CNV_type - CNV deletions (losses), duplications (gains), and all CNVs
## region_set - exclude/include specific genomic regions
## CNV_size - Filter on specific CNV size
## CNV_frequency - Filter on specific CNV frequency

## 2) Burden analysis

## PLINK burden analysis, which is limited to case/control comparisons (no covariates included) 
## - Provides case/control descriptive burden statistics
## - uses permutation to generate p-values (not used)
## - Does not split CNVs into genic and non-genic CNVs

## R burden analysis, which includes covariates in logistic regression model
## - uses asymptotic p-values from glm() model (used)
## - Also splits CNVs into genic and non-genic CNVs
## - glm model specification:
## - aff ~ [CNV burden measure] + SEX + CNV_platform + C1 + C2 + C3 + C4 + C8 (full dataset test)
## - aff ~ [CNV burden measure] + SEX + C1 + C2 + C3 + C4 + C8 (within platform test)

## Burden analyses:
## NSEG - number of CNV segments
## KB - total KB of CNV segments
## COUNT - number of genes in CNV segments (as provided by gene file)
## NGENE - number of genes in CNV segments (as provided by SCORE column in .cnv file)
## GENIC - restricting to only CNVs that cover genes
## NONGENIC - restricting to only CNVs that do not cover genes


## ----- SCRIPT SET UP:

## PART 1: CNV burden filters
## PART 2: Variables: Descriptives and Statistics
## PART 3: CNV burden loop
## PART 3.1: Setting up burden loop
## PART 3.2: Applying PLINK filters
## PART 3.3: Checking CNV burden coverage
## PART 3.4: Reading PLINK burden results
## PART 3.5: CNV burden analysis in R
## PART 3.6: Combining all results together and writing to file


## ----- NOTES:
## To change the filters, make adjustments in PART 1 (make sure to adjust both name and filter)
## To change gene list file, make adjustments in the system commands in PART 3.2 
## To change R glm model, make adjustments in the model in PART 3.5
## To change the .status report output, make adjustments in PART 3.1

## KB burden excludes individuals the do not carry any CNVs



## PART 1: CNV burden filters

  ## dataset
  ## CNV type
  ## CNV frequency
  ## CNV size

## read in command line arguments
args <- commandArgs(TRUE)
wdir <- args[1]
input <- args[2]
gene_file <- args[3]
output <- args[4]

## PLINK permutation count (currently ignoring PLINK p-values)
perm <- 100 

## change to working directory
setwd(wdir)

## create temporary subdirectory within working directory (will not overwrite existing directory)
system('mkdir -p burden_loop')

## read in individual data
phe <- read.table(paste(input,'.phe',sep=''),h=T)

## dataset - sort by sizes
dat <- sort(table(phe$CNV_platform),decreasing=T)
data_name <- c('PGC',names(dat))

## write out PLINK ID files 
for (i in 2:length(data_name)) {
    ID <- phe[phe$CNV_platform==data_name[i],1:2]
    write.table(ID,paste('burden_loop/',data_name[i],'.ID',sep=''),col=F,row=F,quo=F,sep='\t')
}

type_name <- c('allCNV','del','dup')
region_name <- c('allregions','novelregions')
size_name <- c('Allsize','100kb','100-200kb','200-500kb','500kb') 
freq_name <- c('Allfreq','singleton','2-10','11plus')

# CNV burden filters
data <- c('',paste(' --keep burden_loop/',data_name[2:length(data_name)],'.ID',sep=''))

type <- c('',
          '--cnv-del',
          '--cnv-dup')

region <- c('','--cnv-exclude hg18_implicated_CNV.txt --cnv-region-overlap 0.01')

size <- c('',
          '--cnv-max-kb 100',
          '--cnv-kb 100 --cnv-max-kb 500',
          '--cnv-kb 500') 

freq <- c('--cnv-overlap 0.50',
          '--cnv-overlap 0 --cnv-freq-exclude-above 1',
          '--cnv-overlap 0.50 --cnv-freq-exclude-above 10',
          '--cnv-overlap 0.50 --cnv-freq-exclude-below 11')



## PART 2: Variables: Descriptives and Statistics

## Note: we can get gene count numbers with --cnv-count

## Files:
##  - .cnv.indiv
##  - .cnv.grp.summary
##  - .cnv.summary.mperm
##  - .cnv with gene annotation

## FILTER sets:

data_set <- NA; region_set <- NA; CNV_type <- NA; CNV_freq <- NA; CNV_size <- NA

## PLINK descriptives and stats:

plink_AFF_CNV <- NA; plink_AFF_RATE <- NA; plink_AFF_PROP <- NA; plink_AFF_TOTKB <- NA; plink_AFF_AVGKB <- NA; plink_AFF_GRATE <- NA; plink_AFF_GPROP <- NA; plink_AFF_GRICH <- NA
plink_UNAFF_CNV <- NA; plink_UNAFF_RATE <- NA; plink_UNAFF_PROP <- NA; plink_UNAFF_TOTKB <- NA; plink_UNAFF_AVGKB <- NA; plink_UNAFF_GRATE <- NA; plink_UNAFF_GPROP <- NA; plink_UNAFF_GRICH <- NA
plink_CASCON_CNV <- NA; plink_CASCON_RATE <- NA; plink_CASCON_PROP <- NA; plink_CASCON_TOTKB <- NA; plink_CASCON_AVGKB <- NA; plink_CASCON_GRATE <- NA; plink_CASCON_GPROP <- NA; plink_CASCON_GRICH <- NA
plink_CASCON_RATE_EMP1 <- NA; plink_CASCON_PROP_EMP1 <- NA; plink_CASCON_TOTKB_EMP1 <- NA; plink_CASCON_AVGKB_EMP1 <- NA; plink_CASCON_GRATE_EMP1 <- NA; plink_CASCON_GPROP_EMP1 <- NA; plink_CASCON_GRICH_EMP1 <- NA

## R descriptives and stats:
  
  ## NSEG sum
  ## case mean
  ## control mean
  ## pearson resid ~ asymptotic t-value
  ## pearson resid ~ asymptotic p-value
  ## permuted p-value for 10,000 permutations of each t-value

FULL_SAMPLE_SIZE <- NA; CNV_CARRIER_SAMPLE_SIZE <- NA; CNV_GENIC_CARRIER_SAMPLE_SIZE <- NA; CNV_NONGENIC_CARRIER_SAMPLE_SIZE <- NA
NSEG <- NA; NSEG_GENIC <- NA; NSEG_NONGENIC <- NA

NSEG_rate <- NA; NSEG_cas_rate <- NA; NSEG_con_rate <- NA; NSEG_cascon_ratio <- NA; NSEG_glm_OR <- NA; NSEG_glm_se <- NA; NSEG_glm_tval <- NA; NSEG_glm_pval <- NA; NSEG_perm_pval <- NA; NSEG_glm_lowerCI <- NA; NSEG_glm_upperCI <- NA 
KB_rate <- NA; KB_cas_rate <- NA; KB_con_rate <- NA; KB_cascon_ratio <- NA; KB_glm_OR <- NA; KB_glm_se <- NA; KB_glm_tval <- NA; KB_glm_pval <- NA; KB_perm_pval <- NA; KB_glm_lowerCI <- NA; KB_glm_upperCI <- NA
COUNT_rate <- NA; COUNT_cas_rate <- NA; COUNT_con_rate <- NA; COUNT_cascon_ratio <- NA; COUNT_glm_OR <- NA; COUNT_glm_se <- NA; COUNT_glm_tval <- NA; COUNT_glm_pval <- NA; COUNT_perm_pval <- NA; COUNT_glm_lowerCI <- NA; COUNT_glm_upperCI <- NA
NGENE_rate <- NA; NGENE_cas_rate <- NA; NGENE_con_rate <- NA; NGENE_cascon_ratio <- NA; NGENE_glm_OR <- NA; NGENE_glm_se <- NA; NGENE_glm_tval <- NA; NGENE_glm_pval <- NA; NGENE_perm_pval <- NA; NGENE_glm_lowerCI <- NA; NGENE_glm_upperCI <- NA

## GENIC CNV
NSEG_GENIC_rate <- NA; NSEG_GENIC_cas_rate <- NA; NSEG_GENIC_con_rate <- NA; NSEG_GENIC_cascon_ratio <- NA; NSEG_GENIC_glm_OR <- NA; NSEG_GENIC_glm_se <- NA; NSEG_GENIC_glm_tval <- NA; NSEG_GENIC_glm_pval <- NA; NSEG_GENIC_perm_pval <- NA; NSEG_GENIC_glm_lowerCI <- NA; NSEG_GENIC_glm_upperCI <- NA 
KB_GENIC_rate <- NA; KB_GENIC_cas_rate <- NA; KB_GENIC_con_rate <- NA; KB_GENIC_cascon_ratio <- NA; KB_GENIC_glm_OR <- NA; KB_GENIC_glm_se <- NA; KB_GENIC_glm_tval <- NA; KB_GENIC_glm_pval <- NA; KB_GENIC_perm_pval <- NA; KB_GENIC_glm_lowerCI <- NA; KB_GENIC_glm_upperCI <- NA

## NONGENIC CNV
NSEG_NONGENIC_rate <- NA; NSEG_NONGENIC_cas_rate <- NA; NSEG_NONGENIC_con_rate <- NA; NSEG_NONGENIC_cascon_ratio <- NA; NSEG_NONGENIC_glm_OR <- NA; NSEG_NONGENIC_glm_se <- NA; NSEG_NONGENIC_glm_tval <- NA; NSEG_NONGENIC_glm_pval <- NA; NSEG_NONGENIC_perm_pval <- NA; NSEG_NONGENIC_glm_lowerCI <- NA; NSEG_NONGENIC_glm_upperCI <- NA 
KB_NONGENIC_rate <- NA; KB_NONGENIC_cas_rate <- NA; KB_NONGENIC_con_rate <- NA; KB_NONGENIC_cascon_ratio <- NA; KB_NONGENIC_glm_OR <- NA; KB_NONGENIC_glm_se <- NA; KB_NONGENIC_glm_tval <- NA; KB_NONGENIC_glm_pval <- NA; KB_NONGENIC_perm_pval <- NA; KB_NONGENIC_glm_lowerCI <- NA; KB_NONGENIC_glm_upperCI <- NA








## PART 3: CNV burden loop


## getting full count of burden analyses
tot <- length(data_name)*length(type_name)*length(region_name)*length(freq_name)*length(size_name)


## PART 3.1: Setting up burden loop

## adjustments to .phe file for burden testinf
phe$matchID <- paste(phe$FID,':',phe$IID,sep='')
phe$aff <- phe$AFF - 1

X <- 0 ## loop COUNTER

for (a in 1:length(data_name)){
for (b in 1:length(type_name)){
for (c in 1:length(region_name)){
for (d in 1:length(freq_name)){
for (e in 1:length(size_name)){
    
X <- X+1

## assign specific loop parameters
data2 <- data[a]
## combining CNV type and region filtering
type2 <- type[b]
region2 <- region[c]
freq2 <- freq[d]
size2 <- size[e]


data_set[X] <- data_name[a]
CNV_type[X] <- type_name[b]
region_set[X] <- region_name[c]
CNV_freq[X] <- freq_name[d]
CNV_size[X] <- size_name[e]

## write out status report
cat('iteration =',X,'of',tot,'\n','data_set',a,'=',data_name[a],'\n','CNV_type',b,'=',type_name[b],'\n','region_set',c,'=',region_name[c],'\n','CNV_freq',d,'=',freq_name[d],'\n','CNV_size',e,'=',size_name[e],'\n',file=paste(output,".status",sep=''),append=F)



## PART 3.2: Applying PLINK filters


## =============== Filter by data, CNV type, and CNV regions first

  system(paste("plink --noweb --cfile ",input," ",data2," ",type2," ",region2," --cnv-write --out burden_loop/type_loop",sep=""))
  ## gene cnv-make-map
  system("plink --noweb --cfile burden_loop/type_loop --cnv-make-map --out burden_loop/type_loop")



## =============== Non 2-10 count CNV burden
if (freq_name[d]!='2-10'){

  ## --- Frequency pruning    
  system(paste("plink --noweb --cfile burden_loop/type_loop ",freq2," --cnv-write --out burden_loop/freq_loop",sep=""))
  ## gene cnv-make-map
  system("plink --noweb --cfile burden_loop/freq_loop --cnv-make-map --out burden_loop/freq_loop")

  
  ## --- Size pruning
  system(paste("plink --noweb --cfile burden_loop/freq_loop ",size2," --cnv-write --out burden_loop/size_loop",sep=""))
  
  ## Adding gene/exon count separately
  system("plink --noweb --cfile burden_loop/size_loop --cnv-make-map --out burden_loop/size_loop")

  system(paste("plink --noweb --cfile burden_loop/size_loop --cnv-indiv-perm --mperm ",perm," --cnv-count ",gene_file," --out burden_loop/burden_loop",sep=""))

} ## END of non-2-6 count levels



## =============== 2-10 count CNV burden
  
if (freq_name[d]=='2-10'){

  ## -- step 1: Remove the > 10 CNV with 50% overlap

    
  ## --- Frequency pruning    
  system(paste("plink --noweb --cfile burden_loop/type_loop ",freq2," --cnv-write --out burden_loop/freq_loop",sep=""))
  ## gene cnv-make-map
  system("plink --noweb --cfile burden_loop/freq_loop --cnv-make-map --out burden_loop/freq_loop")

  
  ## --- Size pruning
  system(paste("plink --noweb --cfile burden_loop/freq_loop ",size2," --cnv-write --out burden_loop/size_loop",sep=""))
  ## gene cnv-make-map
  system("plink --noweb --cfile burden_loop/size_loop --cnv-make-map --out burden_loop/size_loop")
 
  
  ## -- step 2: Now remove the non-overlapped singletons
  system("plink --noweb --cfile burden_loop/size_loop --cnv-overlap 0 --cnv-freq-exclude-below 2 --cnv-write --out burden_loop/size_loop")

  ## Adding gene count separately
  system("plink --noweb --cfile burden_loop/size_loop --cnv-make-map --out burden_loop/size_loop")

  system(paste("plink --noweb --cfile burden_loop/size_loop --cnv-indiv-perm --mperm ",perm," --cnv-count ",gene_file," --out burden_loop/burden_loop",sep=""))

} ## END of 2-10 count levels




## PART 3.3: Checking CNV burden coverage

## Check if any segments were matched
cnv <- read.table('burden_loop/size_loop.cnv',h=T)

if (nrow(cnv)==0){
  
plink_AFF_CNV[X] <- 0; plink_AFF_RATE[X] <- 0; plink_AFF_PROP[X] <- 0; plink_AFF_TOTKB[X] <- 0; plink_AFF_AVGKB[X] <- 0; plink_AFF_GRATE[X] <- 0 ; plink_AFF_GPROP[X] <- 0; plink_AFF_GRICH[X] <- 0
plink_UNAFF_CNV[X] <- 0; plink_UNAFF_RATE[X] <- 0; plink_UNAFF_PROP[X] <- 0; plink_UNAFF_TOTKB[X] <- 0; plink_UNAFF_AVGKB[X] <- 0; plink_UNAFF_GRATE[X] <- 0 ; plink_UNAFF_GPROP[X] <- 0; plink_UNAFF_GRICH[X] <- 0
plink_CASCON_CNV[X] <- NA; plink_CASCON_RATE[X] <- NA; plink_CASCON_PROP[X] <- NA; plink_CASCON_TOTKB[X] <- NA; plink_CASCON_AVGKB[X] <- NA; plink_CASCON_GRATE[X] <- NA ; plink_CASCON_GPROP[X] <- NA; plink_CASCON_GRICH[X] <- NA
plink_CASCON_RATE_EMP1[X] <- NA; plink_CASCON_PROP_EMP1[X] <- NA; plink_CASCON_TOTKB_EMP1[X] <- NA; plink_CASCON_AVGKB_EMP1[X] <- NA; plink_CASCON_GRATE_EMP1[X] <- NA ; plink_CASCON_GPROP_EMP1[X] <- NA; plink_CASCON_GRICH_EMP1[X] <- NA

FULL_SAMPLE_SIZE[X] <- nrow(cnv); CNV_CARRIER_SAMPLE_SIZE[X] <- 0; CNV_GENIC_CARRIER_SAMPLE_SIZE[X] <- 0; CNV_NONGENIC_CARRIER_SAMPLE_SIZE[X] <- 0 
NSEG[X] <- 0; NSEG_GENIC[X] <- 0; NSEG_NONGENIC[X] <- 0

NSEG_rate[X] <- 0; NSEG_cas_rate[X] <- 0; NSEG_con_rate[X] <- 0; NSEG_cascon_ratio[X] <- NA; NSEG_glm_OR[X] <- NA; NSEG_glm_se[X] <- NA; NSEG_glm_tval[X] <- NA; NSEG_glm_pval[X] <- NA; NSEG_glm_lowerCI[X] <- NA; NSEG_glm_upperCI[X] <- NA
KB_rate[X] <- 0; KB_cas_rate[X] <- 0; KB_con_rate[X] <- 0; KB_cascon_ratio[X] <- NA; KB_glm_OR[X] <- NA; KB_glm_se[X] <- NA; KB_glm_tval[X] <- NA; KB_glm_pval[X] <- NA; KB_glm_lowerCI[X] <- NA; KB_glm_upperCI[X] <- NA
COUNT_rate[X] <- 0; COUNT_cas_rate[X] <- 0; COUNT_con_rate[X] <- 0; COUNT_cascon_ratio[X] <- NA; COUNT_glm_OR[X] <- NA; COUNT_glm_se[X] <- NA; COUNT_glm_tval[X] <- NA; COUNT_glm_pval[X] <- NA; COUNT_glm_lowerCI[X] <- NA; COUNT_glm_upperCI[X] <- NA
NGENE_rate[X] <- 0; NGENE_cas_rate[X] <- 0; NGENE_con_rate[X] <- 0; NGENE_cascon_ratio[X] <- NA; NGENE_glm_OR[X] <- NA; NGENE_glm_se[X] <- NA; NGENE_glm_tval[X] <- NA; NGENE_glm_pval[X] <- NA; NGENE_glm_lowerCI[X] <- NA; NGENE_glm_upperCI[X] <- NA

##GENIC
NSEG_GENIC_rate[X] <- 0; NSEG_GENIC_cas_rate[X] <- 0; NSEG_GENIC_con_rate[X] <- 0; NSEG_GENIC_cascon_ratio[X] <- NA; NSEG_GENIC_glm_OR[X] <- NA; NSEG_GENIC_glm_se[X] <- NA; NSEG_GENIC_glm_tval[X] <- NA; NSEG_GENIC_glm_pval[X] <- NA; NSEG_GENIC_glm_lowerCI[X] <- NA; NSEG_GENIC_glm_upperCI[X] <- NA
KB_GENIC_rate[X] <- 0; KB_GENIC_cas_rate[X] <- 0; KB_GENIC_con_rate[X] <- 0; KB_GENIC_cascon_ratio[X] <- NA; KB_GENIC_glm_OR[X] <- NA; KB_GENIC_glm_se[X] <- NA; KB_GENIC_glm_tval[X] <- NA; KB_GENIC_glm_pval[X] <- NA; KB_GENIC_glm_lowerCI[X] <- NA; KB_GENIC_glm_upperCI[X] <- NA

##NONGENIC
NSEG_NONGENIC_rate[X] <- 0; NSEG_NONGENIC_cas_rate[X] <- 0; NSEG_NONGENIC_con_rate[X] <- 0; NSEG_NONGENIC_cascon_ratio[X] <- NA; NSEG_NONGENIC_glm_OR[X] <- NA; NSEG_NONGENIC_glm_se[X] <- NA; NSEG_NONGENIC_glm_tval[X] <- NA; NSEG_NONGENIC_glm_pval[X] <- NA; NSEG_NONGENIC_glm_lowerCI[X] <- NA; NSEG_NONGENIC_glm_upperCI[X] <- NA
KB_NONGENIC_rate[X] <- 0; KB_NONGENIC_cas_rate[X] <- 0; KB_NONGENIC_con_rate[X] <- 0; KB_NONGENIC_cascon_ratio[X] <- NA; KB_NONGENIC_glm_OR[X] <- NA; KB_NONGENIC_glm_se[X] <- NA; KB_NONGENIC_glm_tval[X] <- NA; KB_NONGENIC_glm_pval[X] <- NA; KB_NONGENIC_glm_lowerCI[X] <- NA; KB_NONGENIC_glm_upperCI[X] <- NA

} ## END of NO mapped segments



## PART 3.4: Reading PLINK burden results

if (nrow(cnv) > 0) {
  
## == PLINK results
plink_sum <- read.table("burden_loop/burden_loop.cnv.grp.summary",h=T)
plink_emp <- read.table("burden_loop/burden_loop.cnv.summary.mperm",h=T)

plink_AFF_CNV[X] <- plink_sum[1,3]  
plink_AFF_RATE[X] <- plink_sum[2,3]
plink_AFF_PROP[X] <- plink_sum[3,3]
plink_AFF_TOTKB[X] <- plink_sum[4,3]
plink_AFF_AVGKB[X] <- plink_sum[5,3]
plink_AFF_GRATE[X] <- plink_sum[6,3]
plink_AFF_GPROP[X] <- plink_sum[7,3]
plink_AFF_GRICH[X] <- plink_sum[8,3]
  
plink_UNAFF_CNV[X] <- plink_sum[1,4]  
plink_UNAFF_RATE[X] <- plink_sum[2,4]
plink_UNAFF_PROP[X] <- plink_sum[3,4]
plink_UNAFF_TOTKB[X] <- plink_sum[4,4]
plink_UNAFF_AVGKB[X] <- plink_sum[5,4]
plink_UNAFF_GRATE[X] <- plink_sum[6,4]
plink_UNAFF_GPROP[X] <- plink_sum[7,4]
plink_UNAFF_GRICH[X] <- plink_sum[8,4]

plink_CASCON_CNV[X] <- plink_sum[1,3]/plink_sum[1,4] 
plink_CASCON_RATE[X] <- plink_sum[2,3]/plink_sum[2,4] 
plink_CASCON_PROP[X] <- plink_sum[3,3]/plink_sum[3,4]
plink_CASCON_TOTKB[X] <- plink_sum[4,3]/plink_sum[4,4]
plink_CASCON_AVGKB[X] <- plink_sum[5,3]/plink_sum[5,4]
plink_CASCON_GRATE[X] <- plink_sum[6,3]/plink_sum[6,4]
plink_CASCON_GPROP[X] <- plink_sum[7,3]/plink_sum[7,4]
plink_CASCON_GRICH[X] <- plink_sum[8,3]/plink_sum[8,4]
  
plink_CASCON_RATE_EMP1[X] <- plink_emp[1,3]
plink_CASCON_PROP_EMP1[X] <- plink_emp[2,3]
plink_CASCON_TOTKB_EMP1[X] <- plink_emp[3,3]
plink_CASCON_AVGKB_EMP1[X] <- plink_emp[4,3]
plink_CASCON_GRATE_EMP1[X] <- plink_emp[5,3]
plink_CASCON_GPROP_EMP1[X] <- plink_emp[6,3]
plink_CASCON_GRICH_EMP1[X] <- plink_emp[7,3]


## PART 3.5: CNV burden analysis in R

indiv <- read.table("burden_loop/burden_loop.cnv.indiv",h=T)
indiv$matchID <- paste(indiv$FID,':',indiv$IID,sep='')

cnv <- read.table("burden_loop/size_loop.cnv",h=T)
cnv$matchID <- paste(cnv$FID,':',cnv$IID,sep='')
cnv$bp <- cnv$BP2 - cnv$BP1

## split into genic and non-genic CNVs
genic <- subset(cnv,cnv$SCORE > 0)
nongenic <- subset(cnv,cnv$SCORE == 0)

# -- sum genes for each individual
gene.cnt <- tapply(cnv$SCORE,factor(cnv$matchID),sum)
indiv$NGENE <- 0
# Apply changes only to indexed subset
indx <- match(names(gene.cnt),indiv$matchID)
indiv$NGENE <- replace(indiv$NGENE,indx,gene.cnt)

## -- sum genic CNV count and KB for each individual
genic.tbl <- table(genic$matchID)
indiv$GENIC_CNV_COUNT <- 0
indx <- match(names(genic.tbl),indiv$matchID)
indiv$GENIC_CNV_COUNT <- replace(indiv$GENIC_CNV_COUNT,indx,genic.tbl)

genic.tbl <- tapply(genic$bp,genic$matchID,sum)/1000
indiv$GENIC_KB <- 0
indx <- match(names(genic.tbl),indiv$matchID)
indiv$GENIC_KB <- replace(indiv$GENIC_KB,indx,genic.tbl)


## -- sum non-genic CNV count and KB for each individual
nongenic.tbl <- table(nongenic$matchID)
indiv$NONGENIC_CNV_COUNT <- 0
indx <- match(names(nongenic.tbl),indiv$matchID)
indiv$NONGENIC_CNV_COUNT <- replace(indiv$NONGENIC_CNV_COUNT,indx,nongenic.tbl)

nongenic.tbl <- tapply(nongenic$bp,nongenic$matchID,sum)/1000
indiv$NONGENIC_KB <- 0
indx <- match(names(nongenic.tbl),indiv$matchID)
indiv$NONGENIC_KB <- replace(indiv$NONGENIC_KB,indx,nongenic.tbl)

## -- merge with phenotype
comrg <- merge(indiv,phe,by='matchID',all.x=T)
comrg_cnv <- comrg[comrg$NSEG > 0,]
comrg_genic <- comrg[comrg$GENIC_CNV_COUNT > 0,]
comrg_nongenic <- comrg[comrg$NONGENIC_CNV_COUNT > 0,]



## ==== get SAMPLE_SIZES
FULL_SAMPLE_SIZE[X] <- nrow(comrg) 

CNV_CARRIER_SAMPLE_SIZE[X] <- nrow(comrg_cnv)
CNV_GENIC_CARRIER_SAMPLE_SIZE[X] <- nrow(comrg_genic)
CNV_NONGENIC_CARRIER_SAMPLE_SIZE[X] <- nrow(comrg_nongenic)

NSEG[X] <- sum(comrg$NSEG)
NSEG_GENIC[X] <- sum(comrg$GENIC_CNV_COUNT)
NSEG_NONGENIC[X] <- sum(comrg$NONGENIC_CNV_COUNT)


## ==== CNV Count
NSEG_rate[X] <- mean(comrg$NSEG)
NSEG_cas_rate[X] <- mean(comrg$NSEG[comrg$PHE==2])
NSEG_con_rate[X] <- mean(comrg$NSEG[comrg$PHE==1])
NSEG_cascon_ratio[X] <- NSEG_cas_rate[X]/NSEG_con_rate[X]

if (data_set[X]=='PGC') { NSEG.lm <- glm(aff ~ NSEG + SEX + CNV_platform + C1 + C2 + C3 + C4 + C8,data=comrg,family='binomial') }
if (data_set[X]!='PGC') { NSEG.lm <- glm(aff ~ NSEG + SEX + C1 + C2 + C3 + C4 + C8,data=comrg,family='binomial') }

NSEG.mod <- summary(NSEG.lm)
if(nrow(NSEG.mod$coefficients)==1){
  NSEG_glm_OR[X] <- NA
  NSEG_glm_se[X] <- NA
  NSEG_glm_tval[X] <- NA
  NSEG_glm_pval[X] <- NA
  NSEG_glm_lowerCI[X] <- NA
  NSEG_glm_upperCI[X] <- NA }
if(nrow(NSEG.mod$coefficients) > 1){
NSEG_glm_OR[X] <- exp(NSEG.mod$coefficients[2,1]*NSEG_rate[X])    
NSEG_glm_se[X] <- NSEG.mod$coefficients[2,2]*NSEG_rate[X]  
NSEG_glm_tval[X] <- NSEG.mod$coefficients[2,3]
NSEG_glm_pval[X] <- NSEG.mod$coefficients[2,4]
NSEG_glm_lowerCI[X] <- NSEG_glm_OR[X] - (1.96*NSEG_glm_se[X])
NSEG_glm_upperCI[X] <- NSEG_glm_OR[X] + (1.96*NSEG_glm_se[X]) }

  
## ==== overall KB burden
KB_rate[X] <- mean(comrg_cnv$KB)
KB_cas_rate[X] <- mean(comrg_cnv$KB[comrg_cnv$PHE==2])
KB_con_rate[X] <- mean(comrg_cnv$KB[comrg_cnv$PHE==1])
KB_cascon_ratio[X] <- KB_cas_rate[X]/KB_con_rate[X]

if(sum(comrg_cnv$aff==1)==0 | sum(comrg_cnv$aff==1)==0){
  KB_glm_OR[X] <- NA
  KB_glm_se[X] <- NA
  KB_glm_tval[X] <- NA
  KB_glm_pval[X] <- NA
  KB_glm_lowerCI[X] <- NA
  KB_glm_upperCI[X] <- NA }

if(sum(comrg_cnv$aff==1) > 0 & sum(comrg_cnv$aff==1) > 0){

if (data_set[X]=='PGC') { KB.lm <- glm(aff ~ KB + SEX + CNV_platform + C1 + C2 + C3 + C4 + C8,data=comrg_cnv,family='binomial') }
if (data_set[X]!='PGC') { KB.lm <- glm(aff ~ KB + SEX + C1 + C2 + C3 + C4 + C8,data=comrg_cnv,family='binomial') }

KB.mod <- summary(KB.lm)
KB_glm_OR[X] <- exp(KB.mod$coefficients[2,1]*KB_rate[X])    
KB_glm_se[X] <- KB.mod$coefficients[2,2]*KB_rate[X]  
KB_glm_tval[X] <- KB.mod$coefficients[2,3]
KB_glm_pval[X] <- KB.mod$coefficients[2,4]
KB_glm_lowerCI[X] <- KB_glm_OR[X] - (1.96*KB_glm_se[X])
KB_glm_upperCI[X] <- KB_glm_OR[X] + (1.96*KB_glm_se[X]) }



## ==== avg. COUNT per CNV (sanity check on NGENE)
COUNT_rate[X] <- mean(comrg_cnv$COUNT)
COUNT_cas_rate[X] <- mean(comrg_cnv$COUNT[comrg_cnv$PHE==2])
COUNT_con_rate[X] <- mean(comrg_cnv$COUNT[comrg_cnv$PHE==1])
COUNT_cascon_ratio[X] <- COUNT_cas_rate[X]/COUNT_con_rate[X]

if(sum(comrg_cnv$aff==1)==0 | sum(comrg_cnv$aff==1)==0){
  COUNT_glm_OR[X] <- NA
  COUNT_glm_se[X] <- NA
  COUNT_glm_tval[X] <- NA
  COUNT_glm_pval[X] <- NA
  COUNT_glm_lowerCI[X] <- NA
  COUNT_glm_upperCI[X] <- NA }


if(sum(comrg_cnv$aff==1) > 0 & sum(comrg_cnv$aff==1) > 0){

if (data_set[X]=='PGC') { COUNT.lm <- glm(aff ~ COUNT + SEX + CNV_platform + C1 + C2 + C3 + C4 + C8,data=comrg_cnv,family='binomial') }
if (data_set[X]!='PGC') { COUNT.lm <- glm(aff ~ COUNT + SEX + C1 + C2 + C3 + C4 + C8,data=comrg_cnv,family='binomial') }

COUNT.mod <- summary(COUNT.lm)
COUNT_glm_OR[X] <- exp(COUNT.mod$coefficients[2,1]*COUNT_rate[X])    
COUNT_glm_se[X] <- COUNT.mod$coefficients[2,2]*COUNT_rate[X]  
COUNT_glm_tval[X] <- COUNT.mod$coefficients[2,3]
COUNT_glm_pval[X] <- COUNT.mod$coefficients[2,4]
COUNT_glm_lowerCI[X] <- COUNT_glm_OR[X] - (1.96*COUNT_glm_se[X])
COUNT_glm_upperCI[X] <- COUNT_glm_OR[X] + (1.96*COUNT_glm_se[X]) }


## ==== Genes covered
NGENE_rate[X] <- mean(comrg_cnv$NGENE)
NGENE_cas_rate[X] <- mean(comrg_cnv$NGENE[comrg_cnv$PHE==2])
NGENE_con_rate[X] <- mean(comrg_cnv$NGENE[comrg_cnv$PHE==1])
NGENE_cascon_ratio[X] <- NGENE_cas_rate[X]/NGENE_con_rate[X]

if(sum(comrg_cnv$aff==1)==0 | sum(comrg_cnv$aff==1)==0){
  NGENE_glm_OR[X] <- NA
  NGENE_glm_se[X] <- NA
  NGENE_glm_tval[X] <- NA
  NGENE_glm_pval[X] <- NA
  NGENE_glm_lowerCI[X] <- NA
  NGENE_glm_upperCI[X] <- NA }

if(sum(comrg_cnv$aff==1) > 0 & sum(comrg_cnv$aff==1) > 0){

if (data_set[X]=='PGC') { NGENE.lm <- glm(aff ~ NGENE + SEX + CNV_platform + C1 + C2 + C3 + C4 + C8,data=comrg_cnv,family='binomial') }
if (data_set[X]!='PGC') { NGENE.lm <- glm(aff ~ NGENE + SEX + C1 + C2 + C3 + C4 + C8,data=comrg_cnv,family='binomial') }

NGENE.mod <- summary(NGENE.lm)
NGENE_glm_OR[X] <- exp(NGENE.mod$coefficients[2,1]*NGENE_rate[X])    
NGENE_glm_se[X] <- NGENE.mod$coefficients[2,2]*NGENE_rate[X]  
NGENE_glm_tval[X] <- NGENE.mod$coefficients[2,3]
NGENE_glm_pval[X] <- NGENE.mod$coefficients[2,4]
NGENE_glm_lowerCI[X] <- NGENE_glm_OR[X] - (1.96*NGENE_glm_se[X])
NGENE_glm_upperCI[X] <- NGENE_glm_OR[X] + (1.96*NGENE_glm_se[X]) }


## ==== GENIC CNV Count
NSEG_GENIC_rate[X] <- mean(comrg$GENIC_CNV_COUNT)
NSEG_GENIC_cas_rate[X] <- mean(comrg$GENIC_CNV_COUNT[comrg$PHE==2])
NSEG_GENIC_con_rate[X] <- mean(comrg$GENIC_CNV_COUNT[comrg$PHE==1])
NSEG_GENIC_cascon_ratio[X] <- NSEG_GENIC_cas_rate[X]/NSEG_GENIC_con_rate[X]

if (data_set[X]=='PGC') { NSEG_GENIC.lm <- glm(aff ~ GENIC_CNV_COUNT + SEX + CNV_platform + C1 + C2 + C3 + C4 + C8,data=comrg,family='binomial') }
if (data_set[X]!='PGC') { NSEG_GENIC.lm <- glm(aff ~ GENIC_CNV_COUNT + SEX + C1 + C2 + C3 + C4 + C8,data=comrg,family='binomial') }

NSEG_GENIC.mod <- summary(NSEG_GENIC.lm)
if(nrow(NSEG_GENIC.mod$coefficients)==1){
  NSEG_GENIC_glm_OR[X] <- NA
  NSEG_GENIC_glm_se[X] <- NA
  NSEG_GENIC_glm_tval[X] <- NA
  NSEG_GENIC_glm_pval[X] <- NA
  NSEG_GENIC_glm_lowerCI[X] <- NA
  NSEG_GENIC_glm_upperCI[X] <- NA }
if(nrow(NSEG_GENIC.mod$coefficients) > 1){
NSEG_GENIC_glm_OR[X] <- exp(NSEG_GENIC.mod$coefficients[2,1]*NSEG_GENIC_rate[X])    
NSEG_GENIC_glm_se[X] <- NSEG_GENIC.mod$coefficients[2,2]*NSEG_GENIC_rate[X]  
NSEG_GENIC_glm_tval[X] <- NSEG_GENIC.mod$coefficients[2,3]
NSEG_GENIC_glm_pval[X] <- NSEG_GENIC.mod$coefficients[2,4]
NSEG_GENIC_glm_lowerCI[X] <- NSEG_GENIC_glm_OR[X] - 1.96*NSEG_GENIC_glm_se[X]
NSEG_GENIC_glm_upperCI[X] <- NSEG_GENIC_glm_OR[X] + 1.96*NSEG_GENIC_glm_se[X] }

  
## ==== GENIC KB burden
KB_GENIC_rate[X] <- mean(comrg_genic$GENIC_KB)
KB_GENIC_cas_rate[X] <- mean(comrg_genic$GENIC_KB[comrg_genic$PHE==2])
KB_GENIC_con_rate[X] <- mean(comrg_genic$GENIC_KB[comrg_genic$PHE==1])
KB_GENIC_cascon_ratio[X] <- KB_GENIC_cas_rate[X]/KB_GENIC_con_rate[X]

if(sum(comrg_genic$aff==1)==0 | sum(comrg_genic$aff==1)==0){
  KB_GENIC_glm_OR[X] <- NA
  KB_GENIC_glm_se[X] <- NA
  KB_GENIC_glm_tval[X] <- NA
  KB_GENIC_glm_pval[X] <- NA
  KB_GENIC_glm_lowerCI[X] <- NA
  KB_GENIC_glm_upperCI[X] <- NA }

if(sum(comrg_genic$aff==1) > 0 & sum(comrg_genic$aff==1) > 0){

if (data_set[X]=='PGC') { KB_GENIC.lm <- glm(aff ~ GENIC_KB + SEX + CNV_platform + C1 + C2 + C3 + C4 + C8,data=comrg_genic,family='binomial') }
if (data_set[X]!='PGC') { KB_GENIC.lm <- glm(aff ~ GENIC_KB + SEX + C1 + C2 + C3 + C4 + C8,data=comrg_genic,family='binomial') }

KB_GENIC.mod <- summary(KB_GENIC.lm)
KB_GENIC_glm_OR[X] <- exp(KB_GENIC.mod$coefficients[2,1]*KB_GENIC_rate[X])    
KB_GENIC_glm_se[X] <- KB_GENIC.mod$coefficients[2,2]*KB_GENIC_rate[X]  
KB_GENIC_glm_tval[X] <- KB_GENIC.mod$coefficients[2,3]
KB_GENIC_glm_pval[X] <- KB_GENIC.mod$coefficients[2,4]
KB_GENIC_glm_lowerCI[X] <- exp(KB_GENIC_glm_OR[X] - 1.96*KB_GENIC_glm_se[X])
KB_GENIC_glm_upperCI[X] <- exp(KB_GENIC_glm_OR[X] + 1.96*KB_GENIC_glm_se[X]) }


## ==== NONGENIC CNV Count
NSEG_NONGENIC_rate[X] <- mean(comrg$NONGENIC_CNV_COUNT)
NSEG_NONGENIC_cas_rate[X] <- mean(comrg$NONGENIC_CNV_COUNT[comrg_nongenic$PHE==2])
NSEG_NONGENIC_con_rate[X] <- mean(comrg$NONGENIC_CNV_COUNT[comrg_nongenic$PHE==1])
NSEG_NONGENIC_cascon_ratio[X] <- NSEG_NONGENIC_cas_rate[X]/NSEG_NONGENIC_con_rate[X]

if(sum(comrg$aff==1)==0 | sum(comrg$aff==1)==0){
  NSEG_NONGENIC_glm_OR[X] <- NA
  NSEG_NONGENIC_glm_se[X] <- NA
  NSEG_NONGENIC_glm_tval[X] <- NA
  NSEG_NONGENIC_glm_pval[X] <- NA
  NSEG_NONGENIC_glm_lowerCI[X] <- NA
  NSEG_NONGENIC_glm_upperCI[X] <- NA }

if(sum(comrg$aff==1) > 0 & sum(comrg$aff==1) > 0){

if (data_set[X]=='PGC') { NSEG_NONGENIC.lm <- glm(aff ~ NONGENIC_CNV_COUNT + SEX + CNV_platform + C1 + C2 + C3 + C4 + C8,data=comrg,family='binomial') }
if (data_set[X]!='PGC') { NSEG_NONGENIC.lm <- glm(aff ~ NONGENIC_CNV_COUNT + SEX + C1 + C2 + C3 + C4 + C8,data=comrg,family='binomial') }

NSEG_NONGENIC.mod <- summary(NSEG_NONGENIC.lm)
NSEG_NONGENIC_glm_OR[X] <- exp(NSEG_NONGENIC.mod$coefficients[2,1]*NSEG_NONGENIC_rate[X])    
NSEG_NONGENIC_glm_se[X] <- NSEG_NONGENIC.mod$coefficients[2,2]*NSEG_NONGENIC_rate[X]  
NSEG_NONGENIC_glm_tval[X] <- NSEG_NONGENIC.mod$coefficients[2,3]
NSEG_NONGENIC_glm_pval[X] <- NSEG_NONGENIC.mod$coefficients[2,4]
NSEG_NONGENIC_glm_lowerCI[X] <- NSEG_NONGENIC_glm_OR[X] - 1.96*NSEG_NONGENIC_glm_se[X]
NSEG_NONGENIC_glm_upperCI[X] <- NSEG_NONGENIC_glm_OR[X] + 1.96*NSEG_NONGENIC_glm_se[X] }

  
## ==== NONGENIC KB burden
KB_NONGENIC_rate[X] <- mean(comrg_nongenic$NONGENIC_KB)
KB_NONGENIC_cas_rate[X] <- mean(comrg_nongenic$NONGENIC_KB[comrg_nongenic$PHE==2])
KB_NONGENIC_con_rate[X] <- mean(comrg_nongenic$NONGENIC_KB[comrg_nongenic$PHE==1])
KB_NONGENIC_cascon_ratio[X] <- KB_NONGENIC_cas_rate[X]/KB_NONGENIC_con_rate[X]

if(sum(comrg_nongenic$aff==1)==0 | sum(comrg_genic$aff==1)==0){
  KB_NONGENIC_glm_OR[X] <- NA
  KB_NONGENIC_glm_se[X] <- NA
  KB_NONGENIC_glm_tval[X] <- NA
  KB_NONGENIC_glm_pval[X] <- NA
  KB_NONGENIC_glm_lowerCI[X] <- NA
  KB_NONGENIC_glm_upperCI[X] <- NA }


if(sum(comrg_nongenic$aff==1) > 0 & sum(comrg_nongenic$aff==1) > 0){

if (data_set[X]=='PGC') { KB_NONGENIC.lm <- glm(aff ~ NONGENIC_KB + SEX + CNV_platform + C1 + C2 + C3 + C4 + C8,data=comrg_nongenic,family='binomial') }
if (data_set[X]!='PGC') { KB_NONGENIC.lm <- glm(aff ~ NONGENIC_KB + SEX + C1 + C2 + C3 + C4 + C8,data=comrg_nongenic,family='binomial') }

KB_NONGENIC.mod <- summary(KB_NONGENIC.lm)
KB_NONGENIC_glm_OR[X] <- exp(KB_NONGENIC.mod$coefficients[2,1]*KB_NONGENIC_rate[X])    
KB_NONGENIC_glm_se[X] <- KB_NONGENIC.mod$coefficients[2,2]*KB_NONGENIC_rate[X]  
KB_NONGENIC_glm_tval[X] <- KB_NONGENIC.mod$coefficients[2,3]
KB_NONGENIC_glm_pval[X] <- KB_NONGENIC.mod$coefficients[2,4]
KB_NONGENIC_glm_lowerCI[X] <- exp(KB_NONGENIC_glm_OR[X] - 1.96*KB_NONGENIC_glm_se[X])
KB_NONGENIC_glm_upperCI[X] <- exp(KB_NONGENIC_glm_OR[X] + 1.96*KB_NONGENIC_glm_se[X]) }


} ## End of mapped segments analysis 


    
} ## e LOOP
} ## d LOOP
} ## c LOOP
} ## b LOOP
} ## a LOOP
## ================= END of LOOPs


## PART 3.6: Combining all results together and writing to file

full_data <- cbind.data.frame(data_set,
                              region_set,
                              CNV_type,
                              CNV_freq,
                              CNV_size,
                              plink_AFF_CNV,
                              plink_AFF_RATE,
                              plink_AFF_PROP,
                              plink_AFF_TOTKB,
                              plink_AFF_AVGKB,
                              plink_AFF_GRATE,
                              plink_AFF_GPROP,
                              plink_AFF_GRICH,
                              plink_UNAFF_CNV,
                              plink_UNAFF_RATE,
                              plink_UNAFF_PROP,
                              plink_UNAFF_TOTKB,
                              plink_UNAFF_AVGKB,
                              plink_UNAFF_GRATE,
                              plink_UNAFF_GPROP,
                              plink_UNAFF_GRICH,
                              plink_CASCON_CNV,
                              plink_CASCON_RATE,
                              plink_CASCON_PROP,
                              plink_CASCON_TOTKB,
                              plink_CASCON_AVGKB,
                              plink_CASCON_GRATE,
                              plink_CASCON_GPROP,
                              plink_CASCON_GRICH,
                              plink_CASCON_RATE_EMP1,
                              plink_CASCON_PROP_EMP1,
                              plink_CASCON_TOTKB_EMP1,
                              plink_CASCON_AVGKB_EMP1,
                              plink_CASCON_GRATE_EMP1,
                              plink_CASCON_GPROP_EMP1,
                              plink_CASCON_GRICH_EMP1,
                              FULL_SAMPLE_SIZE,
                              CNV_CARRIER_SAMPLE_SIZE,
                              CNV_GENIC_CARRIER_SAMPLE_SIZE,
                              CNV_NONGENIC_CARRIER_SAMPLE_SIZE,
                              NSEG,
                              NSEG_GENIC,
                              NSEG_NONGENIC,
                              NSEG_rate,
                              NSEG_cas_rate,
                              NSEG_con_rate,
                              NSEG_cascon_ratio,
                              NSEG_glm_OR,
                              NSEG_glm_se,
                              NSEG_glm_tval,
                              NSEG_glm_pval,
                              NSEG_glm_lowerCI,
                              NSEG_glm_upperCI,
                              KB_rate,
                              KB_cas_rate,
                              KB_con_rate,
                              KB_cascon_ratio,
                              KB_glm_OR,
                              KB_glm_se,
                              KB_glm_tval,
                              KB_glm_pval,
                              KB_glm_lowerCI,
                              KB_glm_upperCI,
                              COUNT_rate,
                              COUNT_cas_rate,
                              COUNT_con_rate,
                              COUNT_cascon_ratio,
                              COUNT_glm_OR,
                              COUNT_glm_se,
                              COUNT_glm_tval,
                              COUNT_glm_pval,
                              COUNT_glm_lowerCI,
                              COUNT_glm_upperCI,
                              NGENE_rate,
                              NGENE_cas_rate,
                              NGENE_con_rate,
                              NGENE_cascon_ratio,
                              NGENE_glm_OR,
                              NGENE_glm_se,
                              NGENE_glm_tval,
                              NGENE_glm_pval,
                              NGENE_glm_lowerCI,
                              NGENE_glm_upperCI,
                              NSEG_GENIC_rate,
                              NSEG_GENIC_cas_rate,
                              NSEG_GENIC_con_rate,
                              NSEG_GENIC_cascon_ratio,
                              NSEG_GENIC_glm_OR,
                              NSEG_GENIC_glm_se,
                              NSEG_GENIC_glm_tval,
                              NSEG_GENIC_glm_pval,
                              NSEG_GENIC_glm_lowerCI,
                              NSEG_GENIC_glm_upperCI,
                              KB_GENIC_rate,
                              KB_GENIC_cas_rate,
                              KB_GENIC_con_rate,
                              KB_GENIC_cascon_ratio,
                              KB_GENIC_glm_OR,
                              KB_GENIC_glm_se,
                              KB_GENIC_glm_tval,
                              KB_GENIC_glm_pval,
                              KB_GENIC_glm_lowerCI,
                              KB_GENIC_glm_upperCI,
                              NSEG_NONGENIC_rate,
                              NSEG_NONGENIC_cas_rate,
                              NSEG_NONGENIC_con_rate,
                              NSEG_NONGENIC_cascon_ratio,
                              NSEG_NONGENIC_glm_OR,
                              NSEG_NONGENIC_glm_se,
                              NSEG_NONGENIC_glm_tval,
                              NSEG_NONGENIC_glm_pval,
                              NSEG_NONGENIC_glm_lowerCI,
                              NSEG_NONGENIC_glm_upperCI,
                              KB_NONGENIC_rate,
                              KB_NONGENIC_cas_rate,
                              KB_NONGENIC_con_rate,
                              KB_NONGENIC_cascon_ratio,
                              KB_NONGENIC_glm_OR,
                              KB_NONGENIC_glm_se,
                              KB_NONGENIC_glm_tval,
                              KB_NONGENIC_glm_pval,
                              KB_NONGENIC_glm_lowerCI,
                              KB_NONGENIC_glm_upperCI)

## write to file
write.table(full_data,paste(sep='',wdir,'/',output,'.burden'),col=T,row=F,quo=F,sep='\t')  


