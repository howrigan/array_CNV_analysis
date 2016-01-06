## CNV_SNP_level_association.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: January 2016


## ----- DESCRIPTION:

## Interactive Script to generate CNV SNP level association results

## - Case/control association: CNV_CaseCon_assoc.R
## - Residual phenotype generation
## - Residual phenotype association: CNV_assoc_Quantitative.R
## - Family-wise error correction
## - Manhattan Plot
## - UCSC .bed and .bedGraph files

## ----- SCRIPT SET UP:

## PART 1: Case/control association
## PART 2: Residual phenotype generation
## PART 2.1: Residual phenotype association
## PART 3: Z-score association and family-wise error correction
## PART 4: Manhattan Plot
## PART 5: UCSC .bed file
## PART 5.1: UCSC .bedGraph file


## ----- NOTES:




## PART 1: Case/control association

## example command for running CNV_assoc_CaseCon.R on server
Rscript CNV_assoc_CaseCon.R \
/humgen/atgu1/fs03/wip/howrigan/pgc_cnv/github \
"bsub -P PGC_CNV -q hour -R 'rusage[mem=4]' -o CNV_assoc_CaseCon.out" \
CNV_data \
1000 \
CNV_data \
no-save

## output files in /assoc directory: CNV_data_cnv.cnv.summary and CNV_data_cnv.cnv.summary.mperm



## PART 2: Residual phenotype generation

phe <- read.table('CNV_data.phe',h=T,stringsAsFactors=F)
phe$aff <- phe$AFF - 1

## Extract Pearson residuals from logistic regression with relevant covariates
resid_pheno <- glm(aff ~ SEX + CNV_platform + C1 + C2 + C3 + C4 + C8, data=phe, family=binomial)
summary(resid_pheno)
phe$AFF_residual <- resid(resid_pheno,'pearson')

summary(phe$AFF_residual)
summary(phe$AFF_residual[phe$AFF==2])
summary(phe$AFF_residual[phe$AFF==1])

## examine residuals for outliers. In the PGC SCZ data, we removed a few individuals with values > 3 or < -3

write.table(phe,'CNV_data_resid.phe',col=T,row=F,quo=F,sep='\t')




## PART 2.1: Residual phenotype association

## example command for running CNV_assoc_Quantitative.R on server
Rscript CNV_assoc_Quantitative.R \
/humgen/atgu1/fs03/wip/howrigan/pgc_cnv/github \
"bsub -P PGC_CNV -q hour -R 'rusage[mem=4]' -o CNV_assoc_Quantitative.out" \
CNV_data \
CNV_data_resid \
1000 \
CNV_data \
save

## output files in /assoc directory: CNV_data_cnv.cnv.qt.summary and CNV_data_cnv.cnv.qt.summary.mperm
## saved permutation file in /assoc directory: CNV_data_cnv.mperm.dump.all




## PART 3: Z-score association and family-wise error correction


## NOTE: this script uses the .mperm.dump.all files from the quantitative association script

type <- c('cnv','del','dup')
pnum <- 1000
z_pval_5pct_FWER <- NA


for (a in 1:length(type)) {

    ## association
    summcc <- read.table(paste("assoc/CNV_data_",type[a],".cnv.summary",sep=""),h=T)
    mpermcc <- read.table(paste("assoc/CNV_data_",type[a],".cnv.summary.mperm",sep=""),h=T)    
    colnames(mpermcc)[3] <- 'casecon_EMP1'
    colnames(mpermcc)[4] <- 'casecon_EMP2'
    summqt <- read.table(paste("assoc/CNV_data_",type[a],".cnv.qt.summary",sep=""),h=T)
    mpermqt <- read.table(paste("assoc/CNV_data_",type[a],".cnv.qt.summary.mperm",sep=""),h=T)
    colnames(mpermqt)[4] <- 'quant_EMP1'
    colnames(mpermqt)[5] <- 'quant_EMP2'

    assoc <- cbind.data.frame(summcc,summqt[,4:6],mpermcc[,3:4],mpermqt[,4:5])

    ## permutation
    pp <- scan(paste("assoc/CNV_data_",type[a],".mperm.dump.all",sep=""),nmax=((nrow(assoc)+1)*pnum))
    perm <- t(matrix(pp,ncol=nrow(assoc)+1,byrow=TRUE))
    perm <- perm[-1,2:ncol(perm)] #remove counter

    ## remove sites with NCNV==0
    filter <- assoc$NCNV==0
    assoc <- assoc[!(filter),]
    perm <- perm[!(filter),]

    perms <- ncol(perm) ## get permutation size

    ## get permuted variance
    perm_sd <- NA

    for (i in 1:nrow(assoc)){
      pt_perm <- perm[i,]
      pt_diff <- pt_perm-assoc$M0[i]
      perm_sd[i] <- sd(pt_diff)
      # print(i)
    } ## END of i LooP

    assoc$perm_sd <- perm_sd

    ## convert to empirical z-score
    assoc$z <- (assoc$M1-assoc$M0)/assoc$perm_sd
    assoc$z_risk_pval <- pnorm(assoc$z,0,1,lower.tail=F)
    assoc$z_pval <- NA
    assoc$z_pval[assoc$z >= 0] <- 2*pnorm(assoc$z[assoc$z >= 0],0,1,lower.tail=F)
    assoc$z_pval[assoc$z < 0] <- 2*pnorm(assoc$z[assoc$z < 0],0,1,lower.tail=T)

    write.table(assoc,paste("assoc/CNV_Zscore_",type[a],"_",perms,"perm.assoc",sep=""),col=T,row=F,quo=F,sep='\t')

    ## get 95% quantile across permutations

    ## Get min p-val in each permutation
    z_min_pval <- NA
    
    for (i in 1:perms){
      single_perm <- perm[,i]
      pt_diff <- single_perm-assoc$M0
      z <- pt_diff/perm_sd
      z_pval <- NA
      z_pval[z >= 0] <- 2*pnorm(z[z >= 0],0,1,lower.tail=F)
      z_pval[z < 0] <- 2*pnorm(z[z < 0],0,1,lower.tail=T)
      z_min_pval[i] <- min(z_pval)      
      # print(paste('perm',i))
    } ## END of i LooP

    ## get 95% quantile of max Z-scores
    z_pval_5pct_FWER[a] <- as.numeric(quantile(z_min_pval,0.05))

    print(a)    
      
 } ## END of a LooP 

## combine 
data <- cbind.data.frame(type,z_pval_5pct_FWER)
write.table(data,paste("assoc/CNV_Zscore_",perms,"perm_FWER.txt",sep=""),col=T,row=F,quo=F,sep='\t')






## PART 4: Manhattan Plot


## graphing directory (where will my graphs be printed?)
graph_dir <- '/home/unix/howrigan/public_html/pgc_cnv/github'

## CNV types and permutation amount
type <- c('cnv','del','dup')
perms <- 999


## read in family-wise error correction
fwer <- read.table('assoc/CNV_Zscore_999perm_FWER.txt',h=T)

## load graphing resources (may need to install ggplot2 using install.packages('ggplot2'))
library(ggplot2)
source('Manhattan.R')



for (a in 1:length(type)) {


  ## read in association results
  res <- read.table(paste("assoc/CNV_Zscore_",type[a],"_",perms,"perm.assoc",sep=""),h=T,stringsAsFactors=F)
  ## restrict to SNPs with >= 6 overlapping CNV
  res2 <- subset(res,res$NCNV >= 6)
  ## Select only relevant columns
  res3 <- res2[,c('CHR','SNP','BP','z_pval')]
  colnames(res3)[4] <- 'P'

  ## chromosome starts (use all available SNPs)
  chr.starts <- min(res$BP[res$CHR==1])
  for (i in 2:23) { chr.starts[i] <- min(res$BP[res$CHR==i]) }
  ## chromosome ends 
  chr.ends <- max(res$BP[res$CHR==1])
  for (i in 2:23) { chr.ends[i] <- max(res$BP[res$CHR==i]) }
  ## chromosome adder
  adder <- c(0,cumsum(as.numeric(chr.ends)))

  ## formatting manhattan plot for sparse values (i.e. when chromosomes don't have qualified SNPs)
  ## create CHR start and endpoints
  CHR <- seq(1,23,1)
  P <- rep(1,23)
  chr.start.snps <- cbind.data.frame(CHR,chr.starts,chr.starts,P)
  names(chr.start.snps) <- c('CHR','SNP','BP','P')
  chr.end.snps <- cbind.data.frame(CHR,chr.ends,chr.ends,P)
  names(chr.end.snps) <- c('CHR','SNP','BP','P')
  chr.snps <- rbind(chr.start.snps,chr.end.snps)

  ## append to results 
  result <- rbind(res3,chr.snps)

  ## format FWER correction
  fwer_pval <- fwer$z_pval_5pct_FWER[a]
  fwer_pval2 <- round(fwer_pval,3)
  fwer_pval2[which(fwer_pval < .001)] <- format(fwer_pval[which(fwer_pval < .001)],scientific=T,digits=2)

  ## get adjusted ymax for both p-values and FWER correction
  ymax_val <- max(c(-log10(fwer_pval),-log10(min(result$P))))

  ## ==== GRAPH
  pdf(paste(sep="",graph_dir,"/CNV_association_",type[a],"_Manhattan.pdf"),height=6,width=12)

  par(mfrow=c(1,1))

  ## PLOT
  manhattan(result,pch=16,cex=0.6,genomewideline=-log10(fwer_pval),colors=c("lightskyblue2","midnightblue"),suggestiveline=FALSE,cex.axis=0.8,main='',xaxt='n',ymax=ymax_val+1)

  axis(side = 1, at = adder, labels = FALSE, tck = -0.01)
  axis(side = 1, at = mid.vec(adder)[seq(1,length(adder),2)], labels = c(seq(1,22),'X')[seq(1,length(adder),2)], tck = 0,cex.axis=0.75)
  axis(side = 1, at = mid.vec(adder)[seq(2,length(adder),2)], labels = c(seq(1,22),'X')[seq(2,length(adder),2)], tck = 0,cex.axis=0.75)

  legend('topleft',legend=c(paste(type[a],' FWER correction = ',fwer_pval2,sep='')),cex=0.75,lty=1,lwd=2,col=c('red'),bty='n')


  dev.off()

  print(a)

} ## END of a LooP






## PART 5: UCSC .bed and .bedGraph file


## creating .bed file of CNVs
system("plink --noweb --cfile CNV_data --cnv-track --out CNV_data_UCSC_hg18")

## creating .bedGraph file of Z-score results
type <- c('cnv','del','dup')
perms <- 999

for (i in 1:length(type)) {

    assoc <- read.table(paste('assoc/CNV_Zscore_',type[i],'_',perms,'perm.assoc',sep=''),h=T)
    summ <- read.table(paste('assoc/CNV_data_',type[i],'.cnv.summary',sep=''),h=T)

    assoc <- merge(assoc,summ[,c('SNP','BP','CHR')],by='SNP',all=T) ## Add in NCNV=0 spots
    assoc <- assoc[order(assoc$CHR.y,assoc$BP.y),]

    chr <- paste('chr',assoc$CHR.y,sep='')
    chr[chr=='chr23'] <- 'chrX'
    bp1 <- assoc$BP.y
    bp2 <- c(bp1[2:nrow(assoc)],bp1[nrow(assoc)]+2)
    bp2 <- bp2-1
    bp2[bp1 > bp2] <- bp1[bp1 > bp2]+1
    log10pval <- -log10(assoc$z_pval) ## choose which p-value to port into the .bedGraph

    dat <- cbind.data.frame(chr,bp1,bp2,log10pval)
    dat$log10pval[is.na(dat$log10pval)] <- 0
    write.table(dat,paste('assoc/CNV_data_UCSC_hg18_Zpval_',type[i],'.cnv.bedGraph',sep=''),col=F,row=F,quo=F,sep='\t')

    print(i)

} ## END of i LooP

## ADD .bedGraph headers 
write(x="track type=bedGraph name=ALL description=\"CNV Zscore pval\" visibility=full autoScale=off viewLimits=0:7 color=0,100,0 altColor=0,100,200 priority=5 graphType=bar",file="assoc/bedGraph_header_CNV.txt")
write(x="track type=bedGraph name=DEL description=\"DEL Zscore pval\" visibility=full autoScale=off viewLimits=0:7 color=100,0,0 altColor=0,100,200 priority=6 graphType=bar",file="assoc/bedGraph_header_DEL.txt")
write(x="track type=bedGraph name=DUP description=\"DUP Zscore pval\" visibility=full autoScale=off viewLimits=0:7 color=0,0,100 altColor=0,100,200 priority=7 graphType=bar",file="assoc/bedGraph_header_DUP.txt")


## Combine headers and .bedGraph results into a single file
system("cat assoc/bedGraph_header_CNV.txt assoc/CNV_data_UCSC_hg18_Zpval_cnv.cnv.bedGraph assoc/bedGraph_header_DEL.txt assoc/CNV_data_UCSC_hg18_Zpval_del.cnv.bedGraph assoc/bedGraph_header_DUP.txt assoc/CNV_data_UCSC_hg18_Zpval_dup.cnv.bedGraph > CNV_data_UCSC_hg18_Zpval.cnv.bedGraph")



























