## CNV_data_descriptives.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: January 2016


## ----- DESCRIPTION:

## Interactive R script to examine CNV datasets / Platforms and Principal Components 

## This script is a short example of dataset properties to investigate potential confounds, limitations, or outliers 

## Main steps:
## - Checking case/control balance
## - Checking CNV burden by dataset / platform
## - Checking CNV burden by PC

## Additional analyses could include:
##  - CNV burden by sex
##  - Splitting CNV burden into autosomes and X chromosome 
##  - Splitting into large (more reliably called) and small (less reliably called) CNVs
##  - Plate effects (if available)
##  - CNV metrics (if available)


## ----- REQUIREMENTS:

## R and PLINK installed
## PLINK .cfile: .cnv / .fam / .map file (need to do --cnv-make-map in PLINK if not provided)
## .phe file with platform/dataset information
## .cfile and .phe file must have the same name


## ----- SCRIPT SET UP:

## PART 1: Checking case/control balance
## PART 2: Checking CNV burden by dataset / platform
## PART 2.1: Per-platform boxplot of CNV burden by Case/Control status
## PART 3: Checking PCs with platform and CNV burden


## ----- NOTES:

## Set up a subdirectory to put graphs

## This is where I put my graphs on the Broad server - which can be viewed online. Change this line to your specific directory  
graph_dir <- '/home/unix/howrigan/public_html/pgc_cnv/github'



## PART 1: Checking case/control balance

phe <- read.table('CNV_data.phe',h=T,stringsAsFactors=F)

## get case/control count
table(phe$AFF)

## get dataset info (strip first character string from FID)
phe$dataset <- sapply(strsplit(phe$FID,'_'),'[',1)

## get per-dataset platform info
tapply(phe$dataset,phe$CNV_platform,table)

## get per-dataset case/control count
tapply(phe$AFF,phe$dataset,table)

## get proportion of controls in each dataset
tapply(phe$AFF,phe$dataset,function(X) {sum(X[X==1])/length(X)} )

## get proportion of controls in each platform
tapply(phe$AFF,phe$CNV_platform,function(X) {sum(X[X==1])/length(X)} )


## potential issues with case-only dataset in A6.0, but somewhat mitigated by control heavy A6.0 dataset
## No single platform is well-matched for cases and controls





## PART 2: Checking CNV burden by dataset / platform


## get summary information for CNV_data in PLINK
system('plink --noweb --cfile CNV_data --out CNV_data')

indiv <- read.table('CNV_data.cnv.indiv',h=T,stringsAsFactors=F)

## merge relevant columns with .phe file
phe2 <- merge(phe,indiv[,c('FID','NSEG','KB','KBAVG')],by='FID')


## Does CNV burden vary by platform?

## Number of CNV segments
tapply(phe2$NSEG,phe2$CNV_platform,summary)
summary(lm(NSEG ~ as.factor(CNV_platform),data=phe2))

## KB burden
tapply(phe2$KB,phe2$CNV_platform,summary)
summary(lm(KB ~ as.factor(CNV_platform),data=phe2))

## KB length per CNV
tapply(phe2$KBAVG,phe2$CNV_platform,summary)
summary(lm(KBAVG ~ as.factor(CNV_platform),data=phe2))




## PART 2.1: Per-platform boxplot of CNV burden by Case/Control status


## define graphing parameters
nms <- levels(as.factor(phe2$CNV_platform))
cols <- c('blue','red')


## ==== Graph
pdf(paste(graph_dir,'/Platform_CaseCon_boxplot.pdf',sep=''),height=11,width=8.5)

par(mfrow=c(3,1)) ## HARD-CODED

## NSEG plot
boxplot(phe2$NSEG ~ phe2$AFF:phe2$CNV_platform,
    ylab="NSEG",
    col=cols,
    cex=1,
    cex.axis=1)
# legend
legend('topright',legend=c('controls','cases'),fill=c('blue','red'),bty='n',cex=1)


## KB plot
boxplot(phe2$KB ~ phe2$AFF:phe2$CNV_platform,
    ylab="KB",
    col=cols,
    cex=1,
    cex.axis=1)
# legend
legend('topright',legend=c('controls','cases'),fill=c('blue','red'),bty='n',cex=1)


## KBAVG plot
boxplot(phe2$KBAVG ~ phe2$AFF:phe2$CNV_platform,
    ylab="KBAVG",
    col=cols,
    cex=1,
    cex.axis=1)
# legend
legend('topright',legend=c('controls','cases'),fill=c('blue','red'),bty='n',cex=1)


dev.off()



## PART 3: Checking PCs with platform and CNV burden


## PCs provided in .phe file
pcs <- c('C1','C2','C3','C4','C8')


## Run through each PC to see how CNV platform varies
for(a in 1:length(pcs)) {
    cat(pcs[a],'\n\n')
    print(eval(parse(text=paste("tapply(phe2$",pcs[a],",phe2$CNV_platform,mean)",sep=""))))
    print(eval(parse(text=paste("summary(lm(",pcs[a]," ~ as.factor(CNV_platform),data=phe2))",sep=""))))
    cat('\n\n\n\n\n')
}

## Run through each PC to see how CNV count (NSEG) varies
for(a in 1:length(pcs)) {
    cat(pcs[a],'\n\n')
    print(eval(parse(text=paste("summary(lm(",pcs[a]," ~ NSEG,data=phe2))",sep=""))))
    cat('\n\n\n\n\n')
}

## Run through each PC to see how CNV count (NSEG) varies controlling for platform
for(a in 1:length(pcs)) {
    cat(pcs[a],'\n\n')
    print(eval(parse(text=paste("summary(lm(",pcs[a]," ~ NSEG + as.factor(CNV_platform),data=phe2))",sep=""))))
    cat('\n\n\n\n\n')
}

## Run through each PC to see how CNV count (NSEG) varies controlling for dataset
for(a in 1:length(pcs)) {
    cat(pcs[a],'\n\n')
    print(eval(parse(text=paste("summary(lm(",pcs[a]," ~ NSEG + as.factor(dataset),data=phe2))",sep=""))))
    cat('\n\n\n\n\n')
}

## PC2 looks to have an effect on CNV count independent of CNV platform (less so for dataset)



## PART 3.1: Plot of PC and CNV burden

## ==== Graph

## define graphing parameters
cnv.colors <- rainbow(length(levels(as.factor(phe2$CNV_platform))))
col.indx <- cnv.colors[1:nlevels(as.factor(phe2$CNV_platform))] # number of colors
cols <- col.indx[as.numeric(as.factor(phe2$CNV_platform))] # apply to datasets

pdf(paste(graph_dir,'/PC2_NSEG_plot.pdf',sep=''),height=10,width=10)

plot(phe2$C2,
    phe2$NSEG,
    col=cols,
    cex=1,
    pch=20)

# legend
legend("topleft",bty='n',legend=nms,title='CNV platform',col=col.indx,pch=20,bg='white',cex=1)

dev.off()















