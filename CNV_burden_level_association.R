
## CNV_burden_level_association.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: January 2016


## ----- DESCRIPTION:

## Interactive R script to read and graph CNV burden results

## - Running CNV burden script
## - specify parameters to view and graph
## - show main results
## - Forest plot


## ----- REQUIREMENTS:

## PLINK and R installed
## PLINK .cfile: .cnv / .fam / .map file (need to do --cnv-make-map in PLINK if not provided)
## .phe file with platform/dataset information
## .cfile and .phe file must have the same name
## gene list file (currently using hg18_refGene_plink.txt)


## ----- SCRIPT SET UP:

## PART 1: Running CNV burden script
## PART 2: Read in burden results and set parameters
## PART 3: Quick view of results
## PART 4: CNV burden KB forest plot


## ----- NOTES:

## This script is currently designed to be run interactively, so there are no input/output arguments. 
## This script is not optimized to handle various datasets, so adjustements may be needed to correctly format results

## Here are the main formatting options to adjust to fit values correctly:

## I created a row with all NAs to separate the PGC results from each dataset. You can remove this in the 'dataset information' setup

## Margins around graph:
# par(oma=c(1,22,1,1),xpd=FALSE)

## width of x-axis in graph is relative to the largest confidence intervals
# xrange <- xlims[2]-xlims[1]

## Margin information parameters
# text(x=rep(xlims[1]-(1*xrange),rows),
#      y=1:rows,as.character(bn1$data_set),
#      pos=4,xpd=NA,cex=1,font=c(2,rep(1,rows-1))) 

## 'rows' is an object I created that has the number of rows graphed in the forest plot
## 'x=' is the relative position outside of the graph for each value, where xlims[1]-(1*xrange) is one graph length to the left of the start of the graph
## 'y=' is the position on the y-axis (1 being the bottom and going up)
## 'cex=' adjusts the size of the text
## 'font=' 1 is regular, 2 is bold





## PART 1: Running CNV burden script


## Example submission to Broad LSF cluster:

bsub -P CNV -q priority -R 'rusage[mem=8]' -o CNV_burden_platform.out Rscript --verbose CNV_burden_platform.R \
/humgen/atgu1/fs03/wip/howrigan/pgc_cnv/github \
CNV_data \
hg18_refGene_plink.txt \
CNV_burden_platform

## NOTE: details of burden script are available at the start the script code



## PART 2: Read in burden results and set parameters

## read in burden results
bn <- read.table('CNV_burden_platform.burden',h=T,stringsAsFactors=F)

## read in phenotype data
phe <- read.table('CNV_data.phe',h=T,stringsAsFactors=F)

## graphing directory (where will my graphs be printed?)
graph_dir <- '/home/unix/howrigan/private_html/pgc_cnv/github'


## dataset information sorted by dataset size
dat <- sort(table(phe$CNV_platform),decreasing=F)
data_set <- c(names(dat),NA,'PGC')
data_size <- c(as.integer(dat),NA,nrow(phe))
cas <- sort(table(phe$CNV_platform[phe$AFF==2]),decreasing=F)
cas_size <- c(as.integer(cas),NA,sum(cas))
con <- sort(table(phe$CNV_platform[phe$AFF==1]),decreasing=F)
con_size <- c(as.integer(con),NA,sum(con))
## You should adjust the colors column to your liking
colors <- c(rep('blue',length(data_set)-2),NA,'black')
data <- cbind.data.frame(data_set,data_size,cas_size,con_size,colors)
data$ordering <- seq(1,nrow(data),1)


## specify parameters to view and graph
type <- c('allCNV','del','dup')
region <- c('allregions','novelregions')
size <- c('Allsize','100kb','100-200kb','200-500kb','500kb') 
freq <- c('Allfreq','singleton','2-10','11plus')


## PART 3: Quick view of results

for (a in 1:length(type)) {
for (b in 1:length(region)) {
for (c in 1:length(size)) {
for (d in 1:length(freq)) {


## Reduce burden file to specific parameters
bn0 <- bn[bn$CNV_type == type[a] &
    bn$region_set == region[b] &
    bn$CNV_size==size[c] & 
    bn$CNV_freq == freq[d],]

options(width=200)

print(bn0[,c('data_set','region_set','CNV_type','CNV_freq','NSEG_rate','plink_AFF_CNV','plink_UNAFF_CNV','plink_CASCON_RATE','NSEG_cas_rate','NSEG_con_rate','NSEG_glm_OR','NSEG_glm_pval')])
cat('\n')
print(bn0[,c('data_set','region_set','CNV_type','CNV_freq','KB_rate','plink_AFF_TOTKB','plink_UNAFF_TOTKB','plink_CASCON_TOTKB','KB_cas_rate','KB_con_rate','KB_glm_OR','KB_glm_pval')])
cat('\n')
print(bn0[,c('data_set','region_set','CNV_type','CNV_freq','COUNT_rate','plink_AFF_GRATE','plink_UNAFF_GRATE','plink_CASCON_GRATE','COUNT_cas_rate','COUNT_con_rate','COUNT_glm_OR','COUNT_glm_pval')])
cat('\n')
print(bn0[,c('data_set','region_set','CNV_type','CNV_freq','NGENE_rate','plink_AFF_GRATE','plink_UNAFF_GRATE','plink_CASCON_GRATE','NGENE_cas_rate','NGENE_con_rate','NGENE_glm_OR','NGENE_glm_pval')])
cat('\n')
print(bn0[,c('data_set','region_set','CNV_type','CNV_freq','NSEG_GENIC_rate','plink_AFF_GPROP','plink_UNAFF_GPROP','plink_CASCON_GPROP','NSEG_GENIC_cas_rate','NSEG_GENIC_con_rate','NSEG_GENIC_glm_OR','NSEG_GENIC_glm_pval')])
cat('\n')
print(bn0[,c('data_set','region_set','CNV_type','CNV_freq','KB_GENIC_rate','KB_GENIC_cas_rate','KB_GENIC_con_rate','KB_GENIC_glm_OR','KB_GENIC_glm_pval')])
cat('\n')
print(bn0[,c('data_set','region_set','CNV_type','CNV_freq','NSEG_NONGENIC_rate','NSEG_NONGENIC_cas_rate','NSEG_NONGENIC_con_rate','NSEG_NONGENIC_glm_OR','NSEG_NONGENIC_glm_pval')])
cat('\n')
print(bn0[,c('data_set','region_set','CNV_type','CNV_freq','KB_NONGENIC_rate','KB_NONGENIC_cas_rate','KB_NONGENIC_con_rate','KB_NONGENIC_glm_OR','KB_NONGENIC_glm_pval')])
cat('\n')

cat(paste0('a=',a,'; b=',b,'; c=',c,'; d=',d,'\n'))

} ## END of d LooP
} ## END of c LooP
} ## END of b LooP
} ## END of a LooP




## PART 4: CNV burden KB forest plot

## Burden type
burden <- c('NSEG','KB','COUNT','NGENE','NSEG_GENIC','KB_GENIC','NSEG_NONGENIC','KB_NONGENIC')
burden_names <- c('CNV','KB','Genes','Genes','CNV','KB','CNV','KB')

# a=1 ; b=1 ; c=1 ; d=1; e=7
# a=1; b=1; c=2; d=4; e=8

for (a in 1:length(type)) {
for (b in 1:length(region)) {
for (c in 1:length(size)) {
for (d in 1:length(freq)) {
for (e in 1:length(burden)) {

## Reduce burden file to specific parameters
bn0 <- bn[bn$CNV_type == type[a] &
    bn$region_set == region[b] &
    bn$CNV_freq == freq[c] &
    bn$CNV_size==size[d],]

## merge with phenotype data
bn1 <- merge(data,bn0,by='data_set',all.x=T)
bn1 <- bn1[order(bn1$ordering,decreasing=T),]
bn1$colors <- as.character(bn1$colors)

## specify burden measurement
con_rate <- eval(parse(text=paste(sep='','bn1$',burden[e],'_con_rate')))
OR <- eval(parse(text=paste(sep='','bn1$',burden[e],'_glm_OR')))
lowerCI <- eval(parse(text=paste(sep='','bn1$',burden[e],'_glm_lowerCI')))
upperCI <- eval(parse(text=paste(sep='','bn1$',burden[e],'_glm_upperCI')))
pval <- eval(parse(text=paste(sep='','bn1$',burden[e],'_glm_pval')))
plot_titles <- paste(c('Platform','Cases','Controls',burden_names[e],'pval','OR'),sep='')

## removing infinite CIs
lowerCI[lowerCI==-Inf] <- NA ; lowerCI[lowerCI==Inf] <- NA
lowerCI[upperCI==-Inf] <- NA ; lowerCI[upperCI==Inf] <- NA
upperCI[upperCI==-Inf] <- NA ; upperCI[upperCI==Inf] <- NA
upperCI[lowerCI==-Inf] <- NA ; upperCI[lowerCI==Inf] <- NA


## ===== GRAPH
pdf(paste(graph_dir,"/CNV_burden_platform_",type[a],"_",region[b],"_",size[c],"_",freq[d],"_",burden[e],"_forestplot.pdf",sep=""),height=8.5,width=11)

op <- par()
par(oma=c(1,22,1,1),xpd=FALSE)

## point sizes
bn1$cexsize <- (bn1$data_size/max(bn1$data_size,na.rm=T))*2
bn1$cexsize[!is.na(bn1$cexsize)] <- bn1$cexsize[!is.na(bn1$cexsize)] + 1.8

#confidence intervals
xlims <- c(min(lowerCI,na.rm=T)-(0.1*max(upperCI,na.rm=T)),
           max(upperCI,na.rm=T)+(0.05*max(upperCI,na.rm=T)))
xrange <- xlims[2]-xlims[1]

## override soft plotting axis in event of non-numerical CIs
if (xrange==-Inf) { xlims <- c(0,2); xrange <- xlims[2]-xlims[1] }

#plot the betas & SE(betas)
plot(OR,seq(1,nrow(bn1),1),axes=F,xlab='Odds Ratio (95% CI)',ylab='',col=bn1$colors,main='',pch=20,cex=0.5,xlim=xlims)

axis(1)

## 95% CI lines
for (x in 2:nrow(bn1)){lines(x=c(lowerCI[x],upperCI[x]),y=c(x,x),col=bn1$colors[x],lwd=3)}

## size oriented Points
points(OR[2:nrow(bn1)],seq(2,nrow(bn1),1),pch=20,col='black',cex=bn1$cexsize[2:nrow(bn1)])

## combined line
diamond.x <- c(lowerCI[1],OR[1],upperCI[1],OR[1],lowerCI[1])
diamond.y <- c(1,1.1,1,0.9,1)
polygon(diamond.x,diamond.y,col='black',border=T)

## Null line
abline(v=1,lty=1,lwd=.5) 

## Combined data line
abline(v=OR[1],lwd=.5,lty=2)


## === Margin Information

rows <- length(bn1$data_set)

## Sample names
text(x=rep(xlims[1]-(1*xrange),rows),
     y=1:rows,as.character(bn1$data_set),
     pos=4,xpd=NA,cex=1,font=c(2,rep(1,rows-1))) 
## cases
text(x=rep(xlims[1]-(0.82*xrange),rows),
     y=1:rows,as.character(bn1$cas_size),
     pos=4,xpd=NA,cex=1,font=c(2,rep(1,rows-1))) 
## controls
text(x=rep(xlims[1]-(0.64*xrange),rows),
     y=1:rows,as.character(bn1$con_size),
     pos=4,xpd=NA,cex=1,font=c(2,rep(1,rows-1))) 
## control burden rate
text(x=rep(xlims[1]-(0.45*xrange),rows),
     y=1:rows,round(con_rate,2),
     pos=4,xpd=NA,cex=1,font=c(2,rep(1,rows-1))) 
## p-values
pvals <- round(pval,3)
pvals[which(pvals < .001)] <- format(pval[which(pvals < .001)],scientific=T,digits=2)
pvals[is.na(con_rate)] <- NA
text(x=rep(xlims[1]-(0.3*xrange),rows),
     y=1:rows,pvals,
     pos=4,xpd=NA,cex=1,font=c(2,rep(1,rows-1))) 
## Odds ratio
text(x=rep(xlims[1]-(0.15*xrange),rows),
     y=1:rows,round(OR,2),
     pos=4,xpd=NA,cex=1,font=c(2,rep(1,rows-1))) 


## --- Titles
xpts1 <- c(xlims[1]-(1*xrange),xlims[1]-(0.82*xrange),xlims[1]-(0.64*xrange),xlims[1]-(0.45*xrange),xlims[1]-(0.3*xrange),xlims[1]-(0.15*xrange))
text(x=xpts1,
     y=rows+0.5,
     labels=plot_titles,
     pos=4,xpd=NA,cex=1,font=2) 


dev.off()

print(paste('a=',a,'; b=',b,'; c=',c,'; d=',d,'; e=',e,sep=''))

} ## END of a LooP
} ## END of b LooP
} ## END of c LooP
} ## END of d LooP
} ## END of e LooP


