###############################################
# R script to perform EBSeq-based DE analysis
#
# Run with: Rscript ebseq.R geneExprFile groupsFile outFile normalisation 
#
# Author: ccole
###############################################################################

Version = '0.3'
ver = packageVersion('EBSeq');

## parse arguments and get files - error is they're missing.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
	if (args[1] == '--version') {
		cat(sprintf("ebseq.R version\t%s\n", Version))
		cat(R.version.string,"\n")
		cat(sprintf("EBSeq version\t%s\n",ver))
		quit('no')
	}
	stop('ERROR - Not enough arguments')
}
library("EBSeq")

countsFile = args[1]
groupsFile = args[2]
outFile = args[3]
normalise = args[4]

## Read gene expression matrix
geneCounts = read.delim(countsFile, header=T, row.names=1)

## Read conditions/groups matrix
groups = read.delim(groupsFile, header=F, row.names=1)

## check they have the numbers of samples
if (length(geneCounts) != nrow(groups)) {
	stop(sprintf("ERROR - no. of samples in '%s' (%d) doesn't agree with those in '%s' (%d)",
	             countsFile, length(geneCounts), groupsFile, nrow(groups)))
}

## strip out zero count genes in all samples
geneCounts = geneCounts[rowSums(geneCounts)>0,]

## start analysis
# normalise libraries
if (normalise == 'package:default') {
   Sizes=MedianNorm(geneCounts)
} else {
   # don't normalise, set factors to 1
   Sizes = c(rep(1, length(geneCounts)))
   names(Sizes) <- names(geneCounts)
}

# do DE test
EBOut=EBTest(Data=as.matrix(geneCounts), Conditions=as.factor(groups$V2), sizeFactors=Sizes,
             maxround=5)
EBDERes=GetDEResults(EBOut, FDR=0.05)
# get posterior probabilities
EBPP=GetPP(EBOut)
# fold-change on normalised data after posterior probabilities
FC=PostFC(EBOut)
log2FC = log2(FC$RealFC)

# get normalised read counts and calc per-condition means for reporting
counts.norm = data.frame(GetNormalizedMat(geneCounts, Sizes))
cond1.name = levels(groups$V2)[1]
cond2.name = levels(groups$V2)[2]
cond1 = counts.norm[groups$V2 == cond1.name]
cond2 = counts.norm[groups$V2 == cond2.name]
cond1.means = apply(cond1,1,mean)
cond2.means = apply(cond2,1,mean)

# munge data and output
dat = data.frame(rownames(geneCounts),
      round(log2FC,digits=2),
	  	signif((1-EBPP),digits=3),
      round(cond1.means,digits=2),
      round(cond2.means,digits=2))
names(dat) <- c('GeneID','log2foldchange','significance',cond1.name,cond2.name)
write.table(dat, file=outFile,sep="\t",row.names=F)

