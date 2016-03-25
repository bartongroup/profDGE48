###############################################################################
# R code to perform limma-based differential gene expression analysis
#
# Run with: Rscript limma.R geneExprFile groupsFile outFile normalisation
#
# Author: Chris Cole
#         edited 2015-11-10 by Pieta Schofield to work with R-3.2.2
###############################################################################

## TODO - fix bug of reporting unnormalised data when doing internal normalisation

Version = '1.0'
ver = packageVersion('limma');

## parse arguments and get files - error if they're missing.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
	if (args[1] == '--version') {
		cat(sprintf("limma.R version\t%s\n", Version))
		cat(R.version.string,"\n")
		cat(sprintf("limma version\t%s\n",ver))
		quit('no')
	}
	stop('ERROR - Not enough arguments')
}
library("limma")

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
	stop(sprintf("ERROR - no. of samples in '%s' (%d) doesn't agree with those in '%s' (%d)", countsFile, length(geneCounts), groupsFile, nrow(groups)))
}

## strip out zero count genes in all samples
geneCounts = geneCounts[rowSums(geneCounts)>0,]

## TODO - what to do with genes with some zeros??
# ignore for time being...

## generate two groups and calc means and fold-changes
## groups$V2 is a factor, so pick the first 'level' (condition) from the factor list (should only be two)
## and select columns which correspond to that condition in the geneCounts.norm data frame
if (length(levels(groups$V2)) != 2) {
	stop(sprintf("ERROR - wrong number of required conditions: %d found where 2 required.",length(levels(groups$V2))))
}

## Do limma analysis
conds<-as.factor(groups$V2)
design<-model.matrix(~conds)
# convert counts to log2-CPM (use edgeR calcNormFactors() or assume data is already normalised)
if (normalise=='package:default') {
	library(edgeR)
	nf <- calcNormFactors(geneCounts)
	y <- voom(geneCounts,design,lib.size=colSums(geneCounts)*nf)
} else if (normalise=='package:deseq') {
	# perform DEseq normalisation instead of default edgeR
	library(DESeq)
	geneCounts.int = round(geneCounts)
	cds <- newCountDataSet(geneCounts.int,conds)
	cds <- estimateSizeFactors(cds)
	y <- voom(geneCounts,design,lib.size=colSums(geneCounts.int)*sizeFactors(cds))
} else {
	# don't do any scaling, but voom() it anyway
	y <- voom(geneCounts,design)
}
fit <- lmFit(y,design)
fit <- eBayes(fit)

## grab all data from topTable and munge columns to return required data
all <- topTable(fit,coef=2,sort.by='none',n=1231231,adjust.method='BH')
all$ID <- rownames(all)
dat <- data.frame(all$ID,round(all$logFC,digits=2),signif(all$adj.P.Val,digits=3))
names(dat) <- c('GeneID','log2foldchange','significance')
write.table(dat,file=outFile,sep="\t",row.names=F)
