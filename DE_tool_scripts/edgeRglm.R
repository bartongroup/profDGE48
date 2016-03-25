###############################################################################
# R code to perform edgeR-based differential gene expression analysis
#
# Run with: Rscript edgeR.R geneExprFile groupsFile outFile normalisation
#
# Author: Pieta Schofield
###############################################################################

Version = '1.0'
ver = packageVersion('edgeR');

## parse arguments and get files - error is they're missing.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
	if (args[1] == '--version') {
		cat(sprintf("edgeR.R version\t%s\n", Version))
		cat(R.version.string,"\n")
		cat(sprintf("edgeR version\t%s\n",ver))
		quit('no')
	}
	stop('ERROR - Not enough arguments')
}
library("edgeR")

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

if("WT" %in% levels(groups$V2)){
  cond1.name <- "WT"
}else{
  cond1.name <- levels(groups$V2)[1]
}
cond2.name <- levels(groups$V2)[which(levels(groups$V2)!=cond1.name)]
groups$V2 <- relevel(groups$V2,cond1.name)

## build the DGEList
design <- model.matrix( ~ groups$V2)
y <- DGEList(counts=geneCounts, group=groups$V2)
normVals <- F
if (normalise == 'package:default'){
  # do default normalisation for edgeR
  y <- calcNormFactors(y, method="upperquartile")
  normVals <- T
}
## estimate dispersion
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
## perform DE test using GLM
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
de.com <- topTags(lrt,n=100000)
## calculate conditional means
counts.norm = data.frame(cpm(y, normalized.lib.sizes=T))
cond1.means = rowMeans(counts.norm[groups$V2 == cond1.name])
cond2.means = rowMeans(counts.norm[groups$V2 == cond2.name])
means <- as.data.frame(cbind(cond1.means,cond2.means))
colnames(means) <- c(cond1.name,cond2.name)
# need to extract table data from de.com object,
# then select only required columns, plus calculate adjust p-value
all=de.com$table
dat = data.frame(cbind(GeneID=rownames(all),log2foldchange=round(all$logFC,digits=2),
								 significance=signif(all$FDR,digits=3)))
dat <- merge(dat,means,by.x="GeneID",by.y="row.names")

write.table(dat, file=outFile,sep="\t",row.names=F)
