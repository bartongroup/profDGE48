###############################################################################
# R code to perform DESeq2-based DE analysis
#
# Performs standard analysis using the Wald method.
#
# Run with: Rscript deseq2.R geneExprFile groupsFile outFile normalization
#
# Author: Chris Cole
###############################################################################

Version = '1.2'
ver = packageVersion('DESeq2');

## parse arguments and get files - error is they're missing.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
	if (args[1] == '--version') {
		cat(sprintf("deseq2.R version\t%s\n", Version))
		cat(R.version.string,"\n")
		cat(sprintf("DESeq2 version\t%s\n",ver))
		quit('no')
	}
	stop('ERROR - Not enough arguments')
}
library("DESeq2")

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

if("WT" %in% levels(groups$V2)){
  cond1.name <- "WT"
}else{
  cond1.name <- levels(groups$V2)[1]
}
cond2.name <- levels(groups$V2)[which(levels(groups$V2)!=cond1.name)]
groups$V2 <- relevel(groups$V2,cond1.name)


dds <- DESeqDataSetFromMatrix(countData = geneCounts, colData = DataFrame(groups), design = ~ V2)
if (normalise == 'package:default') {

   dds <- DESeq(dds)
   res <- results(dds)

} else {
   ## do DESeq() without estimateSizeFactors(). Set size factors to 1.
   sizeFactors(dds) <- c(rep(1,length(geneCounts)))
   dds <- estimateDispersions(dds)
   dds <- nbinomWaldTest(dds)
}

# get normalised read counts and calc per-condition means for reporting


counts.norm = data.frame(counts(dds, normalized=TRUE))
cond1 = counts.norm[groups$V2 == cond1.name]
cond2 = counts.norm[groups$V2 == cond2.name]
cond1.means = apply(cond1,1,mean)
cond2.means = apply(cond2,1,mean)

# munge data and output
dat = data.frame(rownames(res),
      round(res$log2FoldChange,digits=2),
		signif(res$padj,digits=3),
      round(cond1.means,digits=2),
      round(cond2.means,digits=2))
names(dat) <- c('GeneID','log2foldchange','significance',cond1.name,cond2.name)
write.table(dat, file=outFile,sep="\t",row.names=F)

