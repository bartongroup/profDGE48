###############################################################################
# R code to perform PoissonSeq-based differential gene expression analysis
# 
# Run with: Rscript poissonSeq.R geneExprFile groupsFile outFile normalisation
#
# Author: Chris Cole
###############################################################################

Version = '0.3'
ver = packageVersion('PoissonSeq');

## parse arguments and get files - error is they're missing.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
	if (args[1] == '--version') {
		cat(sprintf("poissonSeq.R version\t%s\n", Version))
		cat(R.version.string,"\n")
		cat(sprintf("PoissonSeq version\t%s\n",ver))
		quit('no')
	}
	stop('ERROR - Not enough arguments')
}
library("PoissonSeq")

countsFile = args[1]
groupsFile = args[2]
outFile = args[3]
normalise = args[4] 
# does nothing as 'normalisation' is done by default and can't be 
# controlled without a major edit of the code.

## Read gene expression matrix
geneCounts = read.delim(countsFile, header=T, row.names=1)

## Read conditions/groups matrix
groups = read.delim(groupsFile, header=F,sep="\t", row.names=1)

## check they have the numbers of samples
if (length(geneCounts) != nrow(groups)) {
	stop(paste0("ERROR - no. of samples in ",countsFile, "(",
						 	nrow(geneCounts),") doesn't agree with those in ",
						 	groupsFile," (", nrow(groups),")"))
}

## strip out zero count genes in all samples
geneCounts = geneCounts[rowSums(geneCounts)>0,]

# create required data structure required by PoissonSeq
dat = list()
dat[['n']] = as.matrix(geneCounts)
dat[['y']] = as.numeric(groups[,1])
dat[['type']] = 'twoclass'
dat[['pair']] = FALSE
dat[['gname']] = rownames(geneCounts)

# perform DE analysis.
para = list()

## NB: it's very unclear in the docs and paper how normalisation is controlled,
##     but this seems to be the most obvious place. 
##     It appears 'transform' == 'normalise'
if (normalise == 'package:default') {
	cat("Package default normalisation used: Estimatation of sequencing ",
			"depth based on Poisson goodness-of-fit statistic.", fill=T)
} else {
	cat("Although '",normalise,"' normalisation option has been given, ",
			"it will be ignored and PoissonSeq's default normalisation will be ",
			"used instead.\n", fill=F)
}
para[['ct.sum']] = 0  # force tool to use all genes with sum > 0 (default is 5)
para[['ct.mean']] = 0 # force tool to use all genes with mean > 0 (default is 0.5)
para[['npermu']] = 10000 # number of permutations for p-val calc (default is 100)

res=PS.Main(dat,para)
dat=data.frame(res$gname,
							 round(res$log.fc,digits=4),
							 signif(res$fdr,digits=5),
							 signif(res$pval,digits=5),
							 round(res$tt,digits=4))
names(dat) <- c('GeneID','log2foldchange','significance','rawPval','tt_score')
write.table(dat, file=outFile,sep="\t",row.names=F)
sessionInfo()

