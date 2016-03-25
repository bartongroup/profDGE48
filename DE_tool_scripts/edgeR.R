###############################################################################
# R code to perform edgeR-based differential gene expression analysis
# 
# Run with: Rscript edgeR.R geneExprFile groupsFile outFile normalisation
#
# Author: Chris Cole
###############################################################################

Version = '0.3'
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

## calc means for all genes in each condition as edgeR doesn't give this info
cond1.name = levels(groups$V2)[1]
cond2.name = levels(groups$V2)[2]
cond1 = geneCounts[groups$V2 == cond1.name]
cond2 = geneCounts[groups$V2 == cond2.name]
cond1.means = apply(cond1,1,mean)
cond2.means = apply(cond2,1,mean)

## Start edgeR analysis
d = DGEList(geneCounts,group=groups$V2)

## check what normalisation to do, if any
if (normalise == 'package:default') {
	# do default normalisation for edgeR
	d <- calcNormFactors(d)
   
   # now report normalised read counts at end.
   counts.norm = data.frame(cpm(d, normalized.lib.sizes=T))
   cond1 = counts.norm[groups$V2 == cond1.name]
   cond2 = counts.norm[groups$V2 == cond2.name]
   cond1.means = apply(cond1,1,mean)
   cond2.means = apply(cond2,1,mean)
   
} else {
	# don't do any. Input data is either already normalised or unnormalised.
}

## estimate the common and tagwise dispersion
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

## determine differentially expressed genes (using exact test)
de.com <- exactTest(d)

# need to extract table data from de.com object,
# then select only required columns, plus calculate adjust p-value
all=de.com$table
dat = data.frame(rownames(all),round(all$logFC,digits=2),
								 signif(p.adjust(all$PValue,method='BH'),digits=3),
								 round(cond1.means,digits=2),round(cond2.means,digits=2))
names(dat) <- c('GeneID','log2foldchange','significance',cond1.name,cond2.name)

# optional - plot data with significant genes highlighted
#sig.cut = 0.01
#de.sig = all[all$FDR < sig.cut,]
#plotSmear(d, de.tags = rownames(de.sig), 
#					main = "FC plot using common dispersion")
#abline(h=c(1,-1),col="dodgerblue")
# output significant genes to file
write.table(dat, file=outFile,sep="\t",row.names=F)
