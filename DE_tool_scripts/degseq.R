###############################################################################
# R code to perform DEGseq-based differential gene expression analysis
#
# Run with: Rscript DEGseq.R geneExprFile groupsFile outFile normalise
#
# Author: Chris Cole
###############################################################################

Version = '0.5'
ver = packageVersion('DEGseq');

## parse arguments and get files - error is they're missing.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
	if (args[1] == '--version') {
		cat(sprintf("DEGseq.R version\t%s\n", Version))
		cat(R.version.string,"\n")
		cat(sprintf("DEGseq version\t%s\n",ver))
		quit('no')
	}
	stop('ERROR - Not enough arguments')
}
#suppressPackageStartupMessages(suppressWarnings(library(tcltk)))
library("DEGseq")

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
	stop(sprintf("ERROR - no. of samples in '%s' (%d) doesn't agree with those in '%s' (%d)",countsFile, length(geneCounts), groupsFile, nrow(groups)))
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

## Do DEGseq stuff
## needs to write out result to file by default - provide a tmpdir.
## This is deleted by default when R quits.
outputDir <- file.path(tempdir(), "DEGexpExample")
## needs a gene expression matrix - create it
geneExpMat = as.matrix(cbind(rownames(geneCounts),geneCounts))
## determine the column indices which make up each condition/group
expCol1 = which(colnames(geneExpMat) %in%
								colnames(geneCounts[,groups==cond1.name]))
expCol2 = which(colnames(geneExpMat) %in%
								colnames(geneCounts[,groups==cond2.name]))

## perform diffex now (normalisation is currently identical under all
## conditions as the recommended method is always 'none').
if (normalise=='package:default') {
   cat("Selected default normalisation method: none\n")
	DEGexp(geneExpMatrix1=geneExpMat, geneCol1=1, expCol1=expCol1,
				 groupLabel1=cond1.name, geneExpMatrix2=geneExpMat,
				 geneCol2=1, expCol2=expCol2, groupLabel2=cond2.name,
				 method='MARS', outputDir=outputDir, normalMethod='none')
} else {
   cat("Defined normalisation method as:",normalise, fill=T)
	DEGexp(geneExpMatrix1=geneExpMat, geneCol1=1, expCol1=expCol1,
				 groupLabel1=cond1.name, geneExpMatrix2=geneExpMat, geneCol2=1,
				 expCol2=expCol2, groupLabel2=cond2.name, method='MARS',
				 outputDir=outputDir, normalMethod='none')
}
## as DEGseq outputs to file, need to read in data and select only
## required columns
res = read.delim(sprintf("%s/output_score.txt",outputDir),header=T)
## choosing normalised log2 column and BH corrected p-values for output,
## plus also Storey et al 2003 Q-value and Z-score.
## am reporting 'value1' and 'value2' as the gene expression for each condition.

dat=data.frame(res[,1],round(res[,5],digits=4),signif(res[,9],digits=5),
							 round(res[,2],digits=5), round(res[,3],digits=5),
							 signif(res[,10],digits=5),round(res[,6],digits=4))
names(dat) <- c('GeneID','log2foldchange','significance',
								cond1.name,cond2.name,'Qvalue','Zscore')
write.table(dat, file=outFile,sep="\t",row.names=F)

