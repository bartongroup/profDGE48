###############################################################################
# R code to perform SAM-based differential gene expression analysis
# 
# Run with: Rscript SAMseq.R geneExprFile groupsFile outFile
#
# Author: Pieta Schofield
###############################################################################

Version = '0.1.1'
ver = packageVersion('samr');

## parse arguments and get files - error is they're missing.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
	if (args[1] == '--version') {
		cat(sprintf("samr.R version\t%s\n", Version))
		cat(R.version.string,"\n")
		cat(sprintf("samr version\t%s\n",ver))
		quit('no')
	}
	stop('ERROR - Not enough arguments')
}
library("samr")

countsFile = args[1]
groupsFile = args[2]
outFile = args[3]
normalize = args[4]

## Read gene expression matrix
gCounts = read.delim(countsFile, skip=1,header=F, row.names=1)
geneCounts = gCounts[rowSums(gCounts)>0,]

## Read conditions/groups matrix
gGroups = read.delim(groupsFile, header=F, row.names=1)

## normalize 
if(normalize=='deseq'){
     gc.norm=geneCounts
} else if (normalize=='package:default'){
     gc.norm=samr.norm.data(geneCounts)
} 

## Call SAMSeq
samfit <- SAMseq(round(gc.norm),as.numeric(gGroups$V2),
		resp.type="Two class unpaired",
		geneid=rownames(gc.norm),fdr.output=1)

# output the bits we want
dat <- data.frame(GeneID=rownames(samfit$samr.obj$x),
		log2foldchange=log2(samfit$samr.obj$foldchange),
		significance=samr.pvalues.from.perms(samfit$samr.obj$tt,
				samfit$samr.obj$ttstar),
		Score=samfit$samr.obj$tt)

write.table(dat, file=outFile,sep="\t",row.names=F)
sessionInfo()
