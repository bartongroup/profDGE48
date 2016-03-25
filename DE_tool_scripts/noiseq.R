#!/sw/opt/R/3.0.2/bin/Rscript --vanilla
###############################################################################
# R code to perform edgeR-based differential gene expression analysis
# 
# Run with: Rscript edgeR.R geneExprFile groupsFile outFile normalisation
#
# Author: Pieta Schofield
###############################################################################
require(NOISeq)
Version = '0.2.1'
ver = packageVersion('NOISeq');

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

exprsFile = args[1]
groupsFile = args[2] # not used as this info is already in exprsFile
outFile = args[3]
normalise = args[4]

if(head(unlist(strsplit(normalise,":")),1)=="package"){
	if(tail(unlist(strsplit(normalise,":")),1)=="default"){
		normalise="n"
	}else{
		normalise=tail(unlist(strsplit(normalise,":")),1)
	}
}else{
	normalise="n"
}
groups <- read.delim(groupsFile,head=F,sep="\t",row.names=1)

#Read the raw counts remove null features and sanify column names
rawCounts <- read.delim(exprsFile,h=T,sep="\t",row.names=1)
colnames(rawCounts) <- sapply(colnames(rawCounts),
												function(x){
													paste0(unlist(strsplit(x,"_"))[1:2],collapse="_")
												})
rawCounts <- rawCounts[rowSums(rawCounts)>0,]

#Make NOISeq objects
myFactors <- data.frame(treatment=groups[,1], 
												samples=colnames(rawCounts))
myData <- readData(data=rawCounts,factors=myFactors)

#Run the DE test
mynoiseq <- noiseq(myData,k = 0.1, norm = normalise,
												replicates="biological", 
												factor ="treatment",lc = 1,
												conditions=levels(groups[,1]))

# extract results
res <- degenes(mynoiseq,q=0,M=NULL)

# export the bits we want
dat <- data.frame(GeneID=rownames(res),
								 log2foldchange=res$M,
								 significance=res$prob,
								 Snf2=res[,1],WT=res[,2])
colnames(dat)[4:5] <- levels(groups[,1])
write.table(dat, file=outFile,sep="\t",row.names=F)
sessionInfo()

