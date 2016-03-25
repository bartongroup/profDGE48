#!/sw/opt/R/3.0.2/bin/Rscript --vanilla
###############################################################################
# R code to perform edgeR-based differential gene expression analysis
#
# Run with: Rscript edgeR.R geneExprFile groupsFile outFile normalisation
#
# Author: Pieta Schofield
#         edited 2015-11-10 to work with R-3.2.2
###############################################################################
require(Biobase)
require(baySeq)
require(matrixStats)
cl <- NULL
Version = '1.0.0'
ver = packageVersion('baySeq');

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

countsFile = args[1]
groupsFile = args[2] #Not used as this info is in exprsFile already
outFile = args[3]
normalise = args[4]

if(head(unlist(strsplit(normalise,":")),1)=="package"){
	if(tail(unlist(strsplit(normalise,":")),1)=="default"){
		normalise="total"
	}else{
		normalise=tail(unlist(strsplit(normalise,":")),1)
	}
}else{
	normalise="total"
}
groups <- read.delim(groupsFile,head=F,sep="\t",row.names=1)
if("WT" %in% levels(groups$V2)){
  cond1.name <- "WT"
}else{
  cond1.name <- levels(groups$V2)[1]
}
cond2.name <- levels(groups$V2)[which(levels(groups$V2)!=cond1.name)]
groups$V2 <- relevel(groups$V2,cond1.name)
conds <- groups$V2
#Read the raw counts, remove no count features and sanify columm names
rawCounts <- read.delim(countsFile,h=T,sep="\t",row.names=1)
rawCounts <- rawCounts[rowSums(rawCounts)>0,]
colnames(rawCounts) <- sapply(colnames(rawCounts),
					              function(x){
						              paste0(unlist(strsplit(x,"_"))[1:2],collapse="_")
					              })

#Build a baySeq countData object
CD <- new("countData",data=data.matrix(rawCounts),
					replicates=substr(colnames(rawCounts),1,2),
          groups=list(NDE=rep(1,length(colnames(rawCounts))),
                      DE=conds),
					annotation=data.frame(geneID=rownames(rawCounts),
																row.names=rownames(rawCounts)))

#Run the DE tools in baySeq
libsizes(CD) <- getLibsizes(CD,estimationType=normalise)
CDPriors <- getPriors.NB(CD,
								samplesize=round(length(rownames(CD@annotation))/4,0),
								estimation="QL",cl=NULL)
CDPost <- getLikelihoods(CDPriors,pET="BIC",cl=NULL)

#Get the results and calculate foldchanges
res <- topCounts(CDPost,group="DE",
								 number=length(rownames(CD@annotation)),
								 normal=TRUE)
res$geneID <- rownames(res)
c1Mean <- round(rowMeans(res[,which(grepl(cond1.name,colnames(res)))]),2)
c2Mean <- round(rowMeans(res[,which(grepl(cond2.name,colnames(res)))]),2)
fC <- log2(c2Mean/c1Mean)

#Return data we want
dat <- data.frame(GeneID=res$geneID,
								 log2foldchange=fC,
								 significance=res$FDR,
								 c1Mean,c2Mean)
colnames(dat)[4:5] <- c(cond1.name,cond2.name)
write.table(dat, file=outFile,sep="\t",row.names=F)
sessionInfo()

