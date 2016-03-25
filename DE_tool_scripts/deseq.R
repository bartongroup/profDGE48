###############################################################################
# R code to perform DESeq-based DE analysis
#
# Run with: Rscript deseq.R geneExprFile groupsFile outFile normalization
#
# Author: Alexander Sherstnev
###############################################################################

script.version<-'1.6'
script.testcode<-0

##############################################################################
# main DE analysis routine.
deseq_2cond_Nrepl<-function(datafile,phenofile,outfile,normalization="nothing"){
# test code
#library(DESeq)
#datafile<-"test.dat"
#phenofile<-"pheno.dat"
#outfile<-"test.out"
#normalization="internal"

  if(1==script.testcode){
    cat("Data from ",datafile,", pheno metainfo from ",phenofile,
				", results in ",outfile,"\n",sep="")
    cat("Normalization: ",normalization,"\n",sep="")
  }

# read expression data
  data<-read.csv(datafile,header=TRUE,sep="\t",row.names=1,check.names=F)

##############################################
# round data if there are non-integer columns
  intcol<-vector("logical")
  for(i in 1:ncol(data)){
    intcol<-c(intcol,is.integer(data[,i]))
  }
  if (!all(TRUE==intcol)) {
    warning("WARNING! Non-integer expression levels. Data rounded")
    data<-round(data)
  }

##############################################
# remove genes with 0 in all replicates
  if(1==script.testcode) cat("Before selection: Ngenes = ",nrow(data),"\n")
  data.sel<-data[rowSums(data)>0,]
  if(1==script.testcode) cat("After  selection: Ngenes = ",nrow(data.sel),"\n")

##############################################
# form condition name string
  conds<-read.csv(phenofile,header=FALSE,sep="\t",row.names=1,as.is=T)
  conds$name<-rownames(conds)
  conditions<-vector("character",0)
  for(name in colnames(data)){
    conditions<-c(conditions,conds$V2[conds$name==name])
  }
  if(1==script.testcode){
    print("conditions:")
    print(conditions)
  }

##############################################
# create main deseq object
  cds<-newCountDataSet(data.sel,conditions)

##############################################
# apply normalization if needed
  if ("package:default"==normalization) {
    if(1==script.testcode){
      print("Internal DESeq normalization applied")
    }
    cds <- estimateSizeFactors( cds )
  } else {
    if(1==script.testcode) {
      print("pre-normalized data supplied")
    }
    sizeFactors(cds)<-rep(1,length=ncol(data.sel))
  }

  if(1==script.testcode){
# test code
    print("Raw replicate counts: ")
    print(colSums(counts(cds)))
    print("Normalization factors: ")
    print(sizeFactors(cds))
    print("Normalized replicate counts: ")
    print(colSums(counts(cds,norm=T)))
  }

##############################################
# calculate NORMALIZED data means and fold-change
  condnames<-levels(factor(conditions))
  data.sel.norm<-data.frame(counts(cds,norm=T))
  names(data.sel.norm)<-conditions

# for tests only!
#  write.table(data.sel.norm,file="deseq.norm.dat-corrected",sep="\t")
  
  data.sel.norm.cond1<-data.sel.norm[,colnames(data.sel.norm) %in% condnames[1]]
  data.sel.norm.cond2<-data.sel.norm[,colnames(data.sel.norm) %in% condnames[2]]
  cond1.means<-apply(data.sel.norm.cond1,1,mean)
  cond2.means<-apply(data.sel.norm.cond2,1,mean)
  norm.means<-data.frame(GeneID=names(cond1.means),
												 mean.mine.c1=cond1.means,mean.mine.c2=cond2.means)
  norm.means$mylog2foldchange<-log2(norm.means$mean.mine.c2 /
																		   norm.means$mean.mine.c1)

##############################################
# main analysis routines
#  cds<-estimateDispersions(cds,fitType='local',method='per-condition')
  cds<-estimateDispersions(cds)
  res<-nbinomTest(cds,condnames[1],condnames[2])

##############################################
# prepare and write output object
  res.raw<-data.frame(GeneID=res$id,significance=res$padj,
											log2foldchange=res$log2FoldChange,
											mean.c1=res$baseMeanA,mean.c2=res$baseMeanB)
  res.output<-merge(res.raw,norm.means,by='GeneID',all=T)
  res.output<-res.output[,c("GeneID","log2foldchange",
														"significance","mean.c1","mean.c2")]
  if(1==script.testcode){
    DEnum1<-length(res.output$significance[res.output$significance<0.001])
    DEnum2<-length(res.output$significance[res.output$significance<0.01])
    DEnum3<-length(res.output$significance[res.output$significance<0.05])
    print(paste(DEnum1," genes with p-value<0.001",sel=""))
    print(paste(DEnum2," genes with p-value<0.01",sel=""))
    print(paste(DEnum3," genes with p-value<0.05",sel=""))
  }
  write.table(res.output[order(res.output$significance),],file=outfile,
							sep="\t",row.names=FALSE)
}

###############################################################################
library(Biobase)
library(DESeq)

r.version<-R.version.string
instPkgs<-installed.packages()
tool.version<-instPkgs["DESeq","Version"]
biobase.version<-instPkgs["Biobase","Version"]
if(1==script.testcode){
# test code
  cat("\nVersions: R: ",r.version,", Biobase: ",biobase.version,", 
			DESeq: ",tool.version,"\n\n",sep="")
}

##############################################
# default parameter values
datafile<-"data.csv"
phenofile<-"data.pheno"
outfile<-"out.dat"
normalization<-"nothing"

##############################################
# script parameter values
args<-commandArgs(trailingOnly=TRUE)
if(!is.na(args[1])){
  if("--version"==args[1] | "-version"==args[1]){
    cat("deseq.R","\t",script.version,"\n",sep="")
    cat(r.version,"\n",sep="")
    cat("DESeq","\t",tool.version,"\n",sep="")
    q()
  }
  datafile<-args[1]
  if(!file.exists(datafile)){
    stop(paste("Error! data file",datafile,"does not exist"))
  }
}

if(!is.na(args[2])){
  phenofile<-args[2]
  if(!file.exists(phenofile)){
    stop(paste("Error! pheno file",phenofile,"does not exist"))
  }
}
if(!is.na(args[3])){
  outfile<-args[3]
  if(file.exists(outfile)){
    warning(paste("Warning! output file",outfile,
									"exists and will be overwritten"))
  }
}

if(!is.na(args[4])){
  normalization<-args[4]
}

if(1==script.testcode){
# test code
  cat("Expression data file: ",datafile,"\n",sep="")
  cat("Pheno file: ",phenofile,"\n",sep="")
  cat("Output file: ",outfile,"\n",sep="")
  cat("Normalization: ",normalization,"\n\n",sep="")
}

# DE analysis
deseq_2cond_Nrepl(datafile,phenofile,outfile,normalization)
