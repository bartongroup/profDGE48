###############################################################################
# R code to perform t-test-based differential gene expression analysis
# 
# Assumes log normal distribution
# 
# Run with: Rscript t-test_equal.R geneExprFile groupsFile outFile
#
# Author: Chris Cole
###############################################################################

## function to count the number >0 values present in a vector
countNonZeros = function(x) {
	n = 0
	for (i in 1:length(x)) {
		if (x[i] > 0) {
			n=n+1
		}
	}
	return(n)
}


Version = '0.5'
ver = packageVersion('stats');

## parse arguments and get files - error is they're missing.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
	if (args[1] == '--version') {
		cat(sprintf("t-test_equal.R version\t%s\n", Version))
		cat(R.version.string,"\n")
		cat(sprintf("stats version\t%s\n",ver))
		quit('no')
	}
	stop('ERROR - Not enough arguments')
}

countsFile = args[1]
groupsFile = args[2]
outFile = args[3]

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

## generate two groups and calc means and fold-changes
## groups$V2 is a factor, so pick the first 'level' (condition) from the factor list (should only be two)
## and select columns which correspond to that condition in the geneCounts.norm data frame
if (length(levels(groups$V2)) != 2) {
	stop(sprintf("ERROR - wrong number of required conditions: %d found where 2 required.",length(levels(groups$V2))))
}
cond1.name = levels(groups$V2)[1]
cond2.name = levels(groups$V2)[2]
cond1 = geneCounts[groups$V2 == cond1.name]
cond2 = geneCounts[groups$V2 == cond2.name]

## In order for t-test to not bork with too many non-zero counts remove genes without enough data.
## For each gene in each condition count how many replicates have >0 counts.
## Only retain genes which have at least three non-zero counts per condition.
## TODO - this is still problematic. Sometimes get 'data is essentially the same' errors in t.test()
##      - may need to be more draconian and remove all genes with zero-count data
cond1.expr = rownames(cond1[apply(cond1,1,countNonZeros) > 2,])
cond2.expr = rownames(cond2[apply(cond2,1,countNonZeros) > 2,])
geneCounts.sub = subset(geneCounts, rownames(geneCounts) %in% intersect(cond1.expr,cond2.expr))

cond1 = geneCounts.sub[groups$V2 == cond1.name]
cond2 = geneCounts.sub[groups$V2 == cond2.name]
cond1[cond1==0] <- NA
cond2[cond2==0] <- NA
cat(sprintf("%d genes removed for having too many zero counts in one or other condition\n",nrow(geneCounts) - nrow(geneCounts.sub)))

# calc means and fold-change
cond1.means = apply(cond1,1,mean)
cond2.means = apply(cond2,1,mean)
FC=cond2.means/cond1.means

## perform t-test
res.t = vector()
for (i in 1:nrow(geneCounts.sub)) {
	#print(i)
	r=t.test(log(cond1[i,]),log(cond2[i,]),var.equal=T)
	res.t[i] = r$p.value
}
# Adjust the p-values with Bonferonni (others are available) and add mean and fold-change data.
res.t.test = data.frame(rownames(geneCounts.sub),log2(FC),p.adjust(res.t,method="BH"),cond1.means,cond2.means)
res.t.test = signif(data.frame(log2(FC),p.adjust(res.t,method="BH"),cond1.means,cond2.means),digits=3)
res.t.test = cbind(rownames(geneCounts.sub),res.t.test)
names(res.t.test) <- c('GeneID','log2foldchange','significance',cond1.name,cond2.name)
write.table(res.t.test,file=outFile,sep="\t",row.names=F)
