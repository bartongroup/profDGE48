library(pvclust)
raw = read.table('gene_ranking_sig.txt', sep="\t", header=TRUE)
nr = dim(raw)[2]
dat = raw[c(4:nr)]
row.names(dat) = raw$gene
clust = pvclust(dat, method.dist="cor", method.hclust="comp", nboot=1000)

pdf("cluster_pv.pdf")
plot(clust)
pvrect(clust, alpha=0.95)
dev.off()

