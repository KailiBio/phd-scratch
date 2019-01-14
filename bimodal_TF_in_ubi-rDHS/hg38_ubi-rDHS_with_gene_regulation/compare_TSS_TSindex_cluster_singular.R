
# -- Kaili
# This script is for comparing TSS TSindex between TSS-cluster and singular-TSS.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/")
library(ggplot2)

ts_index = read.table("GRCh38_mergedTSS_TSindex_cluster_singular.txt", row.names = 2)
colnames(ts_index) = c("gene", "ts_index", "type")

a = ts_index[ts_index$type=="cluster",]$ts_index
b = ts_index[ts_index$type=="singular",]$ts_index
wilcox.test(a,b)$p.value
t.test(a,b)$p.value

ggplot(ts_index, aes(x=type, y=ts_index, fill=type)) + geom_boxplot(width = 0.5) + 
  labs(title="TSS-cluster VS. singular-TSS\nTissue-Specificity index",
       subtitle="t-test: p-value < 2.2e-16\nwilcoxon: p-value < 2.2e-16")
ggsave("hg38_mergedTSS_TSindex_cluster_singular_boxplot.pdf")

ggplot(ts_index, aes(x=ts_index, color=type)) + geom_density(size=1.2) +
  labs(title="TSS-cluster VS. singular-TSS\nTissue-Specificity index",
       subtitle="t-test: p-value < 2.2e-16\nwilcoxon: p-value < 2.2e-16")
ggsave("hg38_mergedTSS_TSindex_cluster_singular_densityLine.pdf")


###########
# ubi-rOCR overlapped
###########

ubi = as.vector(read.table("GRCh38_ubi-rOCR_overlapped_merged-TSS_list.txt")[,1])
matrix = ts_index[ubi,]

a = matrix[matrix$type=="cluster",]$ts_index
b = matrix[matrix$type=="singular",]$ts_index
wilcox.test(a,b)$p.value
t.test(a,b)$p.value

ggplot(matrix, aes(x=type, y=ts_index, fill=type)) + geom_boxplot(width = 0.5) + 
  labs(title="ubi-rOCR overlapped\nTSS-cluster VS. singular-TSS\nTissue-Specificity index",
       subtitle="t-test: p-value < 2.2e-16\nwilcoxon: p-value < 2.2e-16")
ggsave("hg38_ubi-rOCR_mergedTSS_TSindex_cluster_singular_boxplot.pdf")

ggplot(matrix, aes(x=ts_index, color=type)) + geom_density(size=1.2) +
  labs(title="ubi-rOCR overlapped\nTSS-cluster VS. singular-TSS\nTissue-Specificity index",
       subtitle="t-test: p-value < 2.2e-16\nwilcoxon: p-value < 2.2e-16")
ggsave("hg38_ubi-rOCR_mergedTSS_TSindex_cluster_singular_densityLine.pdf")

###########
# each TSS in TSS-clusters
###########
each_TSS = data.frame(read.table("hg38_each_TSS_in_mergedTSS_cluster_TSindex.bed"))
colnames(each_TSS) = c("TSSid","ts_index","TSS_cluster","TSindex2","gene")

### scatter plot
pdf("hg38_TSS-cluster_each-TSS_TSindex_scatter.pdf")
smoothScatter(each_TSS$ts_index, each_TSS$TSindex2, xlim=c(0.4,1), ylim=c(0.4,1), 
              xlab = "TSindex of each TSS in TSS-clusters", ylab = "TSindex of TSS-cluster",
              main = "comparing the TSindex of TSS-clusters with each TSS in TSS-clusters")
lines(c(0,1), c(0,1), col="red", lwd=2,lty=2)
dev.off()


## boxplot & density
tmp = each_TSS[,c(5,2)]
tmp = transform(tmp, type="each_TSS")
rownames(tmp) = each_TSS$TSSid

matrix = rbind(ts_index, tmp)

a = matrix[matrix$type=="singular",]$ts_index
b = matrix[matrix$type=="cluster",]$ts_index
c = matrix[matrix$type=="each_TSS",]$ts_index
wilcox.test(a,b)$p.value
t.test(a,b)$p.value
wilcox.test(a,c)$p.value
t.test(a,c)$p.value
wilcox.test(b,c)$p.value
t.test(b,c)$p.value

ggplot(matrix, aes(x=type, y=ts_index, fill=type)) + geom_boxplot(width = 0.5) + 
  labs(title="singular-TSS VS. TSS-cluster VS. each-TSS in TS-cluster\nTissue-Specificity index",
       subtitle = "t-test & wilcoxon siginificnat") +
  scale_x_discrete(limits = c("singular", "cluster", "each_TSS"))
ggsave("hg38_mergedTSS_TSindex_cluster_singular_each_boxplot.pdf")

ggplot(matrix, aes(x=ts_index, color=type)) + geom_density(size=1.2) +
  labs(title="singular-TSS VS. TSS-cluster VS. each-TSS in TS-cluster\nTissue-Specificity index",
       subtitle = "t-test & wilcoxon siginificnat")
ggsave("hg38_mergedTSS_TSindex_cluster_singular_each_densityLine.pdf")

###########
# ubi-rOCRs overlapped each TSS in TSS-clusters
###########
ubi_each_TSS = each_TSS[each_TSS$TSS_cluster %in% ubi,]

### scatter plot
pdf("hg38_ubi-rOCRs_TSS-cluster_each-TSS_TSindex_scatter.pdf")
smoothScatter(ubi_each_TSS$ts_index, ubi_each_TSS$TSindex2, xlim=c(0.4,1), ylim=c(0.4,1), 
              xlab = "TSindex of each TSS in TSS-clusters", ylab = "TSindex of TSS-cluster",
              main = "comparing the TSindex of ubi-rOCRs overlapped\nTSS-clusters with each TSS in TSS-clusters")
lines(c(0,1), c(0,1), col="red", lwd=2,lty=2)
dev.off()

ubi_matrix = rbind(ts_index[ubi,], transform(ubi_each_TSS[,c(5,2)],type="each_TSS"))

a = ubi_matrix[ubi_matrix$type=="singular",]$ts_index
b = ubi_matrix[ubi_matrix$type=="cluster",]$ts_index
c = ubi_matrix[ubi_matrix$type=="each_TSS",]$ts_index
wilcox.test(a,b)$p.value
t.test(a,b)$p.value
wilcox.test(a,c)$p.value
t.test(a,c)$p.value
wilcox.test(b,c)$p.value
t.test(b,c)$p.value

ggplot(ubi_matrix, aes(x=type, y=ts_index, fill=type)) + geom_boxplot(width = 0.5) + 
  labs(title="ubi-rOCR overlapped\nsingular-TSS VS. TSS-cluster VS. each-TSS in TS-cluster\nTissue-Specificity index",
       subtitle = "t-test & wilcoxon siginificnat") +
  scale_x_discrete(limits = c("singular", "cluster", "each_TSS"))
ggsave("hg38_ubi-rOCRs_mergedTSS_TSindex_cluster_singular_each_boxplot.pdf")

ggplot(ubi_matrix, aes(x=ts_index, color=type)) + geom_density(size=1.2) +
  labs(title="ubi-rOCR overlapped\nsingular-TSS VS. TSS-cluster VS. each-TSS in TS-cluster\nTissue-Specificity index",
       subtitle = "t-test & wilcoxon siginificnat")
ggsave("hg38_ubi-rOCRs_mergedTSS_TSindex_cluster_singular_each_densityLine.pdf")
###########
# max TSindex of TSS in TSS_clusters
###########
max = read.table("hg38_each_TSS_in_mergedTSS_cluster_TSindex_max.bed", row.names = 1)
max = transform(max,type="max_TSS")
colnames(max) = colnames(ts_index)

singular = ts_index[ts_index$type=="singular",]
max_matrix = rbind(singular,max)

a = max_matrix[max_matrix$type=="singular",]$ts_index
b = max_matrix[max_matrix$type=="max_TSS",]$ts_index
wilcox.test(a,b)$p.value
t.test(a,b)$p.value

ggplot(max_matrix, aes(x=type, y=ts_index, fill=type)) + geom_boxplot(width = 0.5) + 
  labs(title="singular-TSS VS. max TSS TSindex in TS-cluster\nTissue-Specificity index",
       subtitle = "t-test & wilcoxon siginificnat")
ggsave("hg38_mergedTSS_TSindex_maxInCluster_singular_boxplot.pdf")

ggplot(max_matrix, aes(x=ts_index, color=type)) + geom_density(size=1.2) +
  labs(title="singular-TSS VS. max TSS TSindex in TS-cluster\nTissue-Specificity index",
       subtitle = "t-test & wilcoxon siginificnat")
ggsave("hg38_mergedTSS_TSindex_maxInCluster_singular_densityLine.pdf")

#### ubi overlapped
ubi_max_matrix = max_matrix[ubi,]

a = ubi_max_matrix[ubi_max_matrix$type=="singular",]$ts_index
b = ubi_max_matrix[ubi_max_matrix$type=="max_TSS",]$ts_index
wilcox.test(a,b)$p.value
t.test(a,b)$p.value

ggplot(ubi_max_matrix, aes(x=type, y=ts_index, fill=type)) + geom_boxplot(width = 0.5) + 
  labs(title="ubi-rOCRs overlapped\nsingular-TSS VS. max TSS TSindex in TS-cluster\nTissue-Specificity index",
       subtitle = "t-test & wilcoxon siginificnat")
ggsave("hg38_ubi-rOCR_mergedTSS_TSindex_maxInCluster_singular_boxplot.pdf")

ggplot(ubi_max_matrix, aes(x=ts_index, color=type)) + geom_density(size=1.2) +
  labs(title="ubi-rOCRs overlapped\nsingular-TSS VS. max TSS TSindex in TS-cluster\nTissue-Specificity index",
       subtitle = "t-test & wilcoxon siginificnat")
ggsave("hg38_ubi-rOCR_mergedTSS_TSindex_maxInCluster_singular_densityLine.pdf")





###########