
# -- Kaili
# This script is for comparing TSS TSindex between TSS-cluster and singular-TSS.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/")

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

