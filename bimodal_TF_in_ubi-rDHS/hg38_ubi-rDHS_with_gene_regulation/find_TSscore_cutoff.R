
# -- Kaili
# This script is for plotting tissue specificity index and finding the cut-off.

# Update at Nov21
# make density line plot

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/")

calculate_TSindex = function(x){
  s=0
  if(max(x)==0){return(NA)}
  else{
    for(i in x){
      s=s+((max(x)-i)/max(x))
    }
    return(s/(length(x)-1))
  }
}

######################
# Xiaoou's data
######################
data = read.table("all_ts.txt", header = TRUE, row.names = 1)
pdf("mouse_ts_score_xiaoou.pdf")
hist(data[,1], breaks=200, xlab = "tissue specificity index", 
     main = "Xiaoou's mouse data\n(66 samples)", freq = FALSE)
lines(c(0.52,0.52), c(0,400), col="red", lwd=2, lty=2)
dev.off()

######################
# 
######################
ts_score = read.table("hg38_tissue_gene_exp_TSscore.txt", header = TRUE, row.names = 1)
exp = read.table("hg38_tissue_gene_exp_matrix.txt", header = TRUE, row.names = 1)
hk = as.vector(read.table("hg38_housekeeping_geneID_geneType_geneSymbol.txt")[,1])


matrix = transform(ts_score, type="rest genes")
matrix$type = as.vector(matrix$type)
matrix[hk,]$type = "housekeeping genes"

library(ggplot2)

hk_v = matrix[matrix$type == "housekeeping genes", ]$all
nonhk_v = matrix[matrix$type == "rest genes", ]$all
t.test(hk_v, nonhk_v)$p.value
# p-value < 2.2e-16
wilcox.test(hk_v, nonhk_v)$p.value

ggplot(matrix, aes(all, fill = type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') +
  xlab("tissue_specificity_index") + scale_fill_manual(values =c("#af8dc3", "#7fbf7b") ) +
  labs(title="gene Tissue-Specificity index (RNA-seq)", 
       subtitle="t-test: p-value < 2.2e-16\nwilcoxon: p-value < 2.2e-16")
ggsave("hg38_gene_TSindex_HK.pdf")

# density line
ggplot(matrix, aes(all, color=type)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") + 
  scale_color_manual(values = c("#984ea3", "#4daf4a")) +
  labs(title="Tissue-Specificity index of gene", 
       subtitle="t-test: p-value < 2.2e-16\nwilcoxon: p-value < 2.2e-16")
ggsave("hg38_gene_TSindex_HK_densityLine.pdf")

##########
overlapped = as.vector(read.table("GRCh38_ubi-rOCR_closest_gene_list.bed")[,1])
matrix = transform(matrix, ubi.rOCR = "rest genes")
matrix$ubi.rOCR = as.vector(matrix$ubi.rOCR)
matrix[overlapped,]$ubi.rOCR = "overlapped genes"

overlapped_v = matrix[matrix$ubi.rOCR == "overlapped genes", ]$all
non_overlapped_v = matrix[matrix$ubi.rOCR == "rest genes", ]$all
t.test(overlapped_v, non_overlapped_v)$p.value
# p-value < 2.2e-16
wilcox.test(overlapped_v, non_overlapped_v)$p.value

library(ggplot2)

ggplot(matrix, aes(all, fill = ubi.rOCR)) + 
  geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') +
  xlab("tissue_specificity_index") + scale_fill_manual(values = c("#67a9cf","#ef8a62")) +
  labs(title="gene Tissue-Specificity index (RNA-seq)",
       subtitle="t-test: p-value < 2.2e-16\nwilcoxon: p-value < 2.2e-16")
ggsave("hg38_gene_TSindex_overlapped.pdf")

# density line
ggplot(matrix, aes(all, color=ubi.rOCR)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") + 
  scale_color_manual(values = c("#e41a1c","#4daf4a")) +
  labs(title="Tissue-Specificity index of genes", 
       subtitle="t-test: p-value < 2.2e-16\nwilcoxon: p-value < 2.2e-16")
ggsave("hg38_gene_TSindex_overlapped_densityLine.pdf")


######################

# library(preprocessCore)
# exp_q = normalize.quantiles(as.matrix(exp), copy = TRUE)
# ts_score_q=apply(exp_q,1,fun_tes)
# hist(ts_score_q, breaks=100)


######################
# TSS ts index
######################
tss_ts = read.table("hg38_tissue_TSS_exp_TSscore.txt", header = TRUE, row.names = 1)

overlapped_tss = as.vector(read.table("GRCh38_ubi-rOCR_closest_TSS_list.bed")[,1])
matrix2 = transform(tss_ts, ubi.rOCR="non-overlapped")
matrix2$ubi.rOCR = as.vector(matrix2$ubi.rOCR)
matrix2[overlapped_tss,]$ubi.rOCR = "overlapped"

library(ggplot2)


ggplot(matrix2, aes(all, fill = ubi.rOCR)) + 
  geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + 
  xlab("tissue_specificity_index") + scale_fill_manual(values = c("#af8dc3", "#7fbf7b")) +
  labs(title="Tissue-specificity index for TSSs")
ggsave("TSS_tissue_specificity_index.pdf")

#######################
# Dec 09
## try to normalize gene
#######################
library(preprocessCore)
a = normalize.quantiles(as.matrix(exp))
colnames(a) = colnames(exp)
rownames(a) = rownames(exp)
write.table(a,file = "hg38_tissue_gene_exp_matrix_quantile.txt", sep="\t", quote = FALSE)

ts_quantile = read.table("hg38_tissue_gene_exp_TSscore_quantile.txt", header = TRUE, row.names = 1)
hk = as.vector(read.table("hg38_housekeeping_geneID_geneType_geneSymbol.txt")[,1])
matrix = transform(ts_quantile, type="rest genes")
matrix$type = as.vector(matrix$type)
matrix[hk,]$type = "housekeeping genes"
# density line
ggplot(matrix, aes(all, color=type)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") + 
  scale_color_manual(values = c("#984ea3", "#4daf4a")) +
  labs(title="Tissue-Specificity index of gene", 
       subtitle="quantile normalized signal")
ggsave("hg38_gene_TSindex_HK_densityLine_quantile.pdf")
ggsave("hg38_gene_TSindex_HK_densityLine_quantile.png")
####
overlapped = as.vector(read.table("GRCh38_ubi-rOCR_closest_gene_list.bed")[,1])
matrix = transform(matrix, ubi.rOCR = "rest genes")
matrix$ubi.rOCR = as.vector(matrix$ubi.rOCR)
matrix[overlapped,]$ubi.rOCR = "overlapped genes"
ggplot(matrix, aes(all, color=ubi.rOCR)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") + 
  scale_color_manual(values = c("#e41a1c","#4daf4a")) +
  labs(title="Tissue-Specificity index of genes", 
       subtitle="quantile normalized signal")
ggsave("hg38_gene_TSindex_overlapped_densityLine_quantile.pdf")
ggsave("hg38_gene_TSindex_overlapped_densityLine_quantile.png")

####
summary(a[,1])
hist(log10(a[,1]+0.1))

### max 100
tx_max100 = read.table("hg38_tissue_gene_exp_TSscore_quantile_max100.txt", header = TRUE, row.names = 1)
hk = as.vector(read.table("hg38_housekeeping_geneID_geneType_geneSymbol.txt")[,1])
matrix = transform(tx_max100, type="rest genes")
matrix$type = as.vector(matrix$type)
matrix[hk,]$type = "housekeeping genes"
# density line
ggplot(matrix, aes(all, color=type)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") + 
  scale_color_manual(values = c("#984ea3", "#4daf4a")) +
  labs(title="Tissue-Specificity index of gene", 
       subtitle="quantile & max signal=100")
ggsave("hg38_gene_TSindex_HK_densityLine_quantile_max100.pdf")
ggsave("hg38_gene_TSindex_HK_densityLine_quantile_max100.png")
####
overlapped = as.vector(read.table("GRCh38_ubi-rOCR_closest_gene_list.bed")[,1])
matrix = transform(matrix, ubi.rOCR = "rest genes")
matrix$ubi.rOCR = as.vector(matrix$ubi.rOCR)
matrix[overlapped,]$ubi.rOCR = "overlapped genes"
ggplot(matrix, aes(all, color=ubi.rOCR)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") + 
  scale_color_manual(values = c("#e41a1c","#4daf4a")) +
  labs(title="Tissue-Specificity index of genes", 
       subtitle="quantile & max signal=100")
ggsave("hg38_gene_TSindex_overlapped_densityLine_quantile_max100.pdf")
ggsave("hg38_gene_TSindex_overlapped_densityLine_quantile_max100.png")

### max 20
tx_max100 = read.table("hg38_tissue_gene_exp_TSscore_quantile_max20.txt", header = TRUE, row.names = 1)
hk = as.vector(read.table("hg38_housekeeping_geneID_geneType_geneSymbol.txt")[,1])
matrix = transform(tx_max20, type="rest genes")
matrix$type = as.vector(matrix$type)
matrix[hk,]$type = "housekeeping genes"
# density line
ggplot(matrix, aes(all, color=type)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") + 
  scale_color_manual(values = c("#984ea3", "#4daf4a")) +
  labs(title="Tissue-Specificity index of gene", 
       subtitle="quantile & max signal=20")
ggsave("hg38_gene_TSindex_HK_densityLine_quantile_max20.pdf")
ggsave("hg38_gene_TSindex_HK_densityLine_quantile_max20.png")
####
overlapped = as.vector(read.table("GRCh38_ubi-rOCR_closest_gene_list.bed")[,1])
matrix = transform(matrix, ubi.rOCR = "rest genes")
matrix$ubi.rOCR = as.vector(matrix$ubi.rOCR)
matrix[overlapped,]$ubi.rOCR = "overlapped genes"
ggplot(matrix, aes(all, color=ubi.rOCR)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") + 
  scale_color_manual(values = c("#e41a1c","#4daf4a")) +
  labs(title="Tissue-Specificity index of genes", 
       subtitle="quantile & max signal=20")
ggsave("hg38_gene_TSindex_overlapped_densityLine_quantile_max20.pdf")
ggsave("hg38_gene_TSindex_overlapped_densityLine_quantile_max20.png")
#######################
