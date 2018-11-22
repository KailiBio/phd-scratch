
# -- Kaili
# This script is for analyzing tissue-specificity index of merged-TSS.

# Update on Nov21
# add statistic test
# get density line figures

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/")
library(ggplot2)

index = read.table("hg38_tissue_mergedTSS_exp_TSscore.txt", row.names = 1, header = T)
#overlapped = as.vector(read.table("GRCh38_ubi-rOCR_overlapped_merged-TSS_gene_length.txt")[,4])
gene_overlapped = read.table("hg38_tissue_ubi-rOCR_overlapped_gene_mergedTSS_TSscore.txt", row.names = 2)
overlapped = row.names(gene_overlapped[gene_overlapped[,3]=="overlapped",])

#################
# histogram: HK gene TS index comparison
#################

hk = as.vector(read.table("GRCh38_HK_mergedTSS_list.txt")[,1])

matrix = transform(index, type = "rest TSSs")
colnames(matrix) = c("index", "type")
matrix$type = as.vector(matrix$type)
matrix[hk,]$type = "housekeeping TSSs"

hk_v = matrix[matrix$type == "housekeeping TSSs", ]$index
nonhk_v = matrix[matrix$type == "rest TSSs", ]$index
t.test(hk_v, nonhk_v)$p.value
# p-value < 2.2e-16
wilcox.test(hk_v, nonhk_v)$p.value

ggplot(matrix, aes(index, fill=type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') +
  xlab("tissue_specificity_index") + scale_fill_manual(values = c("#af8dc3", "#7fbf7b")) +
  labs(title="Tissue-specificity index for merged-TSS")
ggsave("mergedTSS_gene_hk_TSindex_comparison.pdf")

ggplot(matrix, aes(index, fill=type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') +
  xlab("tissue_specificity_index") + scale_fill_manual(values = c("#af8dc3", "#7fbf7b")) +
  labs(title="Tissue-specificity index for merged-TSS")
ggsave("mergedTSS_gene_hk_TSindex_comparison_density.pdf")


# density line
ggplot(matrix, aes(index, color=type)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") + 
  scale_color_manual(values = c("#984ea3", "#4daf4a")) +
  labs(title="Tissue-Specificity index of merged-TSS",
       subtitle="t-test: p-value < 2.2e-16\nwilcoxon: p-value < 2.2e-16")
ggsave("mergedTSS_gene_hk_TSindex_comparison_densityLine.pdf")

#################
# histogram: ubi-rOCR overlapped vs. non-overlapped index comparison
#################
matrix = transform(index, type = "non-overlapped")
colnames(matrix) = c("index", "type")
matrix$type = as.vector(matrix$type)
matrix[overlapped,]$type = "overlapped"

ggplot(matrix, aes(index, fill=type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') +
  xlab("tissue_specificity_index") + scale_fill_manual(values = c("#67a9cf", "#ef8a62")) +
  labs(title="Tissue-specificity index for merged-TSS")
ggsave("mergedTSS_gene_ubi-rOCR_overlapped_TSindex_comparison.pdf")

ggplot(matrix, aes(index, fill=type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') +
  xlab("tissue_specificity_index") + scale_fill_manual(values = c("#67a9cf","#ef8a62")) +
  labs(title="Tissue-specificity index for merged-TSS")
ggsave("mergedTSS_gene_ubi-rOCR_overlapped_TSindex_comparison_density.pdf")

#################
# histogram: ubi-rOCR overlapped vs. non-overlapped tss in overlapped genes index comparison
#################
colnames(gene_overlapped) = c("gene","index", "type")

matrix = transform(matrix,three_group = "rest TSSs")
matrix$three_group = as.vector(matrix$three_group)
matrix[rownames(gene_overlapped[gene_overlapped$type == "overlapped",]),]$three_group = "overlapped TSSs"
matrix[rownames(gene_overlapped[gene_overlapped$type == "not_overlapped",]),]$three_group = "TSSs that in overlapped genes"

overlapped_v = matrix[matrix$three_group == "overlapped TSSs", ]$index
non_overlapped_v = matrix[matrix$three_group == "TSSs that in overlapped genes", ]$index
t.test(overlapped_v, non_overlapped_v)$p.value
# p-value < 2.2e-16
wilcox.test(overlapped_v, non_overlapped_v)$p.value


ggplot(matrix, aes(index, fill=three_group)) + 
  geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') +
  xlab("tissue_specificity_index") + scale_fill_manual(values = c("#1a9850","#67a9cf","#ef8a62")) +
  labs(title="Tissue-specificity index for merged-TSS")
ggsave("ubi-rOCR_overlapped_gene_mergedTSS_TSindex_comparison.pdf")

ggplot(matrix, aes(index, fill=three_group)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') +
  xlab("tissue_specificity_index") + scale_fill_manual(values = c("#1a9850","#67a9cf","#ef8a62")) +
  labs(title="Tissue-specificity index for merged-TSS")
ggsave("ubi-rOCR_overlapped_gene_mergedTSS_TSindex_comparison_density.pdf")



# density line
ggplot(matrix, aes(index, color=three_group)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") + 
  scale_color_manual(values = c("#e41a1c","#4daf4a","#377eb8")) +
  labs(title="Tissue-Specificity index of merged-TSS",
       subtitle="t-test: p-value < 2.2e-16\nwilcoxon: p-value < 2.2e-16")
ggsave("ubi-rOCR_overlapped_gene_mergedTSS_TSindex_comparison_densityLine.pdf")

#################
# statictic test
#################
wilcox.test(na.omit(matrix[matrix$three_group=="TSS_overlapped",]$index), 
            na.omit(matrix[matrix$three_group=="gene_overlapped_only",]$index))$p.value

t.test(na.omit(matrix[matrix$three_group=="TSS_overlapped",]$index), 
       na.omit(matrix[matrix$three_group=="gene_overlapped_only",]$index))$p.value

#################
