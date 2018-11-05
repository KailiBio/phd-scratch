
# -- Kaili
# This script is for analyzing tissue-specificity index of merged-TSS.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/")
library(ggplot2)

index = read.table("hg38_tissue_mergedTSS_exp_TSscore.txt", row.names = 1, header = T)
#overlapped = as.vector(read.table("GRCh38_ubi-rOCR_overlapped_merged-TSS_gene_length.txt")[,4])
gene_overlapped = read.table("hg38_tissue_ubi-rOCR_overlapped_gene_mergedTSS_TSscore.txt", row.names = 2)
overlapped = row.names(gene_overlapped[gene_overlapped$type=="overlapped",])

#################
# histogram: HK gene TS index comparison
#################

hk = as.vector(read.table("GRCh38_HK_mergedTSS_list.txt")[,1])

matrix = transform(index, type = "nonHK")
colnames(matrix) = c("index", "type")
matrix$type = as.vector(matrix$type)
matrix[hk,]$type = "housekeeping"


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

matrix = transform(matrix,three_group = "not_overlapped")
matrix$three_group = as.vector(matrix$three_group)
matrix[rownames(gene_overlapped[gene_overlapped$type == "overlapped",]),]$three_group = "TSS_overlapped"
matrix[rownames(gene_overlapped[gene_overlapped$type == "not_overlapped",]),]$three_group = "gene_overlapped_only"



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
