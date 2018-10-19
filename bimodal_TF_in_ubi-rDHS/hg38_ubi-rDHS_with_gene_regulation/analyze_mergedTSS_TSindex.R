
# -- Kaili
# This script is for analyzing tissue-specificity index of merged-TSS.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/")
library(ggplot2)

#################
# histogram: HK gene TS index comparison
#################
index = read.table("hg38_tissue_mergedTSS_exp_TSscore.txt", row.names = 1, header = T)
hk = as.vector(read.table("GRCh38_HK_mergedTSS_list.txt")[,1])

matrix = transform(index, type = "nonHK")
colnames(matrix) = c("index", "type")
matrix$type = as.vector(matrix$type)
matrix[hk,]$type = "housekeeping"


pdf("mergedTSS_gene_hk_TSindex_comparison.pdf", width = 10, height = 10)

dev.off()








