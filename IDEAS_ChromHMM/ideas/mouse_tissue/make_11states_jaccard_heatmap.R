
# -- Kaili
# This script is for making Jaccard heatmap for 11 states in 11 samples.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/state_CTCF_signal/")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

#####################

matrix = read.table("ctcf_states_overlapped.txt", header = TRUE)
rownames(matrix) = colnames(matrix)
matrix_union = read.table("ctcf_states_union.txt", header = TRUE)
rownames(matrix_union) = colnames(matrix_union)


library(pheatmap)
pdf("11_states_between_all_samples.pdf")
pheatmap(matrix/matrix_union, display_numbers = TRUE, 
         main = "11 states in each cell-types",
         breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
         col = brewer.pal(10,"RdYlBu")[10:1])
dev.off()
