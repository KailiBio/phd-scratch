
# -- Kaili
# This script is for making state combination histogram.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/")

# state_combine = read.table("OCR-center_bins_nearest_pairs_state_combine.txt")
# colnames(state_combine) = c("close_OCR", "normal")
# 
# close_OCR_state = data.frame(unclass(table(state_combine$close_OCR)))
# #summary(close_OCR_state[,1])
# #hist(close_OCR_state[,1], breaks=0.5:24151.5, ylim=c(0,10))
# barplot(close_OCR_state[,1])
# 
# close_normal_state = data.frame(unclass(table(state_combine$normal)))
# summary(close_normal_state[,1])
# barplot(close_normal_state[,1])
# 

#########
state_combine = read.table("OCR-center_bins_nearest_pairs_state.txt")
close_ocr_matrix = matrix(0, ncol = 38, nrow=38)
colnames(close_ocr_matrix) = rownames(close_ocr_matrix) = c(0:37)
close_normal_matrix = matrix(0, ncol = 38, nrow=38)
colnames(close_normal_matrix) = rownames(close_normal_matrix) = c(0:37)

for(i in 1:nrow(state_combine)){
  ocr1 = state_combine[i,2]
  ocr2 = state_combine[i,4]
  normal1 = state_combine[i,6]
  normal2 = state_combine[i,8]
  close_ocr_matrix[ocr1,ocr2] = close_ocr_matrix[ocr1,ocr2] + 1
  close_normal_matrix[normal1,normal2] = close_normal_matrix[normal1,normal2] + 1
}


library(pheatmap)
library(RColorBrewer)
pdf("nearby_bins_state_transition_heatmap.pdf")
pheatmap(close_ocr_matrix, cluster_rows=FALSE, cluster_cols=FALSE, 
         breaks = c(0,10,50,100,500,1000,2000,3000,4000,5000),
         col = brewer.pal(11, "Greys"), main = "state in nearby OCR bins")
#
pheatmap(close_normal_matrix, cluster_rows=FALSE, cluster_cols=FALSE, 
         breaks = c(0,10,50,100,500,1000,2000,3000,4000,5000),
         col = brewer.pal(11, "Greys"), main = "state in matched normal bins")
dev.off()