
# -- Kaili
# This script is for plotting tissue specificity index and finding the cut-off.

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

matrix = transform(ts_score, type="nonHK")
matrix$type = as.vector(matrix$type)
matrix[hk,]$type = "HK_published"

library(ggplot2)

ggplot(matrix, aes(all, fill = type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') +
  xlab("tissue_specificity_index") + scale_fill_manual(values = c("#ef8a62", "#67a9cf")) +
  labs(title="Tissue-specificity index for genes")
ggsave("tissue_specificity_index.pdf")


overlapped = as.vector(read.table("GRCh38_ubi-rOCR_closest_gene_list.bed")[,1])
matrix = transform(matrix, ubi.rOCR = "non-overlapped")
matrix$ubi.rOCR = as.vector(matrix$ubi.rOCR)
matrix[overlapped,]$ubi.rOCR = "overlapped"

library(ggplot2)

ggplot(matrix, aes(all, fill = ubi.rOCR)) + 
  geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') +
  xlab("tissue_specificity_index") + scale_fill_manual(values = c("#af8dc3", "#7fbf7b")) +
  labs(title="Tissue-specificity index for genes")
ggsave("tissue_specificity_index2.pdf")


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


