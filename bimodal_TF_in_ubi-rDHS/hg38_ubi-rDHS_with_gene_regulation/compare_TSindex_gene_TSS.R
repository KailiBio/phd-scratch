
# -- Kaili
# This script is for comparing tissue-specificity index between gene and TSS.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/")

data = read.table("hg38_TSindex_TSS_gene.txt", header = TRUE)
ubi = read.table("GRCh38_ubi-rOCRs.bed")

dat = data[ubi[,4],]


pdf("TSindex_all_gene_TSS.pdf", width = 6, height = 6)
smoothScatter(data$ts1, data$ts2, xlim = c(0.4,1), ylim = c(0.4,1), 
              xlab = "tissue-specificity index of all TSSs (n=194,897)", 
              ylab = "tissue-specificity index of all genes",
              main = "tissue-specificity index between all TSSs and genes")
dev.off()

pdf("TSindex_ubi_gene_TSS.pdf", width = 6, height = 6)
smoothScatter(dat$ts1, dat$ts2, xlim = c(0.4,1), ylim = c(0.4,1), 
              xlab = "tissue-specificity index of ubi-rOCR overlapped TSSs (n=11,897)", 
              ylab = "tissue-specificity index of ubi-rOCR overlapped genes",
              main = "tissue-specificity index between\nubi-rOCR overlapped TSSs and genes")
dev.off()








