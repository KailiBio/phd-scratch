
# -- Kaili
# This script for comparing gene TSindex(RNAseq) with RAMPAGE signal weighted TSS TSindex.


setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/")

weighted = read.table("hg38_mergedTSS_gene_weigedt_TSindex.txt", row.names = 1)
gene_ts = read.table("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/hg38_tissue_gene_exp_TSscore.txt", header = TRUE, row.names = 1)

matrix = transform(weighted, normal = gene_ts[rownames(weighted),])
colnames(matrix) = c("weighted", "normal")

pdf("gene_TSindex_weightedTSindex_scatterPlot.pdf")
smoothScatter(matrix$weighted, matrix$normal, xlim = c(0,1), ylim=c(0,1), main = "gene TS index comparison",
              xlab = "sum of weighted TSS TS-index (RAMPAGE)", ylab = "gene TS-index (RNA-seq)")
lines(c(0,1), c(0,1), col="red", lwd=1.5)
dev.off()

# matrix[matrix$normal>0.81 & matrix$weighted<0.2,]

