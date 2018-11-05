
# -- Kaili
# This script is for comparing the TSindex calculate by only tissues or all samples.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/")

gene_TSindex_all = read.table("hg38_all_gene_exp_TSscore.txt", header = T, row.names = 1)
gene_TSindex_tissue = read.table("hg38_tissue_gene_exp_TSscore.txt", header = T, row.names = 1)
ubi_gene = as.vector(read.table("GRCh38_ubi-rOCR_overlapped_gene_id.txt")[,1])

TSS_TSindex_all = read.table("hg38_all_TSS_exp_TSscore.txt", header = T, row.names = 1)
TSS_TSindex_tissue = read.table("hg38_tissue_TSS_exp_TSscore_uniqID.txt", row.names = 1)
ubi_TSS = as.vector(read.table("GRCh38_ubi-rOCR_overlapped_TSS_uniqID.txt")[,1])


#############
# for gene
pdf("TSindex_comparison_all_gene.pdf")
smoothScatter(gene_TSindex_all$all, gene_TSindex_tissue$all, xlab = "all biosamples (n=184)", 
              ylab = "only tissue samples (n=104)", main = "TS index for all genes", 
              xlim=c(0.4,1), ylim=c(0.4,1))
lines(c(0,1), c(0,1), col="red")
dev.off()

#
pdf("TSindex_comparison_ubi-rOCR_overlapped_gene.pdf")
smoothScatter(gene_TSindex_all[ubi_gene,], gene_TSindex_tissue[ubi_gene,], 
              xlab = "all biosamples (n=184)", ylab = "only tissue samples (n=104)", 
              main = "TS index for ubi-rOCR overlapped genes", xlim=c(0.4,1), ylim=c(0.4,1))
lines(c(0,1), c(0,1), col="red")
dev.off()


#############
# for TSS
pdf("TSindex_comparison_all_TSS.pdf")
smoothScatter(TSS_TSindex_all$all, TSS_TSindex_tissue$V2, xlab = "all biosamples (n=155)", 
              ylab = "only tissue samples (n=104)", main = "TS index for all genes", 
              xlim=c(0.7,1), ylim=c(0.7,1))
lines(c(0,1), c(0,1), col="red")
dev.off()

#
pdf("TSindex_comparison_ubi-rOCR_overlapped_TSS.pdf")
smoothScatter(TSS_TSindex_all[ubi_TSS,], TSS_TSindex_tissue[ubi_TSS,], 
              xlab = "all biosamples (n=155)", ylab = "only tissue samples (n=104)", 
              main = "TS index for ubi-rOCR overlapped TSSs", xlim=c(0.65,1), ylim=c(0.65,1))
lines(c(0,1), c(0,1), col="red")
dev.off()





