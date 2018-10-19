
# -- Kaili
# This script for making scatter plot for tissue-specificity index of TSS in genes.


setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/tsindex/")

library(ggplot2)
library(RColorBrewer)

all = read.table("hg38_gene_TSS_TSindex_average.txt", row.names = 1)
overlapped = read.table("hg38_ubi-rOCR_overlapped_gene_TSS_TSindex_average.txt", row.names = 1)
rest = read.table("hg38_non_overlapped_gene_TSS_TSindex_average.txt", row.names = 1)

gene = read.table("hg38_tissue_gene_exp_TSscore.txt", row.names = 1, header = TRUE)
all_max = read.table("hg38_gene_TSS_TSindex_max.txt", row.names = 1, header = TRUE)

#################
# TSS TS index: overlapped vs. all
#################
data = data.frame(cbind(overlapped, all[rownames(overlapped),]))
colnames(data) = c("overlapped", "all")

# ggplot(data, aes(x=overlapped, y=all)) + 
#   stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) + 
#   scale_fill_continuous(low = "white", high = "red") +
#   theme(legend.position = "none")
pdf("hg38_TSS_TSindex_overlapped_all_gene.pdf")
smoothScatter(data[,1], data[,2], xlim=c(0.4,1), ylim=c(0.4,1), nrpoints=0, 
              xlab = "mean TS index for ubi-rOCR overlapped TSSs in gene", 
              ylab = "mean TS index for all TSSs in gene",
              main = "mean TS index in ubi-rOCR overlapped genes")
abline(c(0,0), c(1,1), col="red", lwd=2, lty=2)
dev.off()
#################
# TSS TS index: overlapped vs. rest
#################
data2 = data.frame(cbind(overlapped, rest[rownames(overlapped),]))
colnames(data2) = c("overlapped", "non_overlapped")
data2 = na.omit(data2)
data2 = data2[- which(rownames(data2)=="ENSG00000133980.4"),]
data2 = data2[- which(rownames(data2)=="ENSG00000259867.5"),]

data2 = data2["ENSG00000133980.4",]

pdf("hg38_TSS_TSindex_overlapped_rest_gene.pdf")
smoothScatter(data2[,1], data2[,2], xlim=c(0.4,1), ylim=c(0.4,1), nrpoints=0, 
              xlab = "mean TS index for ubi-rOCR overlapped TSSs in gene", 
              ylab = "mean TS index for non-overlapped TSSs in gene",
              main = "mean TS index in ubi-rOCR overlapped genes")
abline(c(0,0), c(1,1), col="red", lwd=2, lty=2)
dev.off()


#################
# mean TS index, boxplot
#################

pdf("hg38_TSS_meanTSindex_in_gene_boxplot.pdf")
boxplot(data$overlapped, data2$non_overlapped, data$all, ylab = "mean tissue-specificity index",
        names = c("overlapped", "non-overlapped", "all"), col = brewer.pal(3, "Set1"),
        main = "mean TSS tissue-specificity index comparison")
dev.off()

t.test(data2[,1], data2[,2])$p.value


#################
# TS index: overlapped TSS vs. gene
#################
data3 = data.frame(cbind(overlapped, gene[rownames(overlapped),]))
colnames(data3) = c("overlapped", "gene")

pdf("hg38_TSindex_overlapped_gene.pdf")
smoothScatter(data3[,1], data3[,2], xlim=c(0.4,1), ylim=c(0.4,1), nrpoints=0, 
              xlab = "mean TS index for ubi-rOCR overlapped TSSs in gene", 
              ylab = "TS index of gene (RNA-seq)",
              main = "TS index in ubi-rOCR overlapped genes")
abline(c(0,0), c(1,1), col="red", lwd=2, lty=2)
dev.off()

#################
# TS index: mean TSS vs. gene
#################
data4 = data.frame(cbind(all[rownames(overlapped),], gene[rownames(overlapped),]))
colnames(data4) = c("all", "gene")

pdf("hg38_TSindex_meanTSS_gene.pdf")
smoothScatter(data4[,1], data4[,2], xlim=c(0.4,1), ylim=c(0.4,1), nrpoints=0, 
              xlab = "mean TS index for all TSSs in gene (RAMPAGE)", 
              ylab = "TS index of gene (RNA-seq)",
              main = "TS index in ubi-rOCR overlapped genes")
abline(c(0,0), c(1,1), col="red", lwd=2, lty=2)
dev.off()


#################
# TS index: max TSS vs. gene
#################

pdf("hg38_TSindex_maxTSS_gene.pdf")
smoothScatter(all_max[rownames(overlapped),], gene[rownames(overlapped),], xlim=c(0.4,1), ylim=c(0.4,1), 
              nrpoints=0, 
              xlab = "max TS index for all TSSs in gene (RAMPAGE)", 
              ylab = "TS index of gene (RNA-seq)",
              main = "TS index in ubi-rOCR overlapped genes")
abline(c(0,0), c(1,1), col="red", lwd=2, lty=2)
dev.off()

#################
# boxplot: HK gene TS index comparison
#################
hk = as.vector(read.table("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/hg38_housekeeping_geneID_geneType_geneSymbol.txt")[,1])

pdf("gene_hk_TSindex_comparison.pdf", width = 10, height = 10)
boxplot(gene[hk,1], all[hk,1], all_max[hk,1], gene[,1], all[,1], all_max[,1],
        col = c(brewer.pal(3, "Set1"),brewer.pal(3, "Set1")), ylim = c(0.4,1), main = "tissue-specificity index comparison",
        names = c("housekeeping\ngene", "housekeeping\nmean TSS", "housekeeping\nmax TSS", "gene", "mean TSS", "max TSS"),
        ylab = "tissue-specificity index")
dev.off()

#################

