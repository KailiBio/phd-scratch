
# -- Kaili
# This script is for comparing merged-TSS TS-index of gene.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/")

library(ggplot2)
library(RColorBrewer)

all = read.table("GRCh38_mergedTSS_TSindex_average.txt", row.names = 1)
overlapped = read.table("GRCh38_ubi-rOCR_overlapped_mergedTSS_TSindex_average.txt", row.names = 1)
rest = read.table("GRCh38_non_ubi-rOCR_overlapped_mergedTSS_TSindex_average.txt", row.names = 1)

#################
# mergedTSS TS index: overlapped vs. all
#################
data = data.frame(cbind(overlapped, all[rownames(overlapped),]))
colnames(data) = c("overlapped", "all")

# ggplot(data, aes(x=overlapped, y=all)) + 
#   stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) + 
#   scale_fill_continuous(low = "white", high = "red") +
#   theme(legend.position = "none")
pdf("hg38_mergedTSS_TSindex_overlapped_all_gene.pdf")
smoothScatter(data[,1], data[,2], xlim=c(0.4,1), ylim=c(0.4,1), nrpoints=0, 
              xlab = "mean TS index for ubi-rOCR overlapped mergedTSSs in gene", 
              ylab = "mean TS index for all mergedTSSs in gene",
              main = "mean TS index in ubi-rOCR overlapped genes")
abline(c(0,0), c(1,1), col="red", lwd=2, lty=2)
dev.off()
#################
# mergedTSS TS index: overlapped vs. rest
#################
data2 = data.frame(cbind(overlapped, rest[rownames(overlapped),]))
colnames(data2) = c("overlapped", "non_overlapped")
data2 = na.omit(data2)


pdf("hg38_mergedTSS_TSindex_overlapped_rest_gene.pdf")
smoothScatter(data2[,1], data2[,2], xlim=c(0.4,1), ylim=c(0.4,1), nrpoints=0, 
              xlab = "mean TS index for ubi-rOCR overlapped mergedTSSs in gene", 
              ylab = "mean TS index for non-overlapped mergedTSSs in gene",
              main = "mean TS index in ubi-rOCR overlapped genes")
abline(c(0,0), c(1,1), col="red", lwd=2, lty=2)
dev.off()

#################
# mean TS index, boxplot
#################

pdf("hg38_mergedTSS_meanTSindex_in_gene_boxplot.pdf")
boxplot(data$overlapped, data2$non_overlapped, data$all, ylab = "mean tissue-specificity index",
        names = c("overlapped", "non-overlapped", "all"), col = brewer.pal(3, "Set1"),
        main = "mean of mergedTSS tissue-specificity index comparison", ylim=c(0.5,1))
dev.off()

t.test(data2[,1], data2[,2])$p.value

#################




