
# -- Kaili
# This script is for making histogram of merged-TSS length.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/")

data = data.frame(read.table("GRCh38_merged-TSS_gene_length.txt"))
overlapped = data.frame(read.table("GRCh38_ubi-rOCR_overlapped_merged-TSS_gene_length.txt"))

summary(data[,8])
summary(overlapped[,8])

pdf("GRCh38_merged-TSS_length_hist.pdf")
hist(data[,8], breaks = -0.5:272.5, ylim = c(0,1000), xlab = "length of merged-TSS",
     main = "histogram of length of merged-TSS", col = "grey", border = "grey")
hist(overlapped[,8],breaks = -0.5:262.5, col="red", add = T)
dev.off()

