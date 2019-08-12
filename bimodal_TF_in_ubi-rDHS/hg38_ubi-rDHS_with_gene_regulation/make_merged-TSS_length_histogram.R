
# -- Kaili
# This script is for making histogram of merged-TSS length.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/")

data = data.frame(read.table("GRCh38_merged-TSS_gene_length.txt"))
overlapped = data.frame(read.table("GRCh38_ubi-rOCR_overlapped_merged-TSS_gene_length.txt"))

summary(data[,8])
summary(overlapped[,8])

pdf("GRCh38_merged-TSS_length_hist_trimed.pdf")
hist(as.numeric(data[,8]), breaks = -0.5:272.5, ylim = c(0,1000), xlab = "length of merged-TSS",
     main = "histogram of length of merged-TSS", col = "grey", border = "grey", freq = T)
hist(overlapped[,8],breaks = -0.5:262.5, col="red", add = T)
dev.off()

pdf("GRCh38_merged-TSS_length_hist.pdf")
hist(as.numeric(data[,8]), breaks = -0.5:272.5, xlab = "length of merged-TSS",
     main = "histogram of length of merged-TSS", col = "grey", border = "grey", freq = T)
hist(overlapped[,8],breaks = -0.5:262.5, col="red", add = T)
dev.off()


#######################
# Nov 01
# new histogram
library(ggplot2)
matrix = data[,8]
matrix = transform(matrix, type="non_overlapped")
rownames(matrix) = data[,4]
colnames(matrix) = c("length", "type")
matrix$type = as.vector(matrix$type)
matrix[as.vector(overlapped[,4]),]$type = "overlapped"


ggplot(matrix, aes(length, fill = type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') + 
  xlab("length of merged-TSS") + scale_fill_manual(values = c("green", "red")) +
  labs(title="histogram of length of merged-TSS")
ggsave("GRCh38_merged-TSS_length_hist_overlapped_nonoverlapped.pdf")

