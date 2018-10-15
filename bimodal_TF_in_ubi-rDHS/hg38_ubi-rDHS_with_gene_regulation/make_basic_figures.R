
# -- Kaili
# This script is for making figures for basic info about ubi-rOCRs and closest genes.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/")


data = matrix(c(8272, 10959, 11254, 11494, 11888, 8962, 9814, 9988, 10124, 10348, 
                31800, 33087, 33302, 33480, 33754), ncol=3)
colnames(data) = c("ubi.rOCRs", "gene", "tss")
data = transform(distance=c("overlapped", "<2kb", "<5kb", "<10kb", "closest"), data)


########################
# basic num
########################
library(ggplot2)
library(gridExtra)

pdf("closest_gene_number.pdf", width=12)
a = ggplot(data, aes(x=distance, y=ubi.rOCRs, label = ubi.rOCRs)) + 
  geom_bar(stat="identity", fill = "grey") + 
  scale_x_discrete(limits = c("overlapped", "<2kb", "<5kb", "<10kb", "closest")) + theme_minimal() +
  ylab("number of ubi-rOCRs") + geom_text(position="stack")
b = ggplot(data, aes(x=distance, y=gene, label = gene)) + geom_bar(stat="identity", fill = "grey") + 
  scale_x_discrete(limits = c("overlapped", "<2kb", "<5kb", "<10kb", "closest")) + theme_minimal() +
  ylab("number of genes") + geom_text(position="stack")
c = ggplot(data, aes(x=distance, y=tss, label = tss)) + geom_bar(stat="identity", fill = "grey") + 
  scale_x_discrete(limits = c("overlapped", "<2kb", "<5kb", "<10kb", "closest")) + theme_minimal() +
  ylab("number of TSSs") + geom_text(position="stack")
grid.arrange(a, b, c, nrow=1)
dev.off()


########################
# ubi-rOCR overlapped gene/TSS
########################

# tss_count = read.table("GRCh38_ubi-rOCR_overlapped_TSS_count.txt")
# 
# pdf("ubi-rOCRs_overlapped_TSS_count.pdf", width=9)
# hist(tss_count[,1], breaks = 0.5:26.5, xlab = "number of overlapped TSSs",
#      main = "number of TSSs overlapped for each ubi-rOCR", col = "grey")
# dev.off()

tss_count_all = read.table("GRCh38_rOCR_overlapped_TSS_count.txt", row.names = 2)
ubi.rOCR = as.vector(read.table("GRCh38_ubi-rOCR_overlapped_TSS_count.txt")[,2])

matrix = transform(tss_count_all, type = "rOCR")
colnames(matrix) = c("num", "type")
matrix$type = as.vector(matrix$type)
matrix[ubi.rOCR,]$type = "ubi-rOCR"

library(ggplot2)

ggplot(matrix, aes(num, fill = type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity', binwidth = 1) + 
  xlab("number of overlapped TSSs") + scale_fill_manual(values = c("#67a9cf", "#ef8a62")) + 
  labs(title="number of TSSs overlapped for each rOCR/ubi-rOCR")
ggsave("ubi-rOCRs_overlapped_TSS_count.pdf")




########################