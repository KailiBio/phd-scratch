
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

tss_count = read.table("GRCh38_ubi-rOCR_overlapped_TSS_count.txt")

pdf("ubi-rOCRs_overlapped_TSS_count.pdf", width=9)
hist(tss_count[,1], breaks = 0.5:26.5, xlab = "number of overlapped TSSs",
     main = "number of TSSs overlapped for each ubi-rOCR", col = "grey")
dev.off()



########################