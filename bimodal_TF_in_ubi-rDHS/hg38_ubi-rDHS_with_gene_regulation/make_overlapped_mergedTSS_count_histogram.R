
# -- Kaili
# This script is for making histogram of numbers of merged-TSS overlapped with each ubi-rOCR/rOCR.


setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/")

ubi = read.table("GRCh38_ubi-rOCR_overlapped_merged-TSS_count.txt")
ocr = read.table("GRCh38_rOCR_overlapped_merged-TSS_count.txt")

summary(ubi[,2])
summary(ocr[,2])

pdf("ubi-rOCRs_overlapped_mergedTSS_count.pdf")
hist(ocr[,2], breaks = c(0.5:6.5), col = "grey", border = "grey", 
     xlab = "number of overlapped merged-TSS", 
     main = "number of merged-TSSs overlapped for each ubi-rOCR/rOCR")
hist(ubi[,2], breaks = c(0.5:6.5), col = "red", add=T)
dev.off()







