
# -- Kaili
# This script is for checking the length distribution of dhs-bins.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/")


v1 = read.table("mm10_OCR-center_bins_v1.bed", header = FALSE, row.names = 4)
v2 = read.table("mm10_OCR-center_bins_v2.bed", header = FALSE, row.names = 4)
v3 = read.table("mm10_OCR-center_bins_v3.bed", header = FALSE, row.names = 4)

###########

pdf("dhs-bins_length_histogram_v1.pdf")
hist(v1[,3]-v1[,2], freq = FALSE, main = "bins length distribution of V1\n(n=13,279,776)", col = "grey", 
     xlab = "length of bins (V1)", breaks = 99.5:300.5)
dev.off()

pdf("dhs-bins_length_histogram_v1_part.pdf")
hist(v1[,3]-v1[,2], freq = FALSE, main = "bins length distribution of V1\n(n=13,279,776)", col = "grey", 
     xlab = "length of bins (V1)", breaks = 99.5:300.5, ylim = c(0,0.005))
dev.off()

###########

pdf("dhs-bins_length_histogram_v2.pdf")
hist(v2[,3]-v2[,2], freq = FALSE, main = "bins length distribution of V2\n(n=13,426,203)", col = "grey", 
     xlab = "length of bins (V2)", breaks = 0.5:300.5)
dev.off()

pdf("dhs-bins_length_histogram_v2_part.pdf")
hist(v2[,3]-v2[,2], freq = FALSE, main = "bins length distribution of V2\n(n=13,426,203)", col = "grey", 
     xlab = "length of bins (V2)", breaks = 0.5:300.5, ylim = c(0,0.005))
dev.off()

###########

pdf("dhs-bins_length_histogram_v3.pdf")
hist(v3[,3]-v3[,2], freq = FALSE, main = "bins length distribution of V3\n(n=12,913,358)", col = "grey", 
     xlab = "length of bins (V3)", breaks = 99.5:473.5)
dev.off()

pdf("dhs-bins_length_histogram_v3_part.pdf")
hist(v3[,3]-v3[,2], freq = FALSE, main = "bins length distribution of V3\n(n=12,913,358)", col = "grey", 
     xlab = "length of bins (V3)", breaks = 99.5:473.5, ylim = c(0,0.003))
dev.off()


#########
v3_ocr = read.table("mm10-rOCRs_v3_sorted.txt")
