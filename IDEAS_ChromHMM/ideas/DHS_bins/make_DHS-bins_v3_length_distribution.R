

# -- Kaili
# This script is for making length distribution of DHS-center bins.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/")

dhs_bins = read.table("mm10-rOCRs_v3_sorted_3.txt")
gap_bins = read.table("mm10_OCR-center_bins_v3_gap.bed")
bins = rbind(dhs_bins[,1:4], gap_bins)

####
pdf("DHS-bins_v3_length_distrobution.pdf")

summary(dhs_bins[,10])
hist(dhs_bins[,10], breaks = 149.5:473.5, xlab = "length of DHSs", 
     main = "length distribution of DHSs in DHS-bins", col = "grey", border = "grey")

hist(dhs_bins[,10], breaks = 149.5:473.5, xlab = "length of DHSs", 
     main = "length distribution of DHSs in DHS-bins", col = "grey", border = "grey", freq = F)

summary(gap_bins[,3]-gap_bins[,2])
hist(gap_bins[,3]-gap_bins[,2], breaks = 99.5:350.5, xlab = "length of DHS_Gaps", 
     main = "length distribution of DHS_Gaps in DHS-bins", col = "grey", border = "grey")

hist(gap_bins[,3]-gap_bins[,2], breaks = 99.5:350.5, xlab = "length of DHS_Gaps", 
     main = "length distribution of DHS_Gaps in DHS-bins", col = "grey", border = "grey", freq = F)

summary(bins[,3]-bins[,2])
hist(bins[,3]-bins[,2], breaks = 99.5:473.5, xlab = "length of DHS-bins", 
     main = "length distribution of DHS-bins", col = "grey", border = "grey")

hist(bins[,3]-bins[,2], breaks = 99.5:473.5, xlab = "length of DHS-bins", 
     main = "length distribution of DHS-bins", col = "grey", border = "grey", freq = F)

dev.off()


###############
# Feb20
colnames(bins) = c("chr", "s", "e", "id")
bins = transform(bins, length = bins$e-bins$s)

ggplot(bins, aes(x = length)) + 
  geom_histogram(binwidth = 1, aes(y=..density..), fill = "#AEAFAE") +
  theme_minimal() + theme(title = element_text(face = "bold", size=12)) +
  xlab("length of DHS-center bins")

ggsave("density_histogram_of_bin_length.pdf")

###############
