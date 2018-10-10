
# -- Kaili
# This script is for making histogram of mm10 DHSs.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dns_bins/")

####
data = read.table("mm10-rOCRs.bed")
dat = data$V3-data$V2

summary(dat)

pdf("histogram_length_mm10_rOCRs.pdf")
hist(dat, breaks = 149.5:350.5, xlab = "length of rOCRs", main = "histogram of length of mm10 rCORs")
dev.off()

#####

data2 = read.table("mm10_rOCR_gap.bed")
dat2 = data2$V3-data2$V2
summary(dat2)
summary(dat2[dat2<1000])
a=dat2[dat2<1000]
b=dat2[dat2<100]
pdf("histogram_length_mm10_rOCRs_gap_1000.pdf")
hist(a, breaks = 0.5:999.5, xlab = "length of rOCR gaps", 
     main = "histogram of length of GRCh38 rCOR gaps")
dev.off()
#
pdf("histogram_length_mm10_rOCRs_gap_100.pdf")
hist(b, breaks = 0.5:99.5, xlab = "length of rOCR gaps", 
     main = "histogram of length of GRCh38 rCOR gaps")
dev.off()
