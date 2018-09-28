
# -- Kaili
# This script is for making histogram of numbers of high DNase samples for each rOCRs, then get the cutoff for ubi-rOCRs.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/")

data = read.table("EDGE-DNase-Biosample-Counts.txt", row.names = 1)

pdf("hg38_rOCRs_high_DNase_histogram.pdf", width=8, height = 6)
hist(data[,1], breaks = max(data[,1]), col="grey", border = "grey", freq = F, main = "hg38 rOCRs",
     xlab="number of high DNase samples")
lines(c(580, 580), c(0,1), col="red", lwd=2, lty=2)
text(550,0.2,"580",col="red")
dev.off()

pdf("hg38_rOCRs_high_DNase_histogram2.pdf", width=8, height = 6)
hist(data[,1], breaks = max(data[,1]), col="grey", border = "grey", freq = F, ylim=c(0,0.005), 
     main = "hg38 rOCRs", xlab="number of high DNase samples")
lines(c(585, 585), c(0,1), col="red", lwd=2, lty=2)
text(550,0.004,"580",col="red")
dev.off()

####
data2=read.table("ccRE-DNase-Biosample-Counts.txt", row.names = 1)

pdf("hg19_rDHSs_high_DNase_histogram.pdf", width=8, height = 6)
hist(data2[,1], breaks = max(data2[,1]), col="grey", border = "grey", freq = F, main = "hg19 rDHSs",
     xlab="number of high DNase samples")
lines(c(450, 450), c(0,1), col="red", lwd=2, lty=2)
text(430,0.3,"450",col="red")
dev.off()

pdf("hg19_rDHSs_high_DNase_histogram2.pdf", width=8, height = 6)
hist(data2[,1], breaks = max(data2[,1]), col="grey", border = "grey", freq = F, main = "hg19 rDHSs",
     xlab="number of high DNase samples", ylim=c(0,0.005))
lines(c(450, 450), c(0,1), col="red", lwd=2, lty=2)
text(430,0.004,"450",col="red")
dev.off()
