
# -- Kaili
# This scirpt is for making histogram of CTCF signal of different Hi-C peal loci.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/")

data = read.table("GRCh38_peakLoci_GM12878_CTCF_z-score_marked.txt")
colnames(data) = c("chr", "s", "e", "peak", "overlapped", "CTCF_signal")
data2 = read.table("GRCh38_peakLoci_GM12878_CTCF_signal_marked.txt")
colnames(data2) = c("chr", "s", "e", "peak", "overlapped", "CTCF_signal")

library(ggplot2)

###################
# CTCF zscore
###################
a = data[data$overlapped=="ubi-rOCR",]$CTCF_signal
b = data[data$overlapped=="active-rOCR",]$CTCF_signal
t.test(a,b)$p.value
c = data[data$overlapped!="remaining",]$CTCF_signal
d = data[data$overlapped=="remaining",]$CTCF_signal
t.test(c,d)$p.value

ggplot(data, aes(x=CTCF_signal, color=overlapped)) + geom_density(size=1.2) + 
  xlab("CTCF signal (z-score)") +
  labs(title="CTCF signal of peak loci in Hi-C loops")
ggsave("GRCh38_peakLoci_GM12878_CTCF_histogram.pdf")

###################
# CTCF signal
###################
ggplot(data2, aes(x=CTCF_signal, color=overlapped)) + geom_density(size=1.2) + 
  xlab("CTCF signal (raw signal)") +
  labs(title="CTCF signal of peak loci in Hi-C loops")
ggsave("GRCh38_peakLoci_GM12878_CTCF_histogram_rawsignal.pdf")

###################
# CTCF signal in rOCRs
###################
setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/loop/")
dat = read.table("GRCh38_GM12878_peakloci_overlapped_rOCR_marked_signal.txt")
colnames(dat) = c("chr","s","e","ID","type","signal")

ggplot(dat, aes(x=signal, color=type)) + geom_density(size=1.2) + 
  xlab("CTCF signal (z-score)") +
  labs(title="CTCF signal of peak loci overlapped OCRs")
ggsave("GRCh38_GM12878_peakLoci_overlapped_OCR_CTCF_histogram_z-score.pdf")

###################