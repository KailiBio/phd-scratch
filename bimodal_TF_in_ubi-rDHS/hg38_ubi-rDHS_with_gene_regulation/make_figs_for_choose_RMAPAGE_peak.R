
# -- Kaili
# This script is for making figures for choose RAMPAGE peak cut-off.
 
setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/peak_shape/")


library(ggplot2)
library(reshape2)

dat = read.table("num_peak_with_diff_cutoff.txt")
colnames(dat) = c("exp_id","id1", "id2", "sample", "RPM>1", "RPM>2", "RPM>10")

###################
# boxplot
###################
d = melt(dat)
ggplot(d, aes(x=variable, y=value)) + geom_boxplot(width=0.3) +
  theme_minimal() + theme(title = element_text(size=12, face="bold")) +
  scale_y_continuous(labels = scales::scientific) +
  xlab("") + ylab("number of RAMPAGE peaks")
ggsave("num_peak_after_RPM_cutoff.pdf")
ggsave("num_peak_after_RPM_cutoff.png")

median(dat$`RPM>1`)
# 30973

mean(dat$`RPM>1`)
# 33118

###################
# histogram
###################
pdf("histogram_of_different_RPM_cutoff.pdf")

ggplot(dat, aes(x=`RPM>1`)) + geom_histogram(binwidth = 1000) +
  theme_minimal() + theme(title = element_text(size=12, face="bold")) +
  ylab("num of samples") + scale_x_continuous(labels = scales::scientific)


ggplot(dat, aes(x=`RPM>2`)) + geom_histogram(binwidth = 1000) +
  theme_minimal() + theme(title = element_text(size=12, face="bold")) +
  ylab("num of samples") + scale_x_continuous(labels = scales::scientific)

ggplot(dat, aes(x=`RPM>10`)) + geom_histogram(binwidth = 1000) +
  theme_minimal() + theme(title = element_text(size=12, face="bold")) +
  ylab("num of samples") + scale_x_continuous(labels = scales::scientific)
dev.off()

###################
