
# -- Kaili
# This script is for making RAMPAGE peak RPM histogram.
# purpose: for choosing cut-off to filter RAMPAGE peaks.

setwd("/data/zusers/fankaili/ccre/hg38_ubi-rDHS/")
# setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/peak_shape/")

library(ggplot2)
library(gridExtra)

####
list = read.table("hg38_RAMPAGE_list.txt")

pdf("./promoter_shape/RAMPAGE_peak_RPM_histogram.pdf", width = 9, height = 4)
for(i in 1:155){
  sample = list[i,4]
  exp_id = list[i,1]
  #
  dat = read.table(paste("/data/zusers/zhangx/projects/RAMPAGE_peaks/hg38/", exp_id, "_rampage_peaks.txt", sep=""))
  colnames(dat)[13] = "RPM"
  #
  p1 = ggplot(dat, aes(x = RPM)) + geom_histogram(binwidth = 1) +
    theme_minimal() + theme(title=element_text(size=12, face="bold")) +
    coord_cartesian(xlim=c(0,25)) +
    labs(title = sample, subtitle = exp_id) +
    scale_y_continuous(labels = scales::scientific)
  p2 = ggplot(dat, aes(x = RPM)) + geom_histogram(binwidth = 1) +
    theme_minimal() + theme(title=element_text(size=12, face="bold")) +
    coord_cartesian(xlim=c(0,25), ylim=c(0,50000)) +
    labs(subtitle = "ZoomIn: RPM=c(0,25)") +
    scale_y_continuous(labels = scales::scientific)
  p = grid.arrange(p1, p2, ncol=2)
  print(p)
}
dev.off()

