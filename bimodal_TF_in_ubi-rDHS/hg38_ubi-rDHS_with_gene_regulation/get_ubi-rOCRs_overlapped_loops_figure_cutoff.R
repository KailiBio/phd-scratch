
# -- Kaili
# This script is for maing figures for ubi-rOCRs overlapped loops. (running cut-off)

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/loop/")

data = read.table("GRCh38_GM12878_ubi-rOCRs_overlapped_runing_cutoff.txt")
data = transform(data, loops = data[,4]+data[,5], percent = (data[,4]+data[,5])/data[,2])
colnames(data) = c("num_of_loops", "loop_hg38", "ubi-rOCRs", "single", "both", "loops", "percent")


library(ggplot2)
ggplot(data, aes(x=num_of_loops, y=percent)) + geom_line() + geom_point() +
  labs(title="percentage of loops that overlapped with ubi-rOCRs") +
  xlab("top loops (fdr)") + ylab("percentage of loops with ubi-rOCRs")
ggsave("percentage_ubi-rOCRs_overlapped_loops.pdf")

