
# -- Kaili
# This script is for making sanity check scatter plots.

setwd("/data/zusers/fankaili/ccre/mm10_rnaseq/")
#setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/mm10_RNAseq/")


library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
file = args[1]
sample = args[2]

# file="tmp.txt"
# sample="ss"
########
dat = read.table(file, header = TRUE, row.names = 1)
colnames(dat) = c("spikeIn92", "id", "spikeIn96")

r = round(cor(dat$spikeIn92, dat$spikeIn96),4)
name=paste(sample, " R=", r, sep="")

ggplot(dat, aes(x=log10(spikeIn92+0.1), y=log10(spikeIn96+0.1))) + geom_point() +
  labs(title=name) + theme_minimal() + theme(title = element_text(size=12, face = "bold"))
ggsave(paste("./check/",sample,".pdf",sep=""))
