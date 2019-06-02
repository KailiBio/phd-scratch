
# -- Kaili
# This script is for making scatter plot for comparing signal in rep1 and rep2.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/rep2/random_10k_signal/")

library(ggplot2)
library(data.table)

#############
mark_list=c("H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3", "H3K9ac", "H3K27me3", "H3K27ac", "H3K36me3", "ATAC", 
            "CTCF", "DNAme")

pdf("random_10k_signal_reps_scatter.pdf")
for(i in 1:11){
  mark = mark_list[i]
  file_name = paste("kidney_0_",mark, "_random_10k.txt", sep="")
  dat = fread(file_name)
  colnames(dat) = c("id", "rep2", "rep1")
  m = ceiling(max(dat$rep2, dat$rep1))
  p = ggplot(dat, aes(x=rep1, y=rep2)) + geom_point() +
    theme_minimal() +  theme(title = element_text(size=12, face="bold")) +
    labs(title=mark) + coord_cartesian(ylim=c(0,m), xlim=c(0,m)) +
    geom_abline(slope=1, intercept = 0, col="red", linetype = "dashed")
  print (p)
}
dev.off()






