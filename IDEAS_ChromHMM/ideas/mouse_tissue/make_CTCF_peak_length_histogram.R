
# -- Kaili
# This script is for make CTCF peak length histogram.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/ctcf_peaks/")

library(ggplot2)

###########
pdf("CTCF_peak_length_histogram.pdf")

sample_list = c("forebrain_0" , "midbrain_0", "hindbrain_0", "heart_0", "intestine_0", "kidney_0", "liver_0",
                "liver_14.5", "lung_0", "lung_14.5", "stomach_0")

for(i in 1:11){
  sample = sample_list[i]
  dat = read.table(paste(sample, "_ctcf_peak_sorted.bed", sep=""))
  colnames(dat) = c("chr", "start", "end", "peak") 
  
  p = ggplot(dat, aes(x = end-start)) + geom_histogram(binwidth = 1, aes(y=..density..)) +
    theme_minimal() + labs(title=sample) + xlab("CTCF peak length (bp)")
  print(p)
}

dev.off()



