
# -- Kaili
# This script is for making barplot of CTCF peak consistency.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/ctcf_peaks/")

library(ggplot2)

###########

pdf("CTCF_peak_consistent_barplot.pdf")

sample_list = c("forebrain_0" , "midbrain_0", "hindbrain_0", "heart_0", "intestine_0", "kidney_0", "liver_0",
                "lung_0", "stomach_0")

for(i in 1:9){
  sample = sample_list[i]
  dat = read.table(paste(sample, "_specific_CTCFpeaks_count_all.txt", sep=""))
  colnames(dat) = c("peak", "num") 
  dd = data.frame(table(dat$num))
  
  p = ggplot(dd, aes(x = Var1, y=Freq)) + geom_bar(stat = "identity") +
    theme_minimal() + labs(title=sample) + xlab("number of tissue have this peak") +
    ylab("number of CTCF peaks") + geom_text(aes(label=Freq), vjust=-0.4)
  print(p)
}

dev.off()

#################
# for all 11 CTCF samples
#################
pdf("CTCF_peak_consistent_barplot_all11samples.pdf")

sample_list = c("forebrain_0" , "midbrain_0", "hindbrain_0", "heart_0", "intestine_0", "kidney_0", "liver_0",
                "lung_0", "stomach_0", "liver_14.5", "lung_14.5")

for(i in 1:11){
  sample = sample_list[i]
  dat = read.table(paste(sample, "_specific_11CTCFpeaks_count_all.txt", sep=""))
  colnames(dat) = c("peak", "num") 
  dd = data.frame(table(dat$num))
  
  p = ggplot(dd, aes(x = Var1, y=Freq)) + geom_bar(stat = "identity") +
    theme_minimal() + labs(title=sample) + xlab("number of tissue have this peak") +
    ylab("number of CTCF peaks") + geom_text(aes(label=Freq), vjust=-0.4)
  print(p)
}

dev.off()


#################

