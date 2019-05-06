
# -- Kaili
# This script is for making length histogram of merged regions.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/regions_for_testing/")

library(ggplot2)


##########
sample_list = c("forebrain", "midbrain", "hindbrain", "neural-tube", "heart", "facial", "limb", "liver")

pdf("merged-bins_length-histogram.pdf")
for(i in 1:length(sample_list)){
  tissue=sample_list[i]
  data = read.table(paste("./table/e11.5_",tissue,"_merged_regions_Table.txt", sep=""))
  #
  p = ggplot(data, aes(x=V2/1000))+ 
    geom_histogram(binwidth=100,aes(y=..density..)) +
    theme_minimal() + theme(title= element_text(face="bold", size=12)) +
    xlab("size (kb)") + 
    labs(title=paste(tissue,"\nlength distribution of after-merged regions\n(bin-size=100kb)",sep=""))
  print(p)
}
dev.off()
