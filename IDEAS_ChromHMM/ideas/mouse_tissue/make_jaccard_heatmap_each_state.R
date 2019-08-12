
# -- Kaili
# This script is for making Jaccard heapmap for given union/overlapped matrix. 

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/state_jaccard/")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(pheatmap)

################
pdf("state_between_all_samples.pdf")
for(i in 0:42){
  data1 = read.table(paste("state_",as.character(i),"_overlapped.txt",sep=""),
                     header=TRUE)
  rownames(data1) = colnames(data1)
  data2 = read.table(paste("state_",as.character(i),"_union.txt",sep=""),
                     header = TRUE)
  rownames(data2) = colnames(data2)
  #
  pheatmap(data1/data2, display_numbers = TRUE, 
           main=paste("state ",as.character(i),sep=""), 
           breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
           col = brewer.pal(10,"RdYlBu")[10:1])
}
dev.off()




