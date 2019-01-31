
# -- Kaili
# This script is for making percentage of states with CTCF peaks with running cut-off.


setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/state_CTCF_signal/percentage_peak_cutoff/")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)


#####################

pdf("percentage_of_state_with_CTCF_peak_cutoff.pdf")
s = c(21,29,30,40,38,42,27,11,28,26,31)
for(i in 1:42){
  if(i%in%s){
    color = "#f8766d"
  }else{
    color = "#00bfc4"
  }
  #
  data = read.table(paste("liver14.5_state_",as.character(i),
                          "_peak_count.txt",sep=""))
  colnames(data) = c("cutoff", "total","num")
  data = transform(data, percentage = data$num/data$total)
  data$cutoff = factor(data$cutoff, levels=data$cutoff)
  
  p = ggplot(data, aes(x=cutoff, y=percentage)) +
    geom_bar(stat="identity", fill = color) + 
    theme_minimal() + 
    theme(title=element_text(face="bold",size=12),
          axis.text.x = element_blank()) +
    ylab("percentage of bins with CTCF peak") +
    labs(title=paste("state ",as.character(i),
                     "\npercentage of bins with CTCF peak\n(in different cut-off)", 
                     sep="")) +
    geom_vline(xintercept = 32, col="grey", linetype="dashed", size=1) +
    geom_vline(xintercept = 50, col="grey", linetype="dashed", size=1) +
    coord_cartesian(ylim=c(0,0.7))
  print(p)
}
dev.off()

#####################