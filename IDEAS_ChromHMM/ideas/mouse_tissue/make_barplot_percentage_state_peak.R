
# -- Kaili
# This script is for making barplot to show the percentage of bin and peaks.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/peak_percentage/")

library(ggplot2)

##############
sample_list = c("liver_14.5", "lung_14.5", "liver_0", "lung_0", "forebrain_0", "midbrain_0", 
                "hindbrain_0", "heart_0", "intestine_0", "stomach_0", "kidney_0")

pdf("bin-ratio-have-peak.pdf", width=10, onefile = TRUE)
for(i in 1:length(sample_list)){
  sample=sample_list[i]
  print(sample)
  #
  data1 = read.table(paste("bin-ratio-have-peak_9impute11_",sample,".txt",sep=""))
  data2 = read.table(paste("bin-ratio-have-peak_9to11_",sample,".txt",sep=""))
  colnames(data1) = colnames(data2) = c("state","num","percentage")
  data = rbind(data1, data2)
  data = transform(data, type=c(rep("impute",47), rep("no-impute",47)))
  #
  p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
    theme_minimal() + theme(title=element_text(face="bold",size=12)) +
    labs(title=paste("percentage of bins that are peaks in ",sample,sep="")) + 
    coord_cartesian(ylim=c(0,0.8)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
  print(p)
}
dev.off()

pdf("percentage_peaks_in_state.pdf", width=10, onefile = TRUE)
for(i in 1:length(sample_list)){
  sample=sample_list[i]
  print(sample)
  #
  data1 = read.table(paste("peak-ratio-in-state_9impute11_",sample,".txt",sep=""))
  data2 = read.table(paste("peak-ratio-in-state_9to11_",sample,".txt",sep=""))
  colnames(data1) = colnames(data2) = c("state","num","percentage")
  data = rbind(data1, data2)
  data = transform(data, type=c(rep("impute",47), rep("no-impute",47)))
  #
  p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
    theme_minimal() + theme(title=element_text(face="bold",size=12)) +
    labs(title=paste("percentage of peaks in each state, ",sample,sep="")) + 
    coord_cartesian(ylim=c(0,0.35)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
  print(p)
}
dev.off()

##############
# 9to11 9impute11 9impute66
##############
pdf("compare_3groups_barplot.pdf", width=10, onefile = TRUE)
##
sample="liver_14.5"
#
data1 = read.table(paste("peak-ratio-in-state_9to11_",sample,".txt",sep=""))
data2 = read.table(paste("peak-ratio-in-state_9impute11_",sample,".txt",sep=""))
data3 = read.table(paste("peak-ratio-in-state_9impute66_",sample,".txt",sep=""))
colnames(data1) = colnames(data2) = colnames(data3) = c("state","num","percentage")
data = rbind(data1, data2, data3)
data = transform(data, type=c(rep("no-impute",47), rep("9impute11",47), rep("9impute66",47)))
p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) +
  labs(title=paste("percentage of peaks in each state, ",sample,sep="")) + 
  coord_cartesian(ylim=c(0,0.35)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
print(p)
#
data1 = read.table(paste("bin-ratio-have-peak_9to11_",sample,".txt",sep=""))
data2 = read.table(paste("bin-ratio-have-peak_9impute11_",sample,".txt",sep=""))
data3 = read.table(paste("bin-ratio-have-peak_9impute66_",sample,".txt",sep=""))
colnames(data1) = colnames(data2) = colnames(data3) = c("state","num","percentage")
data = rbind(data1, data2, data3)
data = transform(data, type=c(rep("no-impute",47), rep("9impute11",47), rep("9impute66",47)))
p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) +
  labs(title=paste("percentage of bins that are peaks in ",sample,sep="")) + 
  coord_cartesian(ylim=c(0,0.8)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
print(p)
##
sample="lung_14.5"
#
data1 = read.table(paste("peak-ratio-in-state_9to11_",sample,".txt",sep=""))
data2 = read.table(paste("peak-ratio-in-state_9impute11_",sample,".txt",sep=""))
data3 = read.table(paste("peak-ratio-in-state_9impute66_",sample,".txt",sep=""))
colnames(data1) = colnames(data2) = colnames(data3) = c("state","num","percentage")
data = rbind(data1, data2, data3)
data = transform(data, type=c(rep("no-impute",47), rep("9impute11",47), rep("9impute66",47)))
p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) +
  labs(title=paste("percentage of peaks in each state, ",sample,sep="")) + 
  coord_cartesian(ylim=c(0,0.35)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
print(p)
#
data1 = read.table(paste("bin-ratio-have-peak_9to11_",sample,".txt",sep=""))
data2 = read.table(paste("bin-ratio-have-peak_9impute11_",sample,".txt",sep=""))
data3 = read.table(paste("bin-ratio-have-peak_9impute66_",sample,".txt",sep=""))
colnames(data1) = colnames(data2) = colnames(data3) = c("state","num","percentage")
data = rbind(data1, data2, data3)
data = transform(data, type=c(rep("no-impute",47), rep("9impute11",47), rep("9impute66",47)))
p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) +
  labs(title=paste("percentage of bins that are peaks in ",sample,sep="")) + 
  coord_cartesian(ylim=c(0,0.8)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
print(p)
dev.off()

##############