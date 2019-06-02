
# -- Kaili
# This script is for making barplot to show the percentage of bin and peaks between rep1&rep2.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/rep2/peak_percentage/")

library(ggplot2)

sample_list = c("liver_14.5", "lung_14.5", "liver_0", "lung_0", "forebrain_0", "midbrain_0", 
                "hindbrain_0", "heart_0", "intestine_0", "stomach_0", "kidney_0")
##############
# percentage of bins have peaks
##############

pdf("bin-ratio-have-peak_rep1-2.pdf", width=10, onefile = TRUE)
for(i in 1:length(sample_list)){
  sample=sample_list[i]
  print(sample)
  #
  data1 = read.table(paste("bin-ratio-have-peak_rep1_",sample,".txt",sep=""))
  data2 = read.table(paste("bin-ratio-have-peak_rep1-impute-rep1_",sample,".txt",sep=""))
  data3 = read.table(paste("bin-ratio-have-peak_rep1-impute-rep2_",sample,".txt",sep=""))
  data = rbind(rbind(data1, data2), data3)
  data = transform(data, type=c(rep("rep1",43), rep("rep1-impute-rep1",43), rep("rep1-impute-rep2",43)))
  colnames(data) = c("state","num","percentage","type")
  #
  p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
    theme_minimal() + theme(title=element_text(face="bold",size=12)) +
    labs(title=paste("percentage of bins that are peaks in ",sample,sep="")) + 
    coord_cartesian(ylim=c(0,0.8)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
  print(p)
}
dev.off()

##############
# percentage of peaks in each state
##############
pdf("percentage_peaks_in_state_rep1-2.pdf", width=10, onefile = TRUE)
for(i in 1:length(sample_list)){
  sample=sample_list[i]
  print(sample)
  #
  data1 = read.table(paste("peak-ratio-in-state_rep1_",sample,".txt",sep=""))
  data2 = read.table(paste("peak-ratio-in-state_rep1-impute-rep1_",sample,".txt",sep=""))
  data3 = read.table(paste("peak-ratio-in-state_rep1-impute-rep2_",sample,".txt",sep=""))
  data = rbind(data1, data2, data3)
  data = transform(data, type=c(rep("rep1",43), rep("rep1-impute-rep1",43), rep("rep1-impute-rep2",43)))
  colnames(data) = c("state","num","percentage","type")
  #
  p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
    theme_minimal() + theme(title=element_text(face="bold",size=12)) +
    labs(title=paste("percentage of peaks in each state, ",sample,sep="")) + 
    coord_cartesian(ylim=c(0,0.35)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
  print(p)
}
dev.off()

##############
# percentage of bins have peaks in lung14.5
##############
sample="lung_14.5"
#
data1 = read.table(paste("bin-ratio-have-peak_rep1_",sample,".txt",sep=""))
data2 = read.table(paste("bin-ratio-have-peak_rep1-impute-rep1_oneSample_",sample,".txt",sep=""))
data3 = read.table(paste("bin-ratio-have-peak_rep1-impute-rep1_oneSample_",sample,".txt",sep=""))
data = rbind(rbind(data1, data2), data3)
data = transform(data, type=c(rep("rep1",43), rep("rep1-impute-rep1",43), rep("rep1-impute-rep2",43)))
colnames(data) = c("state","num","percentage","type")
#
ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) +
  labs(title=paste("percentage of bins that are peaks in ",sample,sep="")) + 
  coord_cartesian(ylim=c(0,0.8)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))

ggsave("bin-ratio-have-peak_rep1-2_lung_14.5.pdf")
##############
##############
sample="lung_14.5"
#
data1 = read.table(paste("peak-ratio-in-state_rep1_",sample,".txt",sep=""))
data2 = read.table(paste("peak-ratio-in-state_rep1-impute-rep1_oneSample_",sample,".txt",sep=""))
data3 = read.table(paste("peak-ratio-in-state_rep1-impute-rep2_oneSample_",sample,".txt",sep=""))
data = rbind(data1, data2, data3)
data = transform(data, type=c(rep("rep1",43), rep("rep1-impute-rep1",43), rep("rep1-impute-rep2",43)))
colnames(data) = c("state","num","percentage","type")
#
p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) +
  labs(title=paste("percentage of peaks in each state, ",sample,sep="")) + 
  coord_cartesian(ylim=c(0,0.35)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))

ggsave("percentage_peaks_in_state_rep1-2_lung_14.5.pdf")
##############
# percentage of bins have peaks, given 9 samples
##############

pdf("bin-ratio-have-peak_rep1-2_given9samples.pdf", width=12, onefile = TRUE)
for(i in 1:length(sample_list)){
  sample=sample_list[i]
  print(sample)
  #
  data1 = read.table(paste("bin-ratio-have-peak_rep1_",sample,".txt",sep=""))
  data2 = read.table(paste("bin-ratio-have-peak_rep1-impute-rep1_given9samples_",sample,".txt",sep=""))
  data3 = read.table(paste("bin-ratio-have-peak_rep1-impute-rep2_given9samples_",sample,".txt",sep=""))
  data4 = read.table(paste("bin-ratio-have-peak_rep1-impute-rep2_given9samples_10k_",sample,".txt",sep=""))
  data = rbind(rbind(data1, data2), data3, data4)
  data = transform(data, type=c(rep("rep1",43), rep("rep1-impute-rep1",43), rep("rep1-impute-rep2",43), 
                                rep("rep1-impute-rep2_10k",43)))
  colnames(data) = c("state","num","percentage","type")
  #
  p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
    theme_minimal() + theme(title=element_text(face="bold",size=12)) +
    labs(title=paste("percentage of bins that are peaks in ",sample,"\ngiven_CTCF_in_9samples", sep="")) + 
    coord_cartesian(ylim=c(0,0.8)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
  print(p)
}
dev.off()

##############
# percentage of peaks in each state, 9 samples
##############
pdf("percentage_peaks_in_state_rep1-2_given9samples.pdf", width=12, onefile = TRUE)
for(i in 1:length(sample_list)){
  sample=sample_list[i]
  print(sample)
  #
  data1 = read.table(paste("peak-ratio-in-state_rep1_",sample,".txt",sep=""))
  data2 = read.table(paste("peak-ratio-in-state_rep1-impute-rep1_given9samples_",sample,".txt",sep=""))
  data3 = read.table(paste("peak-ratio-in-state_rep1-impute-rep2_given9samples_",sample,".txt",sep=""))
  data4 = read.table(paste("peak-ratio-in-state_rep1-impute-rep2_given9samples_10k_",sample,".txt",sep=""))
  data = rbind(data1, data2, data3, data4)
  data = transform(data, type=c(rep("rep1",43), rep("rep1-impute-rep1",43), rep("rep1-impute-rep2",43), 
                                rep("rep1-impute-rep2_10k",43)))
  colnames(data) = c("state","num","percentage","type")
  #
  p = ggplot(data, aes(x=state, y=percentage, fill=type)) + geom_bar(stat="identity", position="dodge") +
    theme_minimal() + theme(title=element_text(face="bold",size=12)) +
    labs(title=paste("percentage of peaks in each state, ",sample,"\ngiven_CTCF_in_9samples",sep="")) + 
    coord_cartesian(ylim=c(0,0.45)) + scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
  print(p)
}
dev.off()

#########
m = data.frame(cbind(c(77172, 13728, 7203, 5327, 9145, 6096, 6452, 4078, 5212, 5125, 3224),
                 c(79115, 11883, 5876, 3174, 5956, 3748, 3619, 3344, 2741, 3165, 1296),
                 c(83605, 18368, 3607, 322807, 2483, 4005, 4288, 3092, 2152, 3588, 1160),
                 c(83871, 17228, 3291, 323308, 2483, 4095, 4254, 2811, 2022, 3807, 1228)))
colnames(m) = c("rep1", "rep1-impute-rep1", "rep1-impute-rep2", "rep1-impute-rep2_10k")
rownames(m) = c("11","21","26","27","28","29","30","31","38","40","42")
m2 = melt(m, id.vars = NULL)
m3 = cbind(m2, rep(c("11","21","26","27","28","29","30","31","38","40","42"),4))
colnames(m3) = c("type", "count", "state")

ggplot(m3, aes(x=state, y=log10(count), group=type, fill=type)) + 
  geom_bar(stat="identity",position=position_dodge())
ggsave("number_of_CTCF_states.pdf")
ggsave("number_of_CTCF_states.png")
##############