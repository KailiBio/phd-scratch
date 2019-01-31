
# -- Kaili
# This script is for making CTCF signal for all states, in order to decide and choose CTCF states.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/state_CTCF_signal/")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

################
# peak
################
peak_count = read.table("liver14.5_11state_peak.txt")
peak_count = transform(peak_count, percenatge=round(peak_count$V3/peak_count$V2,2))
colnames(peak_count) = c("state", "total", "with_peak", "percentage")
peak_count$state = factor(peak_count$state, levels = peak_count$state)

p1 = ggplot(peak_count, aes(x=state, y=percentage)) +
  geom_bar(stat="identity", fill="#f8766d", width=0.7) +
  theme_minimal() + 
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12)) +
  ylab("percentage of bins with CTCF peak") +
  labs(title="percentage of bins with CTCF peaks in 11 states")
p1

p2 = ggplot(peak_count, aes(x=state, y=total)) +
  geom_bar(stat="identity", fill="#999999", width=0.7) + 
  geom_text(aes(label=total), size=4, position=position_dodge(width=0.9)) +
  theme_minimal() + 
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12)) +
  ylab("number of bins") +
  labs(title="number of bins in 11 states")
p2
################
# motif
################
motif_count = read.table("liver14.5_11state_motif_0.5.txt")
motif_count = transform(motif_count, percenatge=round(motif_count$V3/motif_count$V2,2))
colnames(motif_count) = c("state", "total", "with_motif", "percentage")
motif_count$state = factor(motif_count$state, levels = motif_count$state)

p3 = ggplot(motif_count, aes(x=state, y=percentage)) +
  geom_bar(stat="identity", fill="#00bfc4", width=0.7) +
  theme_minimal() + 
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12)) +
  ylab("percentage of bins with CTCF motif") +
  labs(title="percentage of bins with CTCF motif in 11 states")
p3
################
# ave signal
################
s = c(21,29,30,40,38,42,27,11,28,26,31)
data = read.table("liver14.5_state_CTCF_signal_ave_sorted.txt")
colnames(data) = c("state", "ave")
rownames(data) = data$state
dat = data[as.character(s),]
dat$state = factor(dat$state, levels = dat$state)

p4 = ggplot(dat,aes(x=state, y=ave)) + 
  geom_bar(stat="identity",fill="#B79F00", width=0.7 ) +
  theme_minimal() +
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12)) +
  ylab("average of CTCF signal") +
  labs(title="average of CTCF signal in 11 states")
p4
################
# CTCF signal
################
matrix = matrix(,ncol=4)
s = c(21,29,30,40,38,42,27,11,28,26,31)
for(i in 1:11){
  a = read.table(paste("state_",as.character(s[i]),"_CTCF_signal.txt",sep=""))
  matrix = rbind(matrix, a)
}
matrix2 = matrix[-c(1),]
colnames(matrix2) = c("state", "sample", "bin", "signal")
matrix3 = matrix2[matrix2$sample=="liver_14.5",]
matrix3$state = factor(matrix3$state, levels = s)


p5 = ggplot(matrix3, aes(x=reorder(state, -signal, mean), y=signal, fill=state, 
                         group=state)) + 
  geom_boxplot(width=0.5, outlier.size=0, fill = "#D55E00") +
  theme_minimal() + 
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12)) + 
  labs(title="CTCF signal of bins that in different states", 
       x="statets", y="CTCF signal")
p5
################
# median of CTCF signal
################
s = c(21,29,30,40,38,42,27,11,28,26,31)
m = cbind(as.numeric(as.vector(dat$state)),rep(0,11))
for(i in 1:11){
  m[i,2] = as.numeric(median(matrix3[matrix3$state==s[i],]$signal)[1])
}
m = data.frame(m)
colnames(m) = c("state", "median")
m$state = factor(s, levels = s)

p6 = ggplot(m,aes(x=state, y=median)) + 
  geom_bar(stat="identity",fill="#E69F00", width=0.7) +
  theme_minimal() +
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12)) +
  ylab("median of CTCF signal") +
  labs(title="median of CTCF signal in 11 states")
p6
################
p = grid.arrange(p5, p1, p4, p3, p6, p2, ncol=2)

ggsave("CTCF_11states_barplot.pdf",p,width=12, height = 12)
ggsave("CTCF_11states_barplot.png",p,width=12, height = 12)

################
library(reshape2)
peak_count_all = read.table("liver14.5_11state_peak.txt")
peak_count_32 = read.table("liver14.5_11state_peak32.txt", row.names = 1)
peak_count = transform(peak_count_all, 
                           high_32 = peak_count_32[as.character(peak_count_all[,1]),2])
peak_count = transform(peak_count, low_32 = peak_count$V3-peak_count$high_32)
peak_count = transform(peak_count, p1 = peak_count$V3/peak_count$V2,
                       p2 = peak_count$high_32/peak_count$V2,
                       p3 = peak_count$low_32/peak_count$V2)
colnames(peak_count) = c("state","total","with_peaks","high32","low32",
                         "percentage","percentage_of_high_signal_peaks",
                         "percentage_of_low_signal_peaks")

peak_matrix = melt(peak_count[,c(1,6:8)], id.vars = "state")
peak_matrix$state = factor(peak_matrix$state, 
                           levels = c(21,29,30,40,38,42,27,11,28,26,31))

ggplot(peak_matrix, aes(x = state, y=value, fill = variable)) + 
  geom_bar(stat = "identity",position = "dodge") +
  theme_minimal() + 
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12),
        legend.position = c(0.8,0.85)) +
  ylab("percentage of bins with high CTCF peak\n") +
  labs(title="percentage of bins with CTCF peaks in 11 states")

ggsave("CTCF_peak_percentage_11states_32.pdf")
ggsave("CTCF_peak_percentage_11states_32.png")
################
