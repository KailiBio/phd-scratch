
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
peak_count = read.table("liver14.5_all_states_peak.txt")
peak_count = transform(peak_count, percenatge=round(peak_count$V3/peak_count$V2,2))
colnames(peak_count) = c("state", "total", "with_peak", "percentage")
peak_count$state = factor(peak_count$state, levels = peak_count$state)

p1 = ggplot(peak_count, aes(x=state, y=percentage)) +
  geom_bar(stat="identity", fill="#f8766d", width=0.7) +
  theme_minimal() + 
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12)) +
  ylab("percentage of bins with CTCF peak") +
  labs(title="percentage of bins with CTCF peaks in all states")
p1

p2 = ggplot(peak_count, aes(x=state, y=total/100000)) +
  geom_bar(stat="identity", fill="#999999", width=0.7) + 
  theme_minimal() + 
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12)) +
  ylab("number of bins(*10e5)") +
  labs(title="number of bins in all states")
p2

################
motif_count = read.table("liver14.5_all_states_motif.txt")
motif_count = transform(motif_count, percenatge=round(motif_count$V3/motif_count$V2,2))
colnames(motif_count) = c("state", "total", "with_motif", "percentage")
motif_count$state = factor(motif_count$state, levels = motif_count$state)

p3 = ggplot(motif_count, aes(x=state, y=percentage)) +
  geom_bar(stat="identity", fill="#00bfc4", width=0.7) +
  theme_minimal() + 
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12)) +
  ylab("percentage of bins with CTCF motif") +
  labs(title="percentage of bins with CTCF motif in all states")
p3

################
data = read.table("liver14.5_state_CTCF_signal_ave_sorted.txt")
colnames(data) = c("state", "ave")


p4 = ggplot(data,aes(x=state, y=ave)) + 
  geom_bar(stat="identity", fill="#B79F00") +
  theme_minimal() +
  theme(title=element_text(face="bold",size=12),
        axis.text = element_text(face="bold",size=12)) +
  ylab("average of CTCF signal") +
  labs(title="average of CTCF signal in all states")
p4
################
p = grid.arrange(p4, p1, p2, p3, ncol=2)

ggsave("CTCF_all_states_barplot.pdf",p,width=18, height = 10)
ggsave("CTCF_all_states_barplot.png",p,width=18, height = 10)
