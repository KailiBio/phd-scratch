
# -- Kaili
# This script is for composition of 11states_nopeak/motif bins.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/state_CTCF_signal/")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

####################
# peak
####################
peak_com = read.table("liver14.5_11states_noPeak_composition.txt")
colnames(peak_com) = c("state", "num")
peak_com$state = factor(peak_com$state, levels = peak_com$state)

sum = sum(peak_com$num)
label = paste(as.character(round(peak_com$num/sum,3)*100),"%",sep="")

ggplot(peak_com, aes(x="", y=num, fill=state)) +geom_bar(stat="identity") +
  coord_polar(theta="y", start=0, direction = -1) +
  theme_minimal() + 
  theme(axis.title=element_blank(), axis.ticks=element_blank(),
        axis.text=element_blank()) +
  scale_fill_manual(values =brewer.pal(11,"Set3")) +
  labs(title="percentage of 11states_noPeak bins") +
  geom_text(x=1.5, y=sum(peak_com$num) - cumsum(peak_com$num) + peak_com$num/2,
            label=label)

ggsave("percentage_11states_noPeak_bins_pie.pdf")
ggsave("percentage_11states_noPeak_bins_pie.png")


####################
# motif
####################
motif_com = read.table("liver14.5_11states_noMotif_composition.txt")
colnames(motif_com) = c("state", "num")
motif_com$state = factor(motif_com$state, levels = motif_com$state)

sum = sum(motif_com$num)
label = paste(as.character(round(motif_com$num/sum,3)*100),"%",sep="")


ggplot(motif_com, aes(x="", y=num, fill=state)) +geom_bar(stat="identity") +
  coord_polar(theta="y", start=0, direction = -1) +
  theme_minimal() + 
  theme(axis.title=element_blank(), axis.ticks=element_blank(),
        axis.text=element_blank()) +
  scale_fill_manual(values =brewer.pal(11,"Set3")) +
  labs(title="percentage of 11states_noMotif bins") +
  geom_text(x=1.5, y=sum(motif_com$num) - cumsum(motif_com$num) + motif_com$num/2,
            label=label)

ggsave("percentage_11states_noMotif_bins_pie.pdf")
ggsave("percentage_11states_noMotif_bins_pie.png")

