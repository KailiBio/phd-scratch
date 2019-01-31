
# -- Kaili
# This script is for making CTCF signal for all states, in order to decide and choose CTCF states.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/")
library(ggplot2)
library(RColorBrewer)


##################
matrix = matrix(,ncol=4)
for(i in 1:42){
  a = read.table(paste("./state_CTCF_signal/state_",as.character(i),"_CTCF_signal.txt",sep=""))
  matrix = rbind(matrix, a)
}

matrix2 = matrix[-c(1),]
colnames(matrix2) = c("state", "sample", "bin", "signal")
ggplot(matrix2, aes(x=state, y=signal, fill=state, group=state)) + 
  geom_boxplot(width=0.5, outlier.size=0, fill = "#999999") +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) + 
  labs(title="CTCF signal of bins that in different states", x="statets", y="CTCF signal")

ggsave("CTCF_signal_all_states_boxplot.pdf", width=7)
ggsave("CTCF_signal_all_states_boxplot.png", width=7)



#######################
# barplot
#######################
data = read.table("./state_CTCF_signal/state_CTCF_signal_ave_sorted.txt")
colnames(data) = c("state", "ave")
data$state = factor(data$state, levels = data$state)

ggplot(data,aes(x=state, y=ave)) + geom_bar(stat="identity") +
  theme_minimal() +
  theme(title=element_text(face="bold",size=12)) +
  ylab("average of CTCF signal")
ggsave("ave_CTCF_signal_all_states_barplot.pdf", width=7)
ggsave("ave_CTCF_signal_all_states_barplot.png", width=7)


#######################

