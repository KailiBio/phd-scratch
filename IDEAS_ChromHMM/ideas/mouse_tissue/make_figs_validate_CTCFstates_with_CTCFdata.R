
# -- Kaili
# This script is for making figures for validate_CTCFstates_with_CTCFdata.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/ctcfstate_with_ctcfdata/")

library(ggplot2)
library(gridExtra)

##################
# signal of lung_0 in state 12
##################
signal_state12 = read.table("lung0_state12_signal.txt")
colnames(signal_state12) = c("id", "signal", "mark")

p1 = ggplot(signal_state12, aes(x = mark, y = signal, group = mark, fill = mark)) +
  geom_violin() + geom_boxplot(width=0.1, outlier.size = 0) +
  xlab("") + theme_minimal() +
  theme(title = element_text(face="bold", size=12)) +
  labs(title="signal of 11 marks in state12 (CTCF state)")
p1

p2 = ggplot(signal_state12, aes(x = mark, y = log10(signal+1e-5), group = mark, fill = mark)) +
  geom_violin() + geom_boxplot(width=0.1, outlier.size = 0) +
  xlab("") + theme_minimal() +
  theme(title = element_text(face="bold", size=12))
p2

p = grid.arrange(p1,p2, ncol=1)
p
ggsave("signal_of_lung0_in_state12_violin.pdf", p, width = 12, height=8)
ggsave("signal_of_lung0_in_state12_violin.png", p, width = 12, height=8)


#####
signal_allstates = read.table("lung0_allstates_signal.txt")
colnames(signal_allstates) = c("id", "signal", "mark")

p1 = ggplot(signal_allstates, aes(x = mark, y = signal, group = mark, fill = mark)) +
  geom_violin() + geom_boxplot(width=0.1, outlier.size = 0) +
  xlab("") + theme_minimal() +
  theme(title = element_text(face="bold", size=12)) +
  labs(title="signal of 11 marks")
p1

p2 = ggplot(signal_allstates, aes(x = mark, y = log10(signal+1e-5), group = mark, fill = mark)) +
  geom_violin() + geom_boxplot(width=0.1, outlier.size = 0) +
  xlab("") + theme_minimal() +
  theme(title = element_text(face="bold", size=12))
p2

p = grid.arrange(p1,p2, ncol=1)
p
ggsave("signal_of_lung0_in_allstates_violin.pdf", p, width = 12, height=8)
ggsave("signal_of_lung0_in_allstates_violin.png", p, width = 12, height=8)


##################
# number of CTCF states in different samples
##################
state_num = read.table("sample_state12_count.txt")
colnames(state_num) = c("sample", "num")
a = t(matrix(unlist(strsplit(as.vector(state_num$sample), split="_")), nrow=2))
state_num = data.frame(cbind(state_num, a))
colnames(state_num) = c("sample", "num", "tissue","timepoint")

ggplot(state_num, aes(x = timepoint, y = num, fill = timepoint)) + 
  geom_bar(stat = "identity") + facet_wrap(~tissue) +
  theme_minimal() +
  theme(title = element_text(face="bold", size=12)) +
  labs(title="number of state12 (CTCF state) in each sample") +
  xlab("") + ylab("number of state12 bins")

ggsave("number_of_state12_in_each_sample_barplot.pdf", width=12)
ggsave("number_of_state12_in_each_sample_barplot.png", width=12)

##################
# regression
##################
signal_matrix = read.table("lung0_signal_random100k.txt", header = TRUE)

model1 = lm(signal_matrix$CTCF~signal_matrix$ATAC + signal_matrix$DNAme + signal_matrix$H3K4me1 +
             signal_matrix$H3K4me2 + signal_matrix$H3K4me3 + signal_matrix$H3K9me3 +
             signal_matrix$H3K9ac + signal_matrix$H3K27me3 + signal_matrix$H3K27ac +
             signal_matrix$H3K36me3)
summary(model1)$adj.r.squared

# model2 = lm(log(signal_matrix$CTCF+1e-5)~log(signal_matrix$ATAC+1e-5) + 
#               log(signal_matrix$DNAme+1e-5) + log(signal_matrix$H3K4me1+1e-5) +
#               log(signal_matrix$H3K4me2+1e-5) + log(signal_matrix$H3K4me3+1e-5) + 
#               log(signal_matrix$H3K9me3+1e-5) + log(signal_matrix$H3K9ac+1e-5) + 
#               log(signal_matrix$H3K27me3+1e-5) + log(signal_matrix$H3K27ac+1e-5) +
#               log(signal_matrix$H3K36me3+1e-5))
# summary(model2)$adj.r.squared

model3 = lm(signal_matrix$H3K4me3~signal_matrix$ATAC + signal_matrix$DNAme + signal_matrix$H3K4me1 +
              signal_matrix$H3K4me2 + signal_matrix$H3K9me3 +
              signal_matrix$H3K9ac + signal_matrix$H3K27me3 + signal_matrix$H3K27ac +
              signal_matrix$H3K36me3)
summary(model3)$adj.r.squared

# model4 = lm(log(signal_matrix$H3K4me3+1e-5) ~log(signal_matrix$ATAC+1e-5) + 
#               log(signal_matrix$DNAme+1e-5) + log(signal_matrix$H3K4me1+1e-5) +
#               log(signal_matrix$H3K4me2+1e-5) + 
#               log(signal_matrix$H3K9me3+1e-5) + log(signal_matrix$H3K9ac+1e-5) + 
#               log(signal_matrix$H3K27me3+1e-5) + log(signal_matrix$H3K27ac+1e-5) +
#               log(signal_matrix$H3K36me3+1e-5))
# summary(model4)$adj.r.squared

model5 = lm(signal_matrix$H3K27ac~signal_matrix$ATAC + signal_matrix$DNAme + signal_matrix$H3K4me1 +
              signal_matrix$H3K4me2 + signal_matrix$H3K4me3 + signal_matrix$H3K9me3 +
              signal_matrix$H3K9ac + signal_matrix$H3K27me3 + 
              signal_matrix$H3K36me3)
summary(model5)$adj.r.squared
##################
# ARI in HOX region
##################


##################

