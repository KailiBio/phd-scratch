
# -- Kaili
# This script is for making boxplot for ARI matrix between imputation and real run.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/ARI_between_runs/")

library(ggplot2)
library(reshape2)
library(gridExtra)
#################
# 8sample model
#################
between_imputation_8 = melt(read.table("ARI_8sample_imputation_reproducibility.txt"))
real_data_8 = melt(read.table("ARI_8sample_withdata_reproducibility.txt"))
impuatation_real_8 = melt(read.table("ARI_8sample_imputation_validation1.txt"))
impuatation2_real2_8 = melt(read.table("ARI_8sample_imputation_validation2.txt"))

p1 = ggplot(between_imputation_8, aes(y = value, group = variable, fill = variable)) + geom_boxplot() +
  theme_minimal() + labs(title="ARI between two imputation runs\n(imputation across tissue)") +
  coord_cartesian(ylim = c(0,1)) +
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank()) +
  ylab("ARI")
p1

p2 = ggplot(real_data_8, aes(y = value, group = variable, fill = variable)) + geom_boxplot() +
  theme_minimal() + labs(title="ARI between two runs with real data") +
  coord_cartesian(ylim = c(0,1)) +
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank()) + ylab("ARI")
p2

p3 = ggplot(impuatation_real_8, aes(y = value, group = variable, fill = variable)) + geom_boxplot() +
  theme_minimal() + labs(title="ARI between imputation and real data\n(imputation across tissue)") +
  coord_cartesian(ylim = c(0,1)) +
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank()) + 
  ylab("ARI")
p3

p4 = ggplot(impuatation2_real2_8, aes(y = value, group = variable, fill = variable)) + geom_boxplot() +
  theme_minimal() + labs(title="ARI between imputation and real data\n(imputation across tissue)") +
  coord_cartesian(ylim = c(0,1)) +
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank()) + ylab("ARI")
p4

p = grid.arrange(p1, p2, p3, p4, ncol = 2)
p
ggsave("ARI_boxplot_8sample.pdf", p, width = 12, height=8)
ggsave("ARI_boxplot_8sample.png", p, width = 12, height=8)
#################
# 9sample model
#################
between_imputation_9 = melt(read.table("ARI_9sample_imputation_reproducibility.txt"))
real_data_9 = melt(read.table("ARI_9sample_withdata_reproducibility.txt"))
impuatation_real_9 = melt(read.table("ARI_9sample_imputation_validation1.txt"))
impuatation2_real2_9 = melt(read.table("ARI_9sample_imputation_validation2.txt"))

p1 = ggplot(between_imputation_9, aes(y = value, group = variable, fill = variable)) + geom_boxplot() +
  theme_minimal() + labs(title="ARI between two imputation runs\n(imputation across time points)") +
  coord_cartesian(ylim = c(0,1)) +
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank()) +
  ylab("ARI")
p1

p2 = ggplot(real_data_9, aes(y = value, group = variable, fill = variable)) + geom_boxplot() +
  theme_minimal() + labs(title="ARI between two runs with real data") +
  coord_cartesian(ylim = c(0,1)) +
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank()) + ylab("ARI")
p2

p3 = ggplot(impuatation_real_9, aes(y = value, group = variable, fill = variable)) + geom_boxplot() +
  theme_minimal() + labs(title="ARI between imputation and real data\n(imputation across time points)") +
  coord_cartesian(ylim = c(0,1)) +
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank()) + 
  ylab("ARI")
p3

p4 = ggplot(impuatation2_real2_9, aes(y = value, group = variable, fill = variable)) + geom_boxplot() +
  theme_minimal() + labs(title="ARI between imputation and real data\n(imputation across time points)") +
  coord_cartesian(ylim = c(0,1)) +
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank()) + ylab("ARI")
p4

p = grid.arrange(p1, p2, p3, p4, ncol = 2)
p
ggsave("ARI_boxplot_9sample.pdf", p, width = 12, height=8)
ggsave("ARI_boxplot_9sample.png", p, width = 12, height=8)


#################

