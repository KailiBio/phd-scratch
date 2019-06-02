
# -- Kaili
# This script is for making barplot of pearson correlation between replicates.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/rep2/")

dat = read.table("rep1_rep2_pearson.txt")
colnames(dat) = c("sample", "mark", "r")

library(ggplot2)

##############
# boxplot for each mark (normal bins)
##############
ggplot(dat, aes(x=mark, y=r, group=mark)) + geom_point(size=1) +
  theme_minimal() + theme(title = element_text(face="bold", size=12)) +
  labs(title="pearson correlation of signal in rep1 vs. rep2") + xlab("")

ggsave("reps_pearson_boxplot.pdf", width=10)
ggsave("reps_pearson_boxplot.png", width=10)


##############
# boxplot for each mark (dhs bins)
##############
dat_dhs = read.table("rep1_rep2_pearson_dhs.txt")
colnames(dat_dhs) = c("sample", "mark", "r")

ggplot(dat_dhs, aes(x=mark, y=r, group=mark)) + geom_point(size=1) +
  theme_minimal() + theme(title = element_text(face="bold", size=12)) +
  labs(title="pearson correlation of signal in rep1 vs. rep2") + xlab("")

ggsave("reps_pearson_boxplot_dhs.pdf", width=10)
ggsave("reps_pearson_boxplot_dhs.png", width=10)

##############
# boxplot for each mark (normal bins, spearman)
##############
dat_spearman = read.table("rep1_rep2_spearman.txt")
colnames(dat_spearman) = c("sample", "mark", "r")

ggplot(dat_spearman, aes(x=mark, y=r, group=mark)) + geom_point(size=1) +
  theme_minimal() + theme(title = element_text(face="bold", size=12)) +
  labs(title="spearman correlation of signal in rep1 vs. rep2") + xlab("")

ggsave("reps_spearman_boxplot.pdf", width=10)
ggsave("reps_spearman_boxplot.png", width=10)

##############
# boxplot for each mark (2k bins)
##############
dat_2k = read.table("random_10k_rep1_rep2_pearson.txt")
colnames(dat_2k) = c("sample", "mark", "r")

ggplot(dat_2k, aes(x=mark, y=r, group=mark)) + geom_point(size=1) +
  theme_minimal() + theme(title = element_text(face="bold", size=12)) +
  labs(title="pearson correlation of signal in rep1 vs. rep2 (2k bins)") + xlab("")

ggsave("reps_pearson_boxplot_2k.pdf", width=10)
ggsave("reps_pearson_boxplot_2k.png", width=10)

##############

