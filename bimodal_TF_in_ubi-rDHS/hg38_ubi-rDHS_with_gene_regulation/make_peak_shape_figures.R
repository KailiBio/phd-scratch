
# -- Kaili
# This script is for making peak shape figures.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/peak_shape/")

library(ggplot2)

####################
# peak length distribution
####################
peak_length = read.table("peak_length.txt")
colnames(peak_length) = c("id","score","length")

ggplot(peak_length, aes(x=score, y=length)) + geom_violin(fill = "#00bfc4") + geom_boxplot(width=0.1) +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) +
  labs(title="peak length") + xlab("")
ggsave("peak_length_violin.pdf")
ggsave("peak_length_violin.png")

ggplot(peak_length, aes(x=score, y=length)) + geom_violin(fill = "#00bfc4") + geom_boxplot(width=0.1) +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) +
  labs(title="peak length") + xlab("") + coord_cartesian(ylim=c(0,700))
ggsave("peak_length_violin_zoomin.pdf") 
ggsave("peak_length_violin_zoomin.png")

####################
# sample count
####################
sample_count = read.table("peak_sample_count.txt")
colnames(sample_count) = c("id", "count")

ggplot(sample_count, aes(x=count)) + geom_histogram(bins=155, fill = "#f8766d") +
  theme_classic() + theme(title = element_text(face="bold",size=12)) +
  xlab("number of samples") + labs(title="") +
  geom_vline(xintercept=145, col="grey", linetype = "dashed", size=1)
ggsave("peak_count_histogram.pdf") 
ggsave("peak_count_histogram.png")
  

####################
# venn
####################
library("VennDiagram")
venn.plot <- draw.triple.venn(
  area1 = 6867,
  area2 = 19964,
  area3 = 203420,
  n12 = 5465,
  n23 = 17832,
  n13 = 6690,
  n123 = 5329,
  category = c("ubi-peaks", "ubi-rOCRs", "broad-peaks"),
  fill = c("blue", "red", "green"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "green")
);
grid.draw(venn.plot);
grid.newpage();





####################