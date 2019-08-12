
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
# Jun08
# fisher for each sample (n=16)
####################
dat = read.table("fisher_matrix_each_sample.txt", row.names = 1)
dat2 = read.table("fisher_matrix_each_sample_2.txt", row.names = 1)
dat3 = read.table("fisher_matrix_each_sample_3.txt", row.names = 1)
p = apply(dat,1,function(x) fisher.test(cbind(c(x[1],x[2]),c(x[3],x[4])))$p.value)
p_matrix = data.frame(rownames(dat), p)
colnames(p_matrix) = c("sample", "p")
ggplot(p_matrix, aes(x="",y=-log(p,10))) + geom_boxplot(binwidth=0.7) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
  theme_classic() + theme(title=element_text(face="bold", size=12)) +
  geom_hline(yintercept = 2, linetype="dashed", col="red") + ylab("Fisher.test: -log10(p-value)") +
  labs(title="RAMPAGE peaks overlapping\nubi-rOCRs are enrcihed in\nbroad shape peaks (16 samples)") +
  xlab("")
ggsave("promoter_peak_shapes_fisher_eachsample.pdf", width=4, height=7)
ggsave("promoter_peak_shapes_fisher_eachsample.png", width=4, height=7)
####################