
# -- Kaili
# This script is for making tissue-specificity index of RNA&RAMPAGE.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/fig2/")
library(ggplot2)
library(grid)
library(gridExtra)


#####################
# RNA
#####################
gene_ts = read.table("GRCh38_gene_TSindex_quantile.txt")
colnames(gene_ts) = c("gene", "TSindex", "type")

p1 = ggplot(gene_ts, aes(TSindex, color=type)) +
  geom_density(size=2) + xlab("Tissue-Specificity index") +
  scale_color_manual(values = c("#00bfc4","#f8766d")) +
  theme_minimal() + theme(axis.title.x = element_text(size=12, face="bold")) +
  guides(colour = guide_legend(title="genes",override.aes = list(shape = 19))) +
  theme(legend.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.position = c(0.3,0.8))
p1
ggsave("fig2.TSindex_density_1.pdf")
#####################
# RNAPAGE
#####################
## quantile
exp = read.table("hg38_tissue_TSS_exp_matrix_uniqID.txt", header=TRUE, row.names = 1)
library(preprocessCore)
a = normalize.quantiles(as.matrix(exp))
colnames(a) = colnames(exp)
rownames(a) = rownames(exp)
write.table(a,file = "hg38_tissue_TSS_exp_matrix_uniqID_quantile.txt", sep="\t", quote = FALSE)

#######
tss_ts = read.table("GRCh38_TSS_TSindex_quantile.txt")
colnames(tss_ts) = c("gene", "TSindex", "type")
p2 = ggplot(tss_ts, aes(TSindex, color=type)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") +
  scale_color_manual(values = c("#6d9df8","#53b64c","#e87d72")) +
  theme_minimal() + theme(axis.title.x = element_text(size=12, face="bold")) +
  guides(colour = guide_legend(title="TSSs",override.aes = list(shape = 19)))+
  theme(legend.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.position = c(0.45,0.8))
p2
ggsave("fig2.TSindex_density_2.pdf")
#####################
p = grid.arrange(p1, p2, nrow=1)
ggsave("fig2.TSindex_density.pdf", p, width = 10)
#####################
