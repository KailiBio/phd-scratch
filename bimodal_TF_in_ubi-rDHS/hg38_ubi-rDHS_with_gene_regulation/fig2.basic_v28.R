
# -- Kaili
# This script is for making tissue-specificity index of RNA&RAMPAGE.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/fig2/basic_v28/")
library(ggplot2)
library(grid)
library(gridExtra)


#####################
# RNA
#####################
gene_ts = read.table("gene_exp_TSindex_labeled.txt")
colnames(gene_ts) = c("gene", "TSindex", "type")

p1 = ggplot(gene_ts, aes(TSindex, color=type)) +
  geom_density(size=2) + xlab("Tissue-Specificity index") +
  scale_color_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  theme_minimal() + theme(axis.title.x = element_text(size=12, face="bold")) +
  guides(colour = guide_legend(title="genes",override.aes = list(shape = 19))) +
  theme(legend.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.position = c(0.3,0.8))
p1
ggsave("fig2.TSindex_density_1_basic_v28.pdf")


x = gene_ts[gene_ts$type=="genes_whose_TSSs_overlap_ubi-rOCRs",]$TSindex
y = gene_ts[gene_ts$type=="genes_whose_TSSs_overlap_other_active_rOCRs",]$TSindex
wilcox.test(x,y)$p.value


#####################
# RNAPAGE
#####################
tss_ts = read.table("tss_exp_TSindex_labeled.txt")
colnames(tss_ts) = c("gene", "TSindex", "type")
p2 = ggplot(tss_ts, aes(TSindex, color=type)) +
  geom_density(size=1.5) + xlab("Tissue-Specificity index") +
  scale_color_manual(values = c("#6d9df8","#B79F00","#53b64c","#e87d72")) +
  theme_minimal() + theme(axis.title.x = element_text(size=12, face="bold")) +
  guides(colour = guide_legend(title="TSSs",override.aes = list(shape = 19)))+
  theme(legend.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.position = c(0.45,0.8))
p2
ggsave("fig2.TSindex_density_2_basic_v28.pdf")

x = tss_ts[tss_ts$type=="overlap_with_ubi-rOCRs",]$TSindex
y = tss_ts[tss_ts$type=="in_gene_overlapping_ubi-rOCRs",]$TSindex
z = tss_ts[tss_ts$type=="only_overlape_with_other_rOCRs",]$TSindex
wilcox.test(x,y)$p.value
wilcox.test(x,z)$p.value

#####################
p = grid.arrange(p1, p2, nrow=1)
ggsave("fig2.TSindex_density_basic_v28.pdf", p, width = 10)
#####################
