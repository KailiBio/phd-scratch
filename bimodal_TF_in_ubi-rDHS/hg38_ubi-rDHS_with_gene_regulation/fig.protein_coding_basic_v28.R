
# -- Kaili
# This script is for making figures for protein-coding gene/TSSs in v28 basic.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/v28_basic_protein/")
library(ggplot2)
library(grid)
library(gridExtra)


####################
# fig1
####################
cell = c("A172","Daoy","GM23248","GM23338","hepatocyte","HT1080","LHCN-M2","myotube",
         "NCI-H460","neural_progenitor_cell","RPMI-7951","SJCRH30","SJSA1",
         "skeletal_muscle_myoblast","SK-MEL-5","SK-N-DZ")
#
DNase = read.table("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/fig1/basic_v28/sample_DNase_signal_comparison.txt")
colnames(DNase) = c("OCR", "signal", "ubi", "cell_line", "sum")
DNase=transform(DNase, group="a")
DNase$group = as.vector(DNase$group)
DNase[DNase$ubi=="non-ubi_active-rOCR",]$group="b"
p1 = ggplot(data=DNase, aes(x=cell_line, y=signal, fill=group)) + 
  geom_boxplot(width=0.3, outlier.shape = NA) +
  theme_classic() + ylab("signal of DNase I") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        legend.position = "none") + coord_cartesian(ylim=c(1,6))
p1
#
RNA = read.table("sample_RNA_signal_comparison.txt")
colnames(RNA) = c("OCR", "signal", "ubi", "cell_line", "sum")
RNA=transform(RNA, group="a")
RNA$group = as.vector(RNA$group)
RNA[RNA$ubi=="genes_whose_TSSs_overlap_other_active_rOCRs",]$group="b"
p2 = ggplot(data=RNA, aes(x=cell_line, y=log10(signal+0.1), fill=group)) + 
  geom_boxplot(width=0.3, outlier.shape = NA) +
  theme_classic() + ylab("Expression of transcripts\nlog10(TPM+0.1)") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        legend.position = "none") + coord_cartesian(ylim=c(-1,5))
p2
# RAMPAGE
RAMPAGE = read.table("sample_RAMPAGE_signal_comparison.txt")
colnames(RAMPAGE) = c("OCR", "signal", "ubi", "cell_line", "sum")
RAMPAGE=transform(RAMPAGE, group="a")
RAMPAGE$group = as.vector(RAMPAGE$group)
RAMPAGE[RAMPAGE$ubi=="TSS_overlap_non-ubi_active-rOCRs",]$group="b"
p3 = ggplot(data=RAMPAGE, aes(x=cell_line, y=log10(signal+0.1), fill=group)) + 
  geom_boxplot(width=0.3, outlier.shape = NA) + 
  theme_classic() + ylab("Expression of TSSs\nlog10(TPM+0.1)") +
  theme(axis.text.x = element_text(size=12, face="bold", angle=45, vjust=0.6, hjust=0.6), 
        axis.text.y = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x=element_blank(), legend.position = "none") +
  coord_cartesian(ylim=c(-1,3))
p3
# 
p = grid.arrange(p1, p2, p3, ncol=1,
                 layout_matrix = rbind(c(1),c(2),c(3),c(3)))
ggsave("fig1.basic_v28_protein_coding.pdf", p)

for(i in 1:length(cell)){
  a = as.numeric(RNA[RNA$cell_line==cell[i] & RNA$group=="a",]$signal)
  b = as.numeric(RNA[RNA$cell_line==cell[i] & RNA$group=="b",]$signal)
  p = wilcox.test(a,b)$p.value
  print(p)
}

for(i in 1:length(cell)){
  a = as.numeric(RAMPAGE[RAMPAGE$cell_line==cell[i] & RAMPAGE$group=="a",]$signal)
  b = as.numeric(RAMPAGE[RAMPAGE$cell_line==cell[i] & RAMPAGE$group=="b",]$signal)
  p = wilcox.test(a,b)$p.value
  print(p)
}
####################
# fig2
####################
# gene
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
# TSS
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
#
p = grid.arrange(p1, p2, nrow=1)
ggsave("fig2.basic_v28_protein_coding.pdf", p, width = 10)


####################
# fig3
####################
dat = read.table("nearest_TSS_distance_in_multiple_TSS_genes.txt")
colnames(dat) = c("chr", "s", "e", "tss_id", "score", "strand", "gene_id", "type1", "type2", "distance")
#
p1 = ggplot(dat, aes(distance, fill=type1)) + 
  geom_histogram(aes(y=..density..),binwidth = 1, position="identity", alpha=0.5) +
  scale_fill_manual(values = c("#FFCD00","#00B0F0","#FF0000")) +
  theme_minimal() +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.7,0.85), 
        legend.text = element_text(face="bold", size=12)) +
  labs(title="distance to neasest TSS in same gene")
p1
ggsave("fig3.distance_histogram_basic_v28_protein_coding.pdf", p1)
ggsave("fig3.distance_histogram_basic_v28_protein_coding.png", p1)
p2 = p1+ coord_cartesian(xlim=c(0,500))
p2
ggsave("fig3.distance_histogram_zoomin_basic_v28_protein_coding.pdf", p2)
#
p3 = ggplot(dat, aes(y=distance, x=type1, fill=type1)) + 
  geom_boxplot(width=0.5, outlier.size = 0.1) +
  scale_fill_manual(values = c("#FFCD00","#00B0F0","#FF0000")) +
  theme_minimal() + 
  geom_hline(yintercept = 50, col="grey",linetype="dashed", size=1) +
  theme(title = element_text(face="bold",size=12),
        axis.text.x = element_text(face="bold",size=12)) +
  labs(title="distance to neasest TSS in same gene")
p3
ggsave("fig3.distance_boxplot_basic_v28_protein_coding.pdf", p3)

p4 = p3 + coord_cartesian(ylim = c(0, 150))
p4
ggsave("fig3.distance_boxplot_zoomin_basic_v28_protein_coding.pdf", p4)

####
tss_count = read.table("hg38_v28_basic_gene_labled_TSScount_protein_coding.bed")
colnames(tss_count) = c("gene_id", "type", "tss_count")

ggplot(tss_count, aes(x = type, y=tss_count, fill = type)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) + 
  coord_cartesian(ylim = c(0, 8)) +
  theme_minimal() +
  scale_fill_manual(values = c("#FFCD00","#00B0F0","#FF0000")) +
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank(),
        legend.position = c(0.3,0.8)) +
  ylab("number of TSSs") + xlab("gene") +
  labs(title="number of TSSs in gene")

ggsave("fig3.boxplot_number_of_TSSs_basic_v28_protein_coding.pdf")





####################
