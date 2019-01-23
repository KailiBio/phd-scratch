
# -- Kaili
# This script is for making figures for 14_ubi-rOCRs_overlapped_TSSs_to_clusters.sh.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/fig3/")
library(ggplot2)

gene_TSS = read.table("GRCh38_gene_allTSS_count.txt")
gene_ubi_TSS = read.table("GRCh38_gene_ubi-rOCR_TSS_count.txt")

gene_TSS = transform(gene_TSS, type="TSSs")
gene_ubi_TSS = transform(gene_ubi_TSS, type="ubi-rOCRs_overlapped_TSSs")
gene_matrix = rbind(gene_TSS, gene_ubi_TSS)
colnames(gene_matrix) = c("gene","count", "type")

coding_gene_TSS = read.table("GRCh38_coding-gene_allTSS_count.txt")
coding_gene_ubi_TSS = read.table("GRCh38_coding-gene_ubi-rOCR_TSS_count.txt")

coding_gene_TSS = transform(coding_gene_TSS, type="TSSs")
coding_gene_ubi_TSS = transform(coding_gene_ubi_TSS, type="ubi-rOCRs_overlapped_TSSs")
coding_gene_matrix = rbind(coding_gene_TSS,coding_gene_ubi_TSS)
colnames(coding_gene_matrix) = c("coding_gene","count", "type")

coding_percentage = read.table("GRCh38_coding-gene_ubi-rOCR_TSS_percentage.txt")
colnames(coding_percentage) = c("gene", "TSS", "ubi_TSS", "percentage")
####################
# num of TSS in gene
####################

## histogram
ggplot(gene_matrix,aes(count, fill=type)) + 
  geom_histogram(binwidth = 1, position="identity") +
  xlab("number of TSSs/ubi-rOCRs_overlapped_TSSs\nin each gene") +
  theme_minimal() + scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  guides(colour = guide_legend(title="TSSs",override.aes = list(shape = 19)))+
  theme(axis.title.x = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.position = c(0.7,0.8))
ggsave("num_TSS_in_gene_count_histogram.pdf")
ggsave("num_TSS_in_gene_count_histogram.png")

## density
ggplot(gene_matrix,aes(count, color=type)) + geom_density(n=127, size=1.5) +
  xlab("number of TSSs/ubi-rOCRs_overlapped_TSSs\nin each gene") +
  theme_minimal() + scale_color_manual(values = c("#00bfc4","#f8766d")) +  
  guides(colour = guide_legend(title="TSSs",override.aes = list(shape = 19))) +
  theme(axis.title.x = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.position = c(0.7,0.8))
ggsave("num_TSS_in_gene_count_density.pdf")
ggsave("num_TSS_in_gene_count_density.png")
####################
# num of TSS in coding-gene
####################
## histogram
ggplot(coding_gene_matrix,aes(count, fill=type)) + 
  geom_histogram(binwidth = 1, position="identity") +
  xlab("number of TSSs/ubi-rOCRs_overlapped_TSSs\nin each protein-coding gene") +
  theme_minimal() + scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  guides(colour = guide_legend(title="TSSs",override.aes = list(shape = 19)))+
  theme(axis.title.x = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.position = c(0.7,0.8))
ggsave("num_TSS_in_coding-gene_count_histogram.pdf")
ggsave("num_TSS_in_coding-gene_count_histogram.png")

## density
ggplot(coding_gene_matrix,aes(count, color=type)) + geom_density(n=73, size=1.5) +
  xlab("number of TSSs/ubi-rOCRs_overlapped_TSSs\nin each protein-coding gene") +
  theme_minimal() + scale_color_manual(values = c("#00bfc4","#f8766d")) + 
  guides(colour = guide_legend(title="TSSs",override.aes = list(shape = 19)))+
  theme(axis.title.x = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.position = c(0.7,0.8))
ggsave("num_TSS_in_coding-gene_count_density.pdf")
ggsave("num_TSS_in_coding-gene_count_density.png")


####################
# scatter plot: TSS in gene vs. ubi-TSS in gene
####################

ggplot(coding_percentage, aes(x=TSS, y=ubi_TSS, col="#f8766d")) + 
  geom_point(alpha = 0.1) + theme_minimal() + xlim(0,70) + ylim(0,70) +
  theme(legend.position = "none") + 
  geom_abline(intercept = 0, slope = 1, col="grey", linetype="dashed") +
  geom_abline(intercept = 0, slope = 0.5, col="grey", linetype="dotdash")
ggsave("num_TSS_in_coding-gene_count_scatter.pdf")
ggsave("num_TSS_in_coding-gene_count_scatter.png")

####################
# boxplot:TSS nearby TSS numbers
####################
TSS_count = read.table("TSS_neary_TSS_count_annotated.txt")
colnames(TSS_count) = c("chr","s","e","id","dot","strand","gene","count","type")

ggplot(TSS_count, aes(x=type, y=count, fill=type)) + 
  geom_boxplot(width=0.2, position=position_dodge(1)) +
  scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  ylab("number of TSSs within 50bps") +
  labs(title="number of TSSs within TSS/\nubi-rOCRs_overlapped_TSS 50bp region") +
  theme_minimal() +
  theme(axis.title.x=element_blank(), 
        axis.title.y = element_text(face="bold", size=12),
        title = element_text(face="bold", size=12),
        legend.position = "none")
ggsave("TSS_nearby_TSSnum_boxplot.pdf")
ggsave("TSS_nearby_TSSnum_boxplot.png")
####################
# boxplot:TSS nearby ubi-TSS numbers
####################
ubi_TSS_count = read.table("TSS_neary_ubi-TSS_count_annotated.txt")
colnames(ubi_TSS_count) = c("chr","s","e","id","dot","strand","gene","count","type")

ggplot(ubi_TSS_count, aes(x=type, y=count, fill=type)) + geom_boxplot(width=0.3) +
  scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  labs(title = "number of ubi-rOCRs_overlapped_TSSs\nin ±50bp of each TSS") +
  ylab("number of ubi-rOCRs overlapped TSSs in ±50bp region") +
  theme_minimal() +
  theme(axis.title.x=element_blank(), 
        axis.title.y = element_text(face="bold", size=12),
        title = element_text(face="bold", size=12),
        legend.position = "none")
ggsave("TSS_nearby_ubi-TSSnum_boxplot.pdf")
ggsave("TSS_nearby_ubi-TSSnum_boxplot.png")
####################
# scatter: TSS nearby TSS vs.ubi-TSS
####################
count_matrix = read.table("TSS_neary_TSS_ubi-TSS_count_merge.txt")
colnames(count_matrix) = c("id", "type", "tss_count", "ubi_tss_count")

ggplot(count_matrix, aes(x=tss_count, y=ubi_tss_count, col=type)) +
  geom_point(alpha=0.3)

####################
# histogram: number of TSS in rOCR/ubi-rOCRs
####################
ocr_count = read.table("GRCh38_rOCRs_TSS_count_annotated.txt")
colnames(ocr_count) = c("id", "count", "type")

ggplot(ocr_count, aes(count, fill=type)) + 
  geom_histogram(binwidth = 1, aes(y=..density..), position="identity", alpha=0.7) +
  theme_minimal() + scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  guides(colour = guide_legend(title="rOCRs")) +
  theme(legend.position = c(0.7,0.8), 
        legend.title = element_text(face="bold", size=12),
        legend.text = element_text(face="bold", size=12),
        plot.title = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12)) +
  labs(title="number of TSSs in ubi-rOCRs/non-ubi-rOCRs") +
  xlab("number of TSSs in rOCRs")
ggsave("num_TSS_in_rOCRs_density.pdf")
ggsave("num_TSS_in_rOCRs_density.png")

# hist(ocr_count[ocr_count$type=="non_ubi-rOCRs",]$count, breaks=25,freq = F,col="blue")
# hist(ocr_count[ocr_count$type=="ubi-rOCRs",]$count, breaks=26,freq = F, col="red", add=TRUE)
####################
####################


####################

