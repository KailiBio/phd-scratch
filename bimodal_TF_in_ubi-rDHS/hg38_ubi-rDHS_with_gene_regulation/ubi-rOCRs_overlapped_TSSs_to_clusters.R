
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
# scatter plot: TSS in gene vs. ubi-TSS in gene (percentage)
####################
coding_percentage = read.table("GRCh38_coding-gene_ubi-rOCR_TSS_percentage_with_cluster.txt")
colnames(coding_percentage) = c("gene", "TSS", "ubi_TSS", "percentage","cluster","ubi_cluster")

ggplot(coding_percentage, aes(x=TSS, y=ubi_TSS, col="#f8766d")) + 
  geom_point(alpha = 0.1) + theme_minimal() + xlim(0,70) + ylim(0,70) +
  theme(legend.position = "none") + 
  geom_abline(intercept = 0, slope = 1, col="grey", linetype="dashed") +
  geom_abline(intercept = 0, slope = 0.5, col="grey", linetype="dotdash")
ggsave("num_TSS_in_coding-gene_count_scatter.pdf")
ggsave("num_TSS_in_coding-gene_count_scatter.png")


ggplot(coding_percentage, aes(x=TSS, y=percentage, group=TSS)) + 
  geom_boxplot(width=0.7, outlier.size=0.3, fill="#f8766d", color="#bdbdbd") +
  theme_minimal() +
  ylab("percentage of ubi-rOCRs overlapped TSSs") + xlab("number of TSSs in a gene") +
  labs(title="percentage of ubi-rOCRs overlapped TSSs\nin protein-coding genes") +
  theme(title = element_text(face="bold",size=12))
ggsave("percentatge_ubi-rOCRs_TSS_in_coding-gene_boxplot.pdf", width = 10, height=5)
ggsave("percentatge_ubi-rOCRs_TSS_in_coding-gene_boxplot.png", width = 10, height=5)

a=cbind(as.numeric(names(table(coding_percentage$TSS))),as.vector(table(coding_percentage$TSS)))
gene_TSS_num = cbind(1:73,rep(0,73))
gene_TSS_num[a[,1],2] = a[,2]
gene_TSS_num = data.frame(gene_TSS_num)
colnames(gene_TSS_num) = c("num","count")
ggplot(gene_TSS_num, aes(x=num, y=count)) + geom_bar(stat="identity", fill="#f8766d") +
  theme_minimal() + theme(title = element_text(face="bold", size=12)) +
  xlab("number of TSSs in a gene") + ylab("number of genes") +
  labs(title="number of gene with different TSSs")

ggsave("number_of_gene_with_TSS_count.pdf", width = 10, height=5)
ggsave("number_of_gene_with_TSS_count.png", width = 10, height=5)


b=cbind(as.numeric(names(table(coding_percentage$cluster))),
        as.vector(table(coding_percentage$cluster)))
gene_cluster_num = cbind(1:45,rep(0,45))
gene_cluster_num[b[,1],2] = as.vector(b[,2])
gene_cluster_num = data.frame(gene_TSS_num)
colnames(gene_cluster_num) = c("num","count")
ggplot(gene_cluster_num, aes(x=num, y=count)) + geom_bar(stat="identity", fill="#f8766d") +
  theme_minimal() + theme(title = element_text(face="bold", size=12)) +
  xlab("number of merged-TSSs in a gene") + ylab("number of genes") +
  labs(title="gene with different number of merged-TSSs")

ggsave("number_of_gene_with_merged-TSS_count.pdf", width = 10, height=5)
ggsave("number_of_gene_with_merged-TSS_count.png", width = 10, height=5)


ggplot(coding_percentage, aes(x=cluster, y=ubi_cluster, group=TSS)) + 
  geom_boxplot(width=0.7, outlier.size=0.3, fill="#f8766d", color="#bdbdbd") +
  theme_minimal() +
  ylab("number of merged-TSSs overlapping ubi-rOCRs") + 
  xlab("number of merged-TSSs in a gene") +
  labs(title="number of merged-TSSs that overlapping ubi-rOCRs\nin protein-coding genes") +
  theme(title = element_text(face="bold",size=12))

ggsave("num_merged-TSS_in_coding-gene_count_boxplot.pdf", width = 10, height=5)
ggsave("num_merged-TSS_in_coding-gene_count_boxplot.png", width = 10, height=5)

##--------------


gene_percentage = read.table("GRCh38_gene_ubi-rOCR_TSS_percentage.txt")
colnames(gene_percentage) = c("gene", "TSS", "ubi_TSS", "percentage")

ggplot(gene_percentage, aes(x=TSS, y=percentage, group=TSS)) + 
  geom_boxplot(width=0.5, fill="#f8766d", color="#bdbdbd", outlier.size=0.3) +
  theme_minimal() +
  ylab("percentage of ubi-rOCRs overlapped TSSs") + xlab("number of TSSs") +
  labs(title="percentage of ubi-rOCRs overlapped TSSs\nin genes") +
  theme(title = element_text(face="bold",size=12))
ggsave("percentatge_ubi-rOCRs_TSS_in_gene_boxplot.pdf")
ggsave("percentatge_ubi-rOCRs_TSS_in_gene_boxplot.png")

##---------------

coding_percentage2 = read.table("GRCh38_coding-gene_ubi-rOCR_TSS_percentage_sort_count.txt")
colnames(coding_percentage2) = c("TSS", "ubi_TSS", "count")
#coding_percentage3 = coding_percentage2[-c(1),]
ggplot(coding_percentage2, aes(x=TSS, y=ubi_TSS, size=count, col=count)) +
  geom_point(aes(size=count), alpha=0.5) + 
  scale_color_gradient(low="#bdbdbd",high="#000000") +
  theme_minimal() + xlim(0,70) + ylim(0,70) + 
  geom_abline(intercept = 0, slope = 1, col="grey", linetype="dashed") +
  geom_abline(intercept = 0, slope = 0.5, col="grey", linetype="dotdash") 
ggsave("num_TSS_in_coding-gene_count_scatter2.pdf")
ggsave("num_TSS_in_coding-gene_count_scatter2.png")
  
####################
# TSS-clusters
####################

ggplot(coding_percentage, aes(x=TSS, y=cluster, group=TSS)) + 
  geom_boxplot(width=0.7, outlier.size=0.3, fill="#f8766d", color="#bdbdbd") +
  theme_minimal() +
  ylab("number of ubi-rOCRs overlapped TSS-clusters") + 
  xlab("number of TSSs in a gene") +
  labs(title="number of ubi-rOCRs overlapped TSS-clusters\nin protein-coding genes") +
  theme(title = element_text(face="bold",size=12))
ggsave("percentatge_ubi-rOCRs_TSS-cluster_in_coding-gene_boxplot.pdf", width = 10, height=5)
ggsave("percentatge_ubi-rOCRs_TSS-cluster_in_coding-gene_boxplot.png", width = 10, height=5)

###----------------------------
### ubi-rOCTRs overlapped genes, 
overlapped_TSS_type = read.table("GRCh38_coding-gene_numOfTSS_overlappedTSSLength.txt")
colnames(overlapped_TSS_type) = c("gene", "tss", "num_of_TSS", "type") 

ggplot(overlapped_TSS_type, aes(x=type, y=num_of_TSS, fill=type)) + geom_violin(width=0.5) +
  geom_boxplot(width=0.1, outlier.size = 0) +
  theme_minimal() + theme(title = element_text(face="bold", size=12)) +
  labs(title="number of TSSs of gene VS. ubi-rOCRs overlapped TSSs type",
       subtitle="wilcox.test: p-value=3.3e-134") +
  ylab("number of TSS in a gene")

ggsave("numOfTSS_overlappedTSStype_violin.pdf")
ggsave("numOfTSS_overlappedTSStype_violin.png")

x = overlapped_TSS_type[overlapped_TSS_type$type=="cluster",]$num_of_TSS
y = overlapped_TSS_type[overlapped_TSS_type$type=="singular",]$num_of_TSS
wilcox.test(x,y)$p.value

### protein-coding gene with only 1 TSS-cluster overlapped
overlapped_1TSS_type = read.table("GRCh38_coding-gene_numOfTSS1_overlappedTSSLength.txt")
colnames(overlapped_1TSS_type) = c("gene", "tss", "num_of_TSS", "type") 

ggplot(overlapped_1TSS_type, aes(x=type, y=num_of_TSS, fill=type)) + geom_violin(width=0.5) +
  geom_boxplot(width=0.1, outlier.size = 0) +
  theme_minimal() + theme(title = element_text(face="bold", size=12)) +
  labs(title="number of TSSs of gene VS. ubi-rOCRs overlapped TSSs type",
       subtitle="wilcox.test: p-value=4.3e-201") +
  ylab("number of TSS in a gene")

ggsave("numOfTSS1_overlappedTSStype_violin.pdf")
ggsave("numOfTSS1_overlappedTSStype_violin.png")

x = overlapped_1TSS_type[overlapped_1TSS_type$type=="cluster",]$num_of_TSS
y = overlapped_1TSS_type[overlapped_1TSS_type$type=="singular",]$num_of_TSS
wilcox.test(x,y)$p.value

### protein-coding gene with only 1 TSS-cluster overlapped
overlapped_2TSS_type = read.table("GRCh38_coding-gene_numOfTSS2_overlappedTSSLength_new.txt")
colnames(overlapped_2TSS_type) = c("gene", "num_of_TSS", "type") 

ggplot(overlapped_2TSS_type, aes(x=type, y=num_of_TSS, fill=type)) + geom_violin(width=0.5) +
  geom_boxplot(width=0.1, outlier.size = 0) +
  theme_minimal() + theme(title = element_text(face="bold", size=12)) +
  labs(title="number of TSSs of gene VS. ubi-rOCRs overlapped TSSs type",
       subtitle="wilcox.test:\ncluster&half p-value=2.8e-18\nhalf&singular p-value=8.5e-50") +
  ylab("number of TSS in a gene")

ggsave("numOfTSS2_overlappedTSStype_violin.pdf")
ggsave("numOfTSS2_overlappedTSStype_violin.png")

x = overlapped_2TSS_type[overlapped_2TSS_type$type=="cluster",]$num_of_TSS
y = overlapped_2TSS_type[overlapped_2TSS_type$type=="half",]$num_of_TSS
z = overlapped_2TSS_type[overlapped_2TSS_type$type=="singular",]$num_of_TSS
wilcox.test(x,y)$p.value
wilcox.test(y,z)$p.value

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

# ggplot(ocr_count, aes(count, fill=type, col=type)) +
#   geom_density(alpha=0.5, n=26) +
#   coord_cartesian(ylim = c(0, 5))

# hist(ocr_count[ocr_count$type=="non_ubi-rOCRs",]$count, breaks=25,freq = F,col="blue")
# hist(ocr_count[ocr_count$type=="ubi-rOCRs",]$count, breaks=26,freq = F, col="red", add=TRUE)

####################
# histogram: TSS nearest TSS distance
####################
distance = read.table("TSS_nearest_TSS_forHist_annotated.bed")
colnames(distance) = c("tss","distance", "type")


ggplot(distance, aes(distance, fill=type, col=type)) + 
  geom_histogram(alpha=0.5, binwidth = 1) +
  scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  scale_color_manual(values = c("#00bfc4","#f8766d")) +
  theme_minimal() +
  geom_vline(xintercept = log10(50), col="grey",linetype="dashed", size=1) +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.8,0.85)) +
  labs(title="distance to neasest TSS for all TSSs") +
  coord_cartesian(xlim=c(0,100))
ggsave("TSSs_distance2nearest_TSS_density.pdf")
ggsave("TSSs_distance2nearest_TSS_density.png")

ggplot(distance, aes(log10(distance), fill=type, col=type)) + 
  geom_density(size=1, alpha=0.5) +
  scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  scale_color_manual(values = c("#00bfc4","#f8766d")) +
  theme_minimal() +
  geom_vline(xintercept = log10(50), col="grey",linetype="dashed", size=1) +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.8,0.85)) +
  labs(title="distance to neasest TSS for all TSSs")
ggsave("TSSs_distance2nearest_TSS_density_log.pdf")
ggsave("TSSs_distance2nearest_TSS_density_log.png")


ggplot(distance, aes(y=log10(distance), x=type, fill=type)) + 
  geom_boxplot(width=0.5) +
  scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  theme_minimal() +
  geom_hline(yintercept = log10(50), col="grey",linetype="dashed", size=1) +
  theme(title = element_text(face="bold",size=12)) +
  labs(title="distance to neasest TSS for all TSSs")

ggsave("TSSs_distance2nearest_TSS_boxplot.pdf")
ggsave("TSSs_distance2nearest_TSS_boxplot.png")

ggplot(distance, aes(y=distance, x=type, fill=type)) + 
  geom_boxplot(width=0.5, outlier.size = 0.1) +
  scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  theme_minimal() + coord_cartesian(ylim = c(0, 150)) +
  geom_hline(yintercept = 50, col="grey",linetype="dashed", size=1) +
  theme(title = element_text(face="bold",size=12)) +
  labs(title="distance to neasest TSS for all TSSs")

ggsave("TSSs_distance2nearest_TSS_boxplot_zoomin.pdf")
ggsave("TSSs_distance2nearest_TSS_boxplot_zoomin.png")

#coord_cartesian(ylim = c(0, 500)) +

####################
# TSS nearest TSS distance (same gene)
####################
distance_gene = read.table("TSS_nearest_TSS_distance_sameGene_annotated.bed")
colnames(distance_gene) = c("gene","tss","distance", "type")


a=distance_gene[distance_gene$distance<200,]
ggplot(a, aes(distance, fill=type, col=type)) + 
  geom_density(size=1, alpha=0.5) +
  scale_fill_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  scale_color_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  theme_minimal() +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.8,0.85)) +
  labs(title="distance to neasest TSS for all TSSs\n(same gene, only distance<200)")

ggsave("TSSs_distance2nearest_TSS_density_sameGene_200.pdf")
ggsave("TSSs_distance2nearest_TSS_density_sameGene_200.png")


ggplot(distance_gene, aes(log10(distance), fill=type, col=type)) + 
  geom_density(size=1, alpha=0.5) +
  scale_fill_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  scale_color_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  theme_minimal() +
  geom_vline(xintercept = log10(2), col="yellow",linetype="dashed", size=1) +
  geom_vline(xintercept = log10(50), col="grey",linetype="dashed", size=1) +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.8,0.85)) +
  labs(title="distance to neasest TSS for all TSSs\n(same gene)")
ggsave("TSSs_distance2nearest_TSS_density_sameGene.pdf")
ggsave("TSSs_distance2nearest_TSS_density_sameGene.png")

ggplot(distance_gene, aes(y=log10(distance), x=type, fill=type)) + 
  geom_boxplot(width=0.5) +
  scale_fill_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  theme_minimal() +
  geom_hline(yintercept = log10(50), col="grey",linetype="dashed", size=1) +
  theme(title = element_text(face="bold",size=12)) +
  labs(title="distance to neasest TSS for all TSSs\n(same gene)")

ggsave("TSSs_distance2nearest_TSS_boxplot_sameGene.pdf")
ggsave("TSSs_distance2nearest_TSS_boxplot_sameGene.png")

ggplot(distance_gene, aes(y=distance, x=type, fill=type)) + 
  geom_boxplot(width=0.5, outlier.size = 0.1) +
  scale_fill_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  theme_minimal() + coord_cartesian(ylim = c(0, 150)) +
  geom_hline(yintercept = 50, col="grey",linetype="dashed", size=1) +
  theme(title = element_text(face="bold",size=12)) +
  labs(title="distance to neasest TSS for all TSSs\n(same gene)")

ggsave("TSSs_distance2nearest_TSS_boxplot_sameGene_zoomin.pdf")
ggsave("TSSs_distance2nearest_TSS_boxplot_sameGene_zoomin.png")



x = distance_gene[distance_gene$type=="TSSs_overlapping_ubi-rOCRs",]$distance
y = distance_gene[distance_gene$type=="TSSs_overlapping_rOCRs",]$distance
z = distance_gene[distance_gene$type=="TSSs_not_overlapping_rOCRs",]$distance
wilcox.test(x,y)$p.value
wilcox.test(x,z)$p.value
wilcox.test(y,z)$p.value
length(x)
length(y)
length(z)
####################
# basic, distance, same gene
####################
basic_distance = read.table("hg38_v28_basic_TSS_multiple_distance_sameGene.bed")
colnames(basic_distance) = c("gene","tss","distance", "type1", "type2")

ggplot(basic_distance, aes(log10(distance), fill=type1, col=type1)) + 
  geom_density(size=1, alpha=0.5) +
  scale_fill_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  scale_color_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  theme_minimal() +
  geom_vline(xintercept = log10(50), col="grey",linetype="dashed", size=1) +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.8,0.85)) +
  labs(title="distance to neasest TSS for all TSSs\n(same gene)")
ggsave("basic_TSSs_distance2nearest_TSS_density_sameGene.pdf")
ggsave("basic_TSSs_distance2nearest_TSS_density_sameGene.png")


a=basic_distance[basic_distance$distance<700,]
ggplot(a, aes(distance, fill=type1, col=type1)) + 
  geom_density(size=1, alpha=0.5) +
  scale_fill_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  scale_color_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  theme_minimal() +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.8,0.85)) +
  labs(title="distance to neasest TSS for all TSSs\n(same gene, only distance<700)")

ggsave("basic_TSSs_distance2nearest_TSS_density_sameGene_700.pdf")
ggsave("basic_TSSs_distance2nearest_TSS_density_sameGene_700.png")


ggplot(basic_distance, aes(y=distance, x=type1, fill=type1)) + 
  geom_boxplot(width=0.5, outlier.size = 0.1) +
  scale_fill_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  theme_minimal() + coord_cartesian(ylim = c(0, 150)) +
  geom_hline(yintercept = 50, col="grey",linetype="dashed", size=1) +
  theme(title = element_text(face="bold",size=12)) +
  labs(title="distance to neasest TSS for all TSSs\n(same gene)")

ggsave("basic_TSSs_distance2nearest_TSS_boxplot_sameGene_zoomin.pdf")
ggsave("basic_TSSs_distance2nearest_TSS_boxplot_sameGene_zoomin.png")

x = basic_distance[basic_distance$type1=="overlap_with_ubi-rOCRs",]$distance
y = basic_distance[basic_distance$type1=="overlap_with_not-ubi_active-rOCRs",]$distance
z = basic_distance[basic_distance$type1=="no_overlap",]$distance
wilcox.test(x,y)$p.value
wilcox.test(x,z)$p.value
wilcox.test(y,z)$p.value
length(x)
length(y)
length(z)
####################
# TSS number in gene: ubi-rOCRs vs. rOCRs
####################
tss_count = read.table("./basic_v28/hg38_v28_basic_gene_labled_TSScount.bed")
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

ggsave("boxplot_number_of_TSSs.pdf")

####################