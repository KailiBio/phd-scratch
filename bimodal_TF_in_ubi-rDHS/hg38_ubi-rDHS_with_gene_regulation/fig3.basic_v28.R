
# -- Kaili
# This script is for making ubi-rOCRs TSS cluster figures.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/fig3/basic_v28/")
library(ggplot2)

dat = read.table("nearest_TSS_distance_in_multiple_TSS_genes.txt")
colnames(dat) = c("chr", "s", "e", "tss_id", "score", "strand", "gene_id", "type1", "type2", "distance")
#####################
# histogram of distance
#####################
p1 = ggplot(dat, aes(distance, fill=type1)) + 
  geom_histogram(aes(y=..density..),binwidth = 1, position="identity", alpha=0.5) +
  scale_fill_manual(values = c("#FFCD00","#00B0F0","#FF0000")) +
  theme_minimal() +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.7,0.85), 
        legend.text = element_text(face="bold", size=12)) +
  labs(title="distance to neasest TSS in same gene")
p1
ggsave("fig3.distance_histogram_basic_v28.pdf", p1)
ggsave("fig3.distance_histogram_basic_v28.png", p1)

p2 = p1+ coord_cartesian(xlim=c(0,500))
p2
ggsave("fig3.distance_histogram_zoomin_basic_v28.pdf", p2)


###
p10 = ggplot(dat, aes(log10(distance), fill=type1)) +
  geom_histogram(aes(y=..density..),binwidth = 0.1, position="identity", alpha=0.5) +
  scale_fill_manual(values = c("#B79F00","#00bfc4","#f8766d")) +
  theme_minimal() +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.7,0.85),
        legend.text = element_text(face="bold", size=12)) +
  labs(title="log10(distance) to neasest TSS in same gene\n(binwidth=0.1)") + ylab("log10(distance)")
p10
ggsave("fig3.distance_histogram_log_basic_v28.pdf", p10)
ggsave("fig3.distance_histogram_log_basic_v28.png", p10)


#####################
# density of distance
#####################
p7 = ggplot(dat, aes(distance, col=type1)) + geom_density(bins=100) +
  scale_color_manual(values = c("#FFCD00","#00B0F0","#FF0000")) +
  theme_minimal() +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.7,0.85), 
        legend.text = element_text(face="bold", size=12)) +
  labs(title="distance to neasest TSS in same gene")
p7
ggsave("fig3.distance_densityline_basic_v28.pdf", p7)

p8 = p7 + coord_cartesian(xlim=c(0,10000))
p8
ggsave("fig3.distance_densityline_zoomin_basic_v28.pdf", p8)
ggsave("fig3.distance_densityline_zoomin_basic_v28.png", p8)

hist(dat[dat$type1=="overlap_with_ubi-rOCRs",]$distance, freq = FALSE, breaks=1000)
hist(dat[dat$type1=="overlap_with_not-ubi_active-rOCRs",]$distance, freq = FALSE, breaks=1000)
#####################
# density of distance, log
#####################
p9 = ggplot(dat, aes(log10(distance), col=type1)) + geom_density(bins=1) +
  scale_color_manual(values = c("#FFCD00","#00B0F0","#FF0000")) +
  theme_minimal() +
  theme(title = element_text(face="bold",size=12),
        legend.position = c(0.7,0.85), 
        legend.text = element_text(face="bold", size=12)) +
  labs(title="log0(distance) to neasest TSS in same gene")
p9
ggsave("fig3.distance_densityline_log_basic_v28.pdf", p9)

#####################
# boxplot of distance
#####################
dat$order = c("a")
dat[dat$type1=="overlap_with_not-ubi_active-rOCRs",]$order = "b"
dat[dat$type1=="no_overlap",]$order = "c"
p3 = ggplot(dat, aes(y=distance, x=order, fill=order)) + 
  geom_boxplot(width=0.7, outlier.shape= NA) +
  scale_fill_manual(values = c("#E73A2F","#397AF2","#737373")) +
  theme_classic() + 
  theme(title = element_text(face="bold",size=14),
        axis.text = element_text(face="bold",size=10)) +
  ylab("distance to neasest TSS in same gene") + ylab("") +
  coord_cartesian(ylim=c(0,1300)) 
p3
ggsave("fig3.distance_boxplot_basic_v28_newcolor.pdf", p3, width=4, height=7)

p4 = p3 + coord_cartesian(ylim = c(0, 150)) +
  geom_hline(yintercept = 50, col="grey",linetype="dashed", size=1)
p4
ggsave("fig3.distance_boxplot_zoomin_basic_v28_newcolor.pdf", p4, width=4, height=7)


# ggplot(dat, aes(y=distance, x=order, fill=order)) + 
#   geom_boxplot(width=0.7, outlier.shape= NA) + coord_flip() + coord_flip(ylim=c(0,150)) +
#   theme_classic() + scale_fill_manual(values = c("#E73A2F","#397AF2","#737373")) +
#   theme(title = element_text(face="bold",size=14),
#         axis.text = element_text(face="bold",size=10)) +
#   ylab("distance to neasest TSS in same gene") + ylab("")
# ggsave("fig3.distance_boxplot_zoomin_basic_v28_newcolorH.pdf",width=7, height=4)


x = dat[dat$type1=="overlap_with_ubi-rOCRs",]$distance
y = dat[dat$type1=="overlap_with_not-ubi_active-rOCRs",]$distance
z = dat[dat$type1=="no_overlap",]$distance
wilcox.test(x,y)$p.value
wilcox.test(x,z)$p.value
wilcox.test(y,z)$p.value
length(x)
length(y)
length(z)


#####################
# violin plot of distance
#####################

p5 = ggplot(dat, aes(y=distance, x=type1, fill=type1)) + geom_violin(width=0.5) + 
  scale_fill_manual(values = c("#FFCD00","#00B0F0","#FF0000")) +
  theme_minimal() + 
  theme(title = element_text(face="bold",size=12),
        axis.text.x = element_text(face="bold",size=12)) +
  labs(title="distance to neasest TSS in same gene")
p5
ggsave("fig3.distance_violinplot_basic_v28.pdf", p5)


p6 = ggplot(dat, aes(y=distance, x=type1, fill=type1)) + geom_violin(width=0.5) + 
  scale_fill_manual(values = c("#FFCD00","#00B0F0","#FF0000")) +
  theme_minimal() + 
  theme(title = element_text(face="bold",size=12),
        axis.text.x = element_text(face="bold",size=12)) +
  labs(title="distance to neasest TSS in same gene") +
  coord_cartesian(ylim=c(0,8000))
p6
ggsave("fig3.distance_violinplot_zoomin_basic_v28.pdf", p6)


#####################