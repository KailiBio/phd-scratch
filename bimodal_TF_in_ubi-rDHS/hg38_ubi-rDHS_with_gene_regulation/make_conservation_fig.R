
# -- Kaili
# This script is for 

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/conservation/")

library(ggplot2)
library(gridExtra)

##############
# rOCRs
##############
phastcons = read.table("hg38_rOCRs_phastCons7.txt")
p1 = ggplot(phastcons, aes(y = V3, group = V2, fill = V2)) + geom_boxplot(width = 0.1) +
  theme_minimal() + 
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank(), 
        legend.title = element_blank(), legend.position = "none") +
  ylab("phastCons score") + scale_fill_manual(values = c("#00bfc4", "#f8766d"))
  
phylop = read.table("hg38_rOCRs_phyloP7.txt")
p2 = ggplot(phylop, aes(y = V3, group = V2, fill = V2)) + geom_boxplot(width = 0.1) +
  theme_minimal() + 
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank(), legend.title = element_blank()) +
  ylab("phyloP score") + scale_fill_manual(values = c("#00bfc4", "#f8766d"))

p = grid.arrange(p1, p2, ncol = 2)
p
ggsave("conservation_rOCRs_boxplot.pdf", p )
ggsave("conservation_rOCRs_boxplot.png", p )

x1 = phastcons[phastcons$V2=="ubi-rOCRs",]$V3
y1 = phastcons[phastcons$V2=="rOCRs",]$V3
wilcox.test(x1, y1)$p.value

x2 = phylop[phylop$V2=="ubi-rOCRs",]$V3
y2 = phylop[phylop$V2=="rOCRs",]$V3
wilcox.test(x2, y2)$p.value
##############
# rOCRs overlap TSSs
##############
phastcons_TSS = read.table("hg38_rOCRs_overlap_TSS_phastCons7.txt")
p3 = ggplot(phastcons_TSS, aes(y = V3, group = V2, fill = V2)) + geom_boxplot(width = 0.1) +
  theme_minimal() + 
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank(), 
        legend.title = element_blank(), legend.position = "none") +
  ylab("phastCons score") + scale_fill_manual(values = c("#00bfc4", "#f8766d"))

phylop_TSS = read.table("hg38_rOCRs_overlap_TSS_phyloP7.txt")
p4 = ggplot(phylop_TSS, aes(y = V3, group = V2, fill = V2)) + geom_boxplot(width = 0.1) +
  theme_minimal() + 
  theme(title = element_text(face="bold", size=12), axis.text.x = element_blank(), legend.title = element_blank()) +
  ylab("phyloP score") + scale_fill_manual(values = c("#00bfc4", "#f8766d"))

pp = grid.arrange(p3, p4, ncol = 2)
pp
ggsave("conservation_rOCRs_overlap_TSS_boxplot.pdf", pp)
ggsave("conservation_rOCRs_overlap_TSS_boxplot.png", pp)

x3 = phastcons_TSS[phastcons_TSS$V2=="ubi-rOCRs",]$V3
y3 = phastcons_TSS[phastcons_TSS$V2=="rOCRs",]$V3
wilcox.test(x3, y3)$p.value

x4 = phylop_TSS[phylop_TSS$V2=="ubi-rOCRs",]$V3
y4 = phylop_TSS[phylop_TSS$V2=="rOCRs",]$V3
wilcox.test(x4, y4)$p.value
##############
