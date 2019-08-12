
# -- Kaili
# This script is for making basic figures for promoter shape analysis.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/fig4/")

library(ggplot2)
library(grid)
library(gridExtra)

###############
###############
data = read.table("promoter_peak_master_list.txt", header = TRUE)

ggplot(data, aes(y=total)) + geom

hist(data$total, breaks=0.5:4783.5)


###############
# GC content
###############
gc_dat = read.table("rOCRs_overlap_hg38_v28_basic_TSS_GCcontent_ID.txt")

ggplot(gc_dat, aes(y = V6, group = V5, fill = V5)) + geom_boxplot(width = 0.5) +
  theme_classic() + 
  theme(title = element_text(face="bold", size=12), legend.title = element_blank(),
        axis.text.x = element_blank()) + ylab("GC content") +
  scale_fill_manual(values = c("#397AF2", "#E73A2F"))
ggsave("rOCRs_overlap_basic_TSS_GCcontent_boxplot.pdf", width=3.5, height=5)
ggsave("rOCRs_overlap_basic_TSS_GCcontent_boxplot.png", width=3.5, height=5)

x = gc_dat[gc_dat$V5=="ubi-rOCRs_overlap_TSS",]$V6
y = gc_dat[gc_dat$V5=="rOCRs_overlap_TSS",]$V6
wilcox.test(x,y)$p.value


###############