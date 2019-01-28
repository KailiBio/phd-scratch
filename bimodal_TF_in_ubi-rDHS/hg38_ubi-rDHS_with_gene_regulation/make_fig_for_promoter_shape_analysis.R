# -- Kaili
# This script is for making figures for 15_promoter_shape.sh.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/fig4/")
library(ggplot2)


####################
# GC content
####################
gc = read.table("TSS_overlapped_rOCRs_GCcontent_ID_annotated.txt")
colnames(gc) = c("chr","s", "e", "id", "GC", "type")

ggplot(gc, aes(x=type, y=GC, fill=type)) + geom_boxplot(width=0.3) +
  theme_minimal() + scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  theme(title=element_text(face="bold",size=12)) +
  labs(title="GC content of all TSS overlapping rOCRs",
       subtitle = "wilcox.txt: p-value=0\n44,267vs8272") +
  ylab("GC content") + xlab("")

ggsave("GC_content_TSSoverlapping_rOCRs_boxplot.pdf")
ggsave("GC_content_TSSoverlapping_rOCRs_boxplot.png")


x = gc[gc$type=="ubi-rOCRs",]$GC
y = gc[gc$type=="remaining_rOCRs",]$GC
wilcox.test(x,y)$p.value

####################


####################
####################

