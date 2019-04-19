
# -- Kaili
# This script is for making figures for analyzing lincRNA whose TSS overlap ubi-rOCRs.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/lincRNA/")

library(ggplot2)
###############
# distance
###############
dis = read.table("ubi-lincRNA_ubi-PC_TSS_distance.txt")

ggplot(dis, aes(x=V3)) + geom_histogram(binwidth=1)


ggplot(dis, aes(y=log10(V3))) + geom_violin()

summary(dis$V3)
hist(dis$V3, breask = 0.5:243837.5)


boxplot(log10(dis$V3))

#################
# RAMPAGE signal of linc-ubi-rOCRs
#################

RAMPAGE = read.table("sample_RAMPAGE_signal_comparison_linc_in_ubi-rOCRs.txt")
colnames(RAMPAGE) = c("OCR", "signal", "ubi", "cell_line", "sum","PC","type","ubi_type")

ggplot(data=RAMPAGE, aes(x=cell_line, y=log10(signal+0.1), fill=ubi_type)) + 
  geom_boxplot(width=0.5, outlier.shape = NA) + 
  theme_minimal() + ylab("Expression of TSSs\nlog10(TPM+0.1)") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12)) +
  scale_fill_manual(values = c("#e41059","#11c5ee")) +
  coord_cartesian(ylim=c(-1,1.3))

ggsave("ubi-lincRNA_sungroup_RAMPAGE.pdf", width = 12)
ggsave("ubi-lincRNA_subgroup_RAMPAGE.png", width = 12)


#################