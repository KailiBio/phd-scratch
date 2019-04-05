
# -- Kaili
# This script is for making TF density figures.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/TFmotif/")

library(ggplot2)
library(gridExtra)

###############
data = read.table("rOCRs_TF_density_withlable.txt")
colnames(data) = c("chr","s","e","id","length","coverage","percentage","count","type")

p1 = ggplot(data, aes(x = type , y=coverage, fill=type)) + 
  geom_violin() + geom_boxplot(width=0.03) +
  theme_minimal() + 
  theme(title = element_text(face="bold",size=12), legend.position = "none",
        axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_manual(values = c("#397AF2", "#E73A2F")) +
  ylab("TF covergae (bp)") + xlab("")
p1

p2 = ggplot(data, aes(x = type, y=percentage, fill=type)) + 
  geom_violin() + geom_boxplot(width=0.03) +
  theme_minimal() + 
  theme(title = element_text(face="bold",size=12), legend.position = "none",
        axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_manual(values = c("#397AF2", "#E73A2F")) +
  ylab("TF density (%)") + xlab("")
p2

p3 = ggplot(data, aes(x = type, y=count, fill=type)) + 
  geom_violin() + geom_boxplot(width=0.03) +
  theme_minimal() + 
  theme(title = element_text(face="bold",size=12), legend.position = "none",
        axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_manual(values = c("#397AF2", "#E73A2F")) +
  ylab("number of TFs") + xlab("")
p3

p = grid.arrange(p1,p2,p3, ncol=3)
p
ggsave("encode_TF_density.pdf", p, width=10)
ggsave("encode_TF_density.png", p, width=10)

#########
a = data[data$type=="ubi-rOCRs",]$coverage
b = data[data$type=="rOCRs",]$coverage
wilcox.test(a,b)$p.value

a = data[data$type=="ubi-rOCRs",]$percentage
b = data[data$type=="rOCRs",]$percentage
wilcox.test(a,b)$p.value

a = data[data$type=="ubi-rOCRs",]$count
b = data[data$type=="rOCRs",]$count
wilcox.test(a,b)$p.value
