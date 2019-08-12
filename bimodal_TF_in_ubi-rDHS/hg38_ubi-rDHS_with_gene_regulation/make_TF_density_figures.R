
# -- Kaili
# This script is for making TF density figures.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/TFmotif/")

library(ggplot2)
library(gridExtra)

#####################
# TF motif 
#####################
data = read.table("rOCRs_TF_density_withlable_553.txt")
colnames(data) = c("chr","s","e","id","length","coverage","percentage","count","type")

p1 = ggplot(data, aes(x = type , y=coverage, fill=type)) + 
  geom_violin() + geom_boxplot(width=0.1) +
  theme_minimal() + 
  theme(title = element_text(face="bold",size=12), legend.position = "none",
        axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_manual(values = c("#397AF2", "#E73A2F")) +
  ylab("TF covergae (bp)") + xlab("")
p1

p2 = ggplot(data, aes(x = type, y=percentage, fill=type)) + 
  geom_violin() + geom_boxplot(width=0.1) +
  theme_minimal() + 
  theme(title = element_text(face="bold",size=12), legend.position = "none",
        axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_manual(values = c("#397AF2", "#E73A2F")) +
  ylab("TF density (%)") + xlab("")
p2

p3 = ggplot(data, aes(x = type, y=count, fill=type)) + 
  geom_violin() + geom_boxplot(width=0.1) +
  theme_minimal() + 
  theme(title = element_text(face="bold",size=12), legend.position = "none",
        axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_manual(values = c("#397AF2", "#E73A2F")) +
  ylab("number of TFs") + xlab("")
p3

p = grid.arrange(p1,p2,p3, ncol=3)
p
ggsave("encode_TF_density_553.pdf", p, width=10)
ggsave("encode_TF_density_553.png", p, width=10)

#####################
# TF motif with peak
#####################
data = read.table("rOCRs_TF_density_withlable.txt")
colnames(data) = c("chr","s","e","id","length","coverage","percentage","count","type")

p1 = ggplot(data, aes(x = type , y=coverage, fill=type)) + 
  geom_violin() + geom_boxplot(width=0.03) + theme_classic() + 
  theme(title = element_text(face="bold",size=12), legend.position = "none",
        axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_manual(values = c("#397AF2", "#E73A2F"), name = c("ss1", "ss2")) +
  ylab("TF covergae (bp)") + xlab("") +
  scale_x_discrete(breaks=c("rOCRs","ubi-rOCRs"), labels=c("non-ubi rOCRs", "ubi-rOCRs"))
p1

p2 = ggplot(data, aes(x = type, y=percentage, fill=type)) + 
  geom_violin() + geom_boxplot(width=0.03) + theme_classic() + 
  theme(title = element_text(face="bold",size=12), legend.position = "none",
        axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_manual(values = c("#397AF2", "#E73A2F")) +
  ylab("TF density (%)") + xlab("") +
  scale_x_discrete(breaks=c("rOCRs","ubi-rOCRs"), labels=c("non-ubi rOCRs", "ubi-rOCRs"))
p2

p3 = ggplot(data, aes(x = type, y=count, fill=type)) + 
  geom_violin() + geom_boxplot(width=0.03) + theme_classic() + 
  theme(title = element_text(face="bold",size=12), legend.position = "none",
        axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_manual(values = c("#397AF2", "#E73A2F")) +
  ylab("number of TFs") + xlab("") +
  scale_x_discrete(breaks=c("rOCRs","ubi-rOCRs"),
                   labels=c("non-ubi rOCRs", "ubi-rOCRs"))
p3

p = grid.arrange(p1,p2,p3, ncol=3)
p
ggsave("encode_TF_density.pdf", p, width=12)
ggsave("encode_TF_density.png", p, width=12)

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
#####################
# enrichment
#####################
contigency = read.table("motif_contigency_table.txt", header = TRUE)

do_fisher <- function(line){
  name = as.character(line[1])
  a = as.integer(line[2])
  c = as.integer(line[3])
  ubi = as.integer(line[4])
  non_ubi = as.integer(line[5])
  #
  b = ubi-a
  d = non_ubi-c
  #
  p1 = fisher.test(matrix(c(a,b,c,d),2,2),alternative = "two.sided")$p.value
  p2 = fisher.test(matrix(c(a,b,c,d),2,2),alternative = "greater")$p.value
  p3 = fisher.test(matrix(c(a,b,c,d),2,2),alternative = "less")$p.value
  x = a/ubi
  y = c/non_ubi
  r = x/y
  return(c(p1,p2,p3,x,y,r))
  }

p_vector = data.frame(t(apply(contigency,1,do_fisher)))
row.names(p_vector) = contigency[,1]
colnames(p_vector) = c("two_sided","pre_ubi","pre_non_ubi","ubi_ratio","non_ubi_ratio","FC")
p_vector$q = p.adjust(p_vector$two_sided, method = "fdr")
p_vector$log10_pvalue = -log10(p_vector$two_sided)
p_vector$log10_pvalue = ifelse(p_vector$log10_pvalue==Inf, 350, p_vector$log10_pvalue)

p1 = p_vector[p_vector$two_sided!=0,]
p2 = p1[order(-log10(p1$two_sided),decreasing=TRUE),]

pdf("TF_enrichment_in_rOCRs.pdf")
plot(log2(p_vector$FC), p_vector$log10_pvalue, pch=20, 
     xlab = "log2FC", ylab = "-log10(p-value)", 
     main = "the enrichment of TF motifs in ubi-rOCRs")
text(log2(p2[1,]$FC), y=-log10(p2[1,]$two_sided)+6.5, rownames(p2)[1], col="blue", cex=0.8)
text(log2(p2[2,]$FC)+0.2, y=-log10(p2[2,]$two_sided)+4, rownames(p2)[2], col="blue", cex=0.8)
text(log2(p2[3,]$FC), y=-log10(p2[3,]$two_sided)+6, rownames(p2)[3], col="blue", cex=0.8)
text(log2(p2[4,]$FC), y=-log10(p2[4,]$two_sided)+6, rownames(p2)[4], col="blue", cex=0.8)
text(log2(p2[5,]$FC)-0.35, y=-log10(p2[5,]$two_sided), rownames(p2)[5], col="blue", cex=0.8)
text(log2(p2[6,]$FC), y=-log10(p2[6,]$two_sided)+6, rownames(p2)[6], col="blue", cex=0.8)
text(log2(p2[7,]$FC), y=-log10(p2[7,]$two_sided)+6, rownames(p2)[7], col="blue", cex=0.8)
text(log2(p2[8,]$FC)+0.25, y=-log10(p2[8,]$two_sided)+2, rownames(p2)[8], col="blue", cex=0.8)
text(log2(p2[9,]$FC)-0.25, y=-log10(p2[9,]$two_sided), rownames(p2)[9], col="blue", cex=0.8)
text(log2(p2[10,]$FC), y=-log10(p2[10,]$two_sided)+6, rownames(p2)[10], col="blue", cex=0.8)
text(log2(p2[11,]$FC)+0.25, y=-log10(p2[11,]$two_sided)-2, rownames(p2)[11], col="blue", cex=0.8)
dev.off()

## what's the top TFs?



#####################