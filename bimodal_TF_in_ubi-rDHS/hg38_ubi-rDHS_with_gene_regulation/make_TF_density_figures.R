
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
  x = a/7543
  y = c/26635
  r = x/y
  return(c(p1,p2,p3,x,y,r))
  }

p_vector = data.frame(t(apply(contigency,1,do_fisher)))
row.names(p_vector) = contigency[,1]
colnames(p_vector) = c("two_sided","pre_ubi","pre_non_ubi","ubi_ratio","non_ubi_ratio","FC")

pdf("TF_enrichment_in_rOCRs.pdf")
plot(log10(p_vector$FC), -log10(p_vector$two_sided), pch=20, xlim=c(-1,1), 
     xlab = "log10(ubi/non-ubi)", ylab = "-log10(p-value)", 
     main = "TF motif enrichment in ubi-rOCRs/non-ubi-rOCRs")
dev.off()


# ggplot(p_vector, aes(x=log10(FC), y=-log10(two_sided))) + geom_point() +
#   coord_cartesian(xlim = c(-1,1), ylim = c(0,300)) +
#   theme_minimal() + xlab("log10(ubi/non-ubi)") + ylab("-log10(p-value)")



#####################