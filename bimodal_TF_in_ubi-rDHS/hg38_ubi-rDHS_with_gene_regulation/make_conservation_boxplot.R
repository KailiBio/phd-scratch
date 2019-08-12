
# -- Kaili
# This script is for making conservation boxplot comparing ubi-rOCRs and non-ubi-rOCRs.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/conservation/")

library(ggplot2)
library(gridExtra)

############
# rOCRs - phastCon
############
dat = read.table("hg38_rOCRs_phastCons100way.txt")
colnames(dat) = c("id", "value", "rOCRs")

ggplot(dat, aes(x=rOCRs, y=value, fill = rOCRs)) + geom_boxplot() +
  theme_classic() + theme(title = element_text(face="bold", size=12)) +
  scale_fill_manual(values=c("#397AF2", "#E73A2F")) +
  ylab("100-way phastCons") + xlab("") +
  labs(title="all rOCRs")

ggsave("hg38_rOCRs_phastCons100way.pdf", width = 4, height = 7)
ggsave("hg38_rOCRs_phastCons100way.png", width = 4, height = 7)

a = dat[dat$rOCRs=="ubi-rOCRs",]$value
b = dat[dat$rOCRs=="non_ubi-rOCRs",]$value
wilcox.test(a,b)$p.value

############
# rOCRs - phyloP
############
dat = read.table("hg38_rOCRs_phyloP100way.txt")
colnames(dat) = c("id", "value", "rOCRs")

ggplot(dat, aes(x=rOCRs, y=value, fill = rOCRs)) + geom_boxplot() +
  theme_classic() + theme(title = element_text(face="bold", size=12)) +
  scale_fill_manual(values=c("#397AF2", "#E73A2F")) +
  ylab("100-way phyloP") + xlab("") +
  labs(title="all rOCRs")

ggsave("hg38_rOCRs_phyloP100way.pdf", width = 4, height = 7)
ggsave("hg38_rOCRs_phyloP100way.png", width = 4, height = 7)

a = dat[dat$rOCRs=="ubi-rOCRs",]$value
b = dat[dat$rOCRs=="non_ubi-rOCRs",]$value
wilcox.test(a,b)$p.value

############
# rOCRs TSS - phastCon
############
dat = read.table("hg38_rOCRs_TSS_phastCons100way.txt")
colnames(dat) = c("id", "value", "rOCRs")

a = dat[dat$rOCRs=="ubi-rOCRs",]$value
b = dat[dat$rOCRs=="non_ubi-rOCRs",]$value
p = wilcox.test(a,b)$p.value


ggplot(dat, aes(x=rOCRs, y=value, fill = rOCRs)) + geom_boxplot(outlier.shape = NA) +
  theme_classic() + theme(title = element_text(face="bold", size=12)) +
  scale_fill_manual(values=c("#397AF2", "#E73A2F")) +
  ylab("100-way phastCons") + xlab("") +
  labs(title="rOCRs overlapping TSS")

ggsave("hg38_rOCRs_TSS_phastCons100way.pdf", width = 4, height = 7)
ggsave("hg38_rOCRs_TSS_phastCons100way.png", width = 4, height = 7)


############
# rOCRs TSS - phyloP
############
dat = read.table("hg38_rOCRs_TSS_phyloP100way.txt")
colnames(dat) = c("id", "value", "rOCRs")

ggplot(dat, aes(x=rOCRs, y=value, fill = rOCRs)) + geom_boxplot() +
  theme_classic() + theme(title = element_text(face="bold", size=12)) +
  scale_fill_manual(values=c("#397AF2", "#E73A2F")) +
  ylab("100-way phyloP") + xlab("") +
  labs(title="rOCRs overlapping TSS")

ggsave("hg38_rOCRs_TSS_phyloP100way.pdf", width = 4, height = 7)
ggsave("hg38_rOCRs_TSS_phyloP100way.png", width = 4, height = 7)

a = dat[dat$rOCRs=="ubi-rOCRs",]$value
b = dat[dat$rOCRs=="non_ubi-rOCRs",]$value
wilcox.test(a,b)$p.value

############
# Aug11
############
dat = read.table("hg38_rOCRs_TSS_phastCons100way.txt")
colnames(dat) = c("id", "value", "rOCRs")

a = dat[dat$rOCRs=="ubi-rOCRs",]$value
b = dat[dat$rOCRs=="non_ubi-rOCRs",]$value
p = wilcox.test(a,b)$p.value

m = round(median(a)/median(b),3)


p1 = ggplot(dat, aes(x=rOCRs, y=value, fill = rOCRs)) + geom_boxplot(outlier.shape = NA, width=0.3) +
  theme_classic() + theme(title = element_text(face="bold", size=12), legend.position = "none") +
  scale_fill_manual(values=c("#397AF2", "#E73A2F")) +
  ylab("100-way phastCons") + xlab("") +
  labs(subtitle = paste("wilcox.test:\np-value<2e-16\nmedian FC = ",m,sep=""))
p1

###

dat = read.table("hg38_rOCRs_TSS_phyloP100way.txt")
colnames(dat) = c("id", "value", "rOCRs")

a = dat[dat$rOCRs=="ubi-rOCRs",]$value
b = dat[dat$rOCRs=="non_ubi-rOCRs",]$value
pvalue = format(wilcox.test(a,b)$p.value,scientific = TRUE)

m = round(median(a)/median(b),3)

p2 = ggplot(dat, aes(x=rOCRs, y=value, fill = rOCRs)) + geom_boxplot(outlier.shape = NA, width=0.3) +
  theme_classic() + theme(title = element_text(face="bold", size=12), legend.position = "none") +
  scale_fill_manual(values=c("#397AF2", "#E73A2F")) +
  ylab("100-way phyloP") + xlab("") +
  labs(subtitle = paste("wilcox.test:\np-value=",pvalue,"\nmedian FC = ",m,sep=""))
p2

p = grid.arrange(p1,p2, ncol=2)
p
ggsave("hg38_rOCRs_TSS_conservation.pdf", p, width = 6, height = 7)
ggsave("hg38_rOCRs_TSS_conservation.png", p, width = 6, height = 7)


# dat1 = read.table("hg38_rOCRs_TSS_phastCons100way.txt", row.names = 1)
# colnames(dat1) = c("phastCons", "rOCRs")
# 
# dat2 = read.table("hg38_rOCRs_TSS_phyloP100way.txt", row.names = 1)
# colnames(dat2) = c("phyloP", "rOCRs")
# 
# dat = data.frame(dat1, dat2[rownames(dat1),]$phyloP)
# colnames(dat) = c("phastCons", "rOCRs", "phyloP")
# d = melt(dat)
# 
# ggplot(d, aes(x = variable, y=value, fill = rOCRs)) + geom_boxplot(width = 0.3, outlier.shape = NA) +
#   theme_classic() + theme(title = element_text(face="bold", size=12), legend.position = "none") +
#   scale_fill_manual(values=c("#397AF2", "#E73A2F")) +
#   ylab("100-way phyloP") + xlab("") +
#   labs(subtitle = paste("wilcox.test, p-value=",pvalue,"\nFC=",sep=""))


############


