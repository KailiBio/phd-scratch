
# -- Kaili
# This script is for making figures for validate_CTCFstates_with_CTCFdata.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/ctcfstate_with_ctcfdata/")

library(ggplot2)
library(gridExtra)
library(RColorBrewer)

##################
# signal of lung_0 in state 12
##################
signal_state12 = read.table("lung0_state12_signal.txt")
colnames(signal_state12) = c("id", "signal", "mark")

p1 = ggplot(signal_state12, aes(x = mark, y = signal, group = mark, fill = mark)) +
  geom_violin() + geom_boxplot(width=0.1, outlier.size = 0) +
  xlab("") + theme_minimal() +
  theme(title = element_text(face="bold", size=12)) +
  labs(title="signal of 11 marks in state12 (CTCF state)")
p1

p2 = ggplot(signal_state12, aes(x = mark, y = log10(signal+1e-5), group = mark, fill = mark)) +
  geom_violin() + geom_boxplot(width=0.1, outlier.size = 0) +
  xlab("") + theme_minimal() +
  theme(title = element_text(face="bold", size=12))
p2

p = grid.arrange(p1,p2, ncol=1)
p
ggsave("signal_of_lung0_in_state12_violin.pdf", p, width = 12, height=8)
ggsave("signal_of_lung0_in_state12_violin.png", p, width = 12, height=8)


#####
signal_allstates = read.table("lung0_allstates_signal.txt")
colnames(signal_allstates) = c("id", "signal", "mark")

p1 = ggplot(signal_allstates, aes(x = mark, y = signal, group = mark, fill = mark)) +
  geom_violin() + geom_boxplot(width=0.1, outlier.size = 0) +
  xlab("") + theme_minimal() +
  theme(title = element_text(face="bold", size=12)) +
  labs(title="signal of 11 marks")
p1

p2 = ggplot(signal_allstates, aes(x = mark, y = log10(signal+1e-5), group = mark, fill = mark)) +
  geom_violin() + geom_boxplot(width=0.1, outlier.size = 0) +
  xlab("") + theme_minimal() +
  theme(title = element_text(face="bold", size=12))
p2

p = grid.arrange(p1,p2, ncol=1)
p
ggsave("signal_of_lung0_in_allstates_violin.pdf", p, width = 12, height=8)
ggsave("signal_of_lung0_in_allstates_violin.png", p, width = 12, height=8)


##################
# number of CTCF states in different samples
##################
state_num = read.table("sample_state12_count.txt")
colnames(state_num) = c("sample", "num")
a = t(matrix(unlist(strsplit(as.vector(state_num$sample), split="_")), nrow=2))
state_num = data.frame(cbind(state_num, a))
colnames(state_num) = c("sample", "num", "tissue","timepoint")

ggplot(state_num, aes(x = timepoint, y = num, fill = timepoint)) + 
  geom_bar(stat = "identity") + facet_wrap(~tissue) +
  theme_minimal() +
  theme(title = element_text(face="bold", size=12)) +
  labs(title="number of state12 (CTCF state) in each sample") +
  xlab("") + ylab("number of state12 bins")

ggsave("number_of_state12_in_each_sample_barplot.pdf", width=12)
ggsave("number_of_state12_in_each_sample_barplot.png", width=12)

##################
# regression
##################
signal_matrix = read.table("lung0_signal_random100k.txt", header = TRUE, row.names = 1)

r_squared_matrix = matrix(rep(0,11),ncol=11)
colnames(r_squared_matrix) = colnames(signal_matrix)
rownames(r_squared_matrix) = c("r_square")

model1 = lm(signal_matrix$CTCF~signal_matrix$ATAC + signal_matrix$DNAme + signal_matrix$H3K4me1 +
             signal_matrix$H3K4me2 + signal_matrix$H3K4me3 + signal_matrix$H3K9me3 +
             signal_matrix$H3K9ac + signal_matrix$H3K27me3 + signal_matrix$H3K27ac +
             signal_matrix$H3K36me3)
summary(model1)$adj.r.squared
r_squared_matrix[1,"CTCF"]=summary(model1)$adj.r.squared

# model2 = lm(log(signal_matrix$CTCF+1e-5)~log(signal_matrix$ATAC+1e-5) + 
#               log(signal_matrix$DNAme+1e-5) + log(signal_matrix$H3K4me1+1e-5) +
#               log(signal_matrix$H3K4me2+1e-5) + log(signal_matrix$H3K4me3+1e-5) + 
#               log(signal_matrix$H3K9me3+1e-5) + log(signal_matrix$H3K9ac+1e-5) + 
#               log(signal_matrix$H3K27me3+1e-5) + log(signal_matrix$H3K27ac+1e-5) +
#               log(signal_matrix$H3K36me3+1e-5))
# summary(model2)$adj.r.squared

model3 = lm(signal_matrix$H3K4me3~signal_matrix$ATAC + signal_matrix$DNAme + signal_matrix$H3K4me1 +
              signal_matrix$H3K4me2 + signal_matrix$H3K9me3 +
              signal_matrix$H3K9ac + signal_matrix$H3K27me3 + signal_matrix$H3K27ac +
              signal_matrix$H3K36me3 + signal_matrix$CTCF)
summary(model3)$adj.r.squared
r_squared_matrix[1,"H3K4me3"]=summary(model3)$adj.r.squared

# model4 = lm(log(signal_matrix$H3K4me3+1e-5) ~log(signal_matrix$ATAC+1e-5) + 
#               log(signal_matrix$DNAme+1e-5) + log(signal_matrix$H3K4me1+1e-5) +
#               log(signal_matrix$H3K4me2+1e-5) + 
#               log(signal_matrix$H3K9me3+1e-5) + log(signal_matrix$H3K9ac+1e-5) + 
#               log(signal_matrix$H3K27me3+1e-5) + log(signal_matrix$H3K27ac+1e-5) +
#               log(signal_matrix$H3K36me3+1e-5))
# summary(model4)$adj.r.squared

model5 = lm(signal_matrix$H3K27ac~signal_matrix$ATAC + signal_matrix$DNAme + signal_matrix$H3K4me1 +
              signal_matrix$H3K4me2 + signal_matrix$H3K4me3 + signal_matrix$H3K9me3 +
              signal_matrix$H3K9ac + signal_matrix$H3K27me3 + 
              signal_matrix$H3K36me3 + signal_matrix$CTCF)
summary(model5)$adj.r.squared
r_squared_matrix[1,"H3K27ac"]=summary(model5)$adj.r.squared


model6 = lm(signal_matrix$ATAC~signal_matrix$H3K27ac + signal_matrix$DNAme + signal_matrix$H3K4me1 +
              signal_matrix$H3K4me2 + signal_matrix$H3K4me3 + signal_matrix$H3K9me3 +
              signal_matrix$H3K9ac + signal_matrix$H3K27me3 + 
              signal_matrix$H3K36me3 + signal_matrix$CTCF)
summary(model6)$adj.r.squared
r_squared_matrix[1,"ATAC"]=summary(model6)$adj.r.squared

model6 = lm(signal_matrix$DNAme~signal_matrix$H3K27ac + signal_matrix$ATAC + signal_matrix$H3K4me1 +
              signal_matrix$H3K4me2 + signal_matrix$H3K4me3 + signal_matrix$H3K9me3 +
              signal_matrix$H3K9ac + signal_matrix$H3K27me3 + 
              signal_matrix$H3K36me3 + signal_matrix$CTCF)
summary(model6)$adj.r.squared
r_squared_matrix[1,"DNAme"]=summary(model6)$adj.r.squared

model6 = lm(signal_matrix$H3K4me1~signal_matrix$H3K27ac + signal_matrix$ATAC + signal_matrix$DNAme +
              signal_matrix$H3K4me2 + signal_matrix$H3K4me3 + signal_matrix$H3K9me3 +
              signal_matrix$H3K9ac + signal_matrix$H3K27me3 + 
              signal_matrix$H3K36me3 + signal_matrix$CTCF)
summary(model6)$adj.r.squared
r_squared_matrix[1,"H3K4me1"]=summary(model6)$adj.r.squared

model6 = lm(signal_matrix$H3K4me2~signal_matrix$H3K27ac + signal_matrix$ATAC + signal_matrix$DNAme +
              signal_matrix$H3K4me1 + signal_matrix$H3K4me3 + signal_matrix$H3K9me3 +
              signal_matrix$H3K9ac + signal_matrix$H3K27me3 + 
              signal_matrix$H3K36me3 + signal_matrix$CTCF)
summary(model6)$adj.r.squared
r_squared_matrix[1,"H3K4me2"]=summary(model6)$adj.r.squared

model6 = lm(signal_matrix$H3K9me3~signal_matrix$H3K27ac + signal_matrix$ATAC + signal_matrix$DNAme +
              signal_matrix$H3K4me1 + signal_matrix$H3K4me3 + signal_matrix$H3K4me2 +
              signal_matrix$H3K9ac + signal_matrix$H3K27me3 + 
              signal_matrix$H3K36me3 + signal_matrix$CTCF)
summary(model6)$adj.r.squared
r_squared_matrix[1,"H3K9me3"]=summary(model6)$adj.r.squared


model6 = lm(signal_matrix$H3K9ac~signal_matrix$H3K27ac + signal_matrix$ATAC + signal_matrix$DNAme +
              signal_matrix$H3K4me1 + signal_matrix$H3K4me3 + signal_matrix$H3K4me2 +
              signal_matrix$H3K9me3 + signal_matrix$H3K27me3 + 
              signal_matrix$H3K36me3 + signal_matrix$CTCF)
summary(model6)$adj.r.squared
r_squared_matrix[1,"H3K9ac"]=summary(model6)$adj.r.squared

model6 = lm(signal_matrix$H3K27me3~signal_matrix$H3K27ac + signal_matrix$ATAC + signal_matrix$DNAme +
              signal_matrix$H3K4me1 + signal_matrix$H3K4me3 + signal_matrix$H3K4me2 +
              signal_matrix$H3K9me3 + signal_matrix$H3K9ac + 
              signal_matrix$H3K36me3 + signal_matrix$CTCF)
summary(model6)$adj.r.squared
r_squared_matrix[1,"H3K27me3"]=summary(model6)$adj.r.squared

model6 = lm(signal_matrix$H3K36me3~signal_matrix$H3K27ac + signal_matrix$ATAC + signal_matrix$DNAme +
              signal_matrix$H3K4me1 + signal_matrix$H3K4me3 + signal_matrix$H3K4me2 +
              signal_matrix$H3K9me3 + signal_matrix$H3K9ac + 
              signal_matrix$H3K27me3 + signal_matrix$CTCF)
summary(model6)$adj.r.squared
r_squared_matrix[1,"H3K36me3"]=summary(model6)$adj.r.squared

r = cbind(r_squared_matrix[1,],r_squared_matrix[1,])
library(pheatmap)
pdf("random100k_regression_rsquare.pdf", width = 3, height=7)
pheatmap(r, display_numbers = TRUE, breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
col = c(brewer.pal(9,"YlOrBr"),"white"), fontsize = 12, main = "r-square")
dev.off()

#################
# scatter plot
#################
#pdf("random100k_scatter.pdf")
library("PerformanceAnalytics")
chart.Correlation(signal_matrix, histogram=TRUE, pch=19)
library(psych)
pairs.panels(signal_matrix, scale=TRUE)



d=read.table("Random-Summary.txt")
df <- data.frame(x = signal_matrix$ATAC, y = signal_matrix$CTCF, 
                 signal_matrix = densCols(log(signal_matrix$ATAC+0.1,10), log(signal_matrix$CTCF+0.1,10), colramp = colorRampPalette(rev(rainbow(15, end = 4/6)))))
ggplot(df, aes(x=log(x+0.1,10), y=log(y+0.1,10), color=signal_matrix))+
  geom_point()+scale_color_identity()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border=element_rect(colour = "black", fill=NA), legend.position="none")+
  xlab("Log ATAC Signal")+ylab("Log CTCF Signal")


smoothScatter(log2(signal_matrix$ATAC+0.01), log2(signal_matrix$CTCF+0.01), nrpoints=1)

ggplot(signal_matrix, aes(x=log10(ATAC+0.01), y=log10(CTCF+0.01))) + geom_point(col="grey") + 
  stat_density_2d(aes(fill = stat(level)), geom="polygon") +
  scale_fill_gradient(low="blue", high="red") + 
  theme_minimal() + theme(title=element_text(size=12, face="bold")) 
ggsave("random100k_ATAC_CTCF_scatter.pdf")
ggsave("random100k_ATAC_CTCF_scatter.png")


ggplot(signal_matrix, aes(x=ATAC, y=CTCF)) + stat_bin_hex()


ggplot(signal_matrix, aes(x=H3K4me3, y=H3K9ac)) + geom_point(col="grey", size=0.5) +
  stat_density_2d(aes(fill = ..level..), geom="polygon") +
  scale_fill_gradient(low="blue", high="red") + 
  theme_minimal() + theme(title=element_text(size=12, face="bold")) 
ggsave("random100k_H3K4me3_H3K9ac_scatter.pdf")
ggsave("random100k_H3K4me3_H3K9ac_scatter.png")

ggplot(signal_matrix, aes(x=H3K4me3, y=H3K27ac)) + geom_point(col="grey", size=0.5) +
  stat_density_2d(aes(fill = ..level..), geom="polygon") +
  scale_fill_gradient(low="blue", high="red") + 
  theme_minimal() + theme(title=element_text(size=12, face="bold")) 
ggsave("random100k_H3K4me3_H3K27ac_scatter.pdf")
ggsave("random100k_H3K4me3_H3K27ac_scatter.png")

##################
# correlation
##################
pcor = cor(signal_matrix, method="pearson")
scor = cor(signal_matrix, method="spearman")

library(corrplot)

pdf("random100k_pearson_cor.pdf")
corrplot.mixed(pcor,number.cex=0.8, order = "alphabet", upper.col = brewer.pal(11,"RdYlBu")[11:1],
               lower.col = brewer.pal(11,"RdYlBu")[11:1], tl.cex=0.7, tl.col="black",
               main="\n\npearson correlation\n(100k random bins)")
dev.off()

pdf("random100k_spearman_cor.pdf")
corrplot.mixed(scor,number.cex=0.8, order = "alphabet", upper.col = brewer.pal(11,"RdYlBu")[11:1],
               lower.col = brewer.pal(11,"RdYlBu")[11:1], tl.cex=0.7, tl.col="black",
               main="\n\nspearman correlation\n(100k random bins)")
dev.off()


##################
# all bins
##################
dat = read.table("lung0_signal_allbins.txt", header=TRUE, row.names = 1)

pcor_all = cor(dat, method="pearson")
scor_all = cor(dat, method="spearman")

pdf("allbins_pearson_cor.pdf")
corrplot.mixed(pcor_all,number.cex=0.8, order = "alphabet", upper.col = brewer.pal(11,"RdYlBu")[11:1],
               lower.col = brewer.pal(11,"RdYlBu")[11:1], tl.cex=0.7, tl.col="black",
               main="\n\npearson correlation")
dev.off()

pdf("allbins_spearman_cor.pdf")
corrplot.mixed(scor_all,number.cex=0.8, order = "alphabet", upper.col = brewer.pal(11,"RdYlBu")[11:1],
               lower.col = brewer.pal(11,"RdYlBu")[11:1], tl.cex=0.7, tl.col="black",
               main="\n\nspearman correlation")
dev.off()


model1 = lm(dat$CTCF~dat$ATAC + dat$DNAme + dat$H3K4me1 + dat$H3K4me2 + dat$H3K4me3 + dat$H3K9me3 +
              dat$H3K9ac + dat$H3K27me3 + dat$H3K27ac + dat$H3K36me3)
summary(model1)$adj.r.squared


######## 
#violin
library(reshape2)
library(ggplot2)
dat_melt = melt(dat)
ggplot(dat_melt, aes(x=variable, y =value, group = variable)) + geom_violin() +
  theme_minimal() + theme(title=element_text(size=12, face="bold"))

dat_random = melt(signal_matrix)
ggplot(dat_random, aes(x=variable, y =value, group = variable)) + geom_violin() +
  theme_minimal() + theme(title=element_text(size=12, face="bold"))

##################
# 3bins
##################
dat_3bins = read.table("lung0_signal_random100k_3bins.txt", header = TRUE, row.names = 1)

pcor = cor(dat_3bins, method="pearson")
scor = cor(dat_3bins, method="spearman")

pdf("random100k_3bins_pearson_cor.pdf")
corrplot.mixed(pcor,number.cex=0.8, order = "alphabet", upper.col = brewer.pal(11,"RdYlBu")[11:1],
               lower.col = brewer.pal(11,"RdYlBu")[11:1], tl.cex=0.7, tl.col="black",
               main="\n\npearson correlation\n(random 100k, ±1bin)")
dev.off()

pdf("random100k_3bins_spearman_cor.pdf")
corrplot.mixed(scor,number.cex=0.8, order = "alphabet", upper.col = brewer.pal(11,"RdYlBu")[11:1],
               lower.col = brewer.pal(11,"RdYlBu")[11:1], tl.cex=0.7, tl.col="black",
               main="\n\nspearman correlation\n(random 100k, ±1bin)")
dev.off()

##################
# 5bins
##################
dat_5bins = read.table("lung0_signal_random100k_5bins.txt", header = TRUE, row.names = 1)

pcor = cor(dat_5bins, method="pearson")
scor = cor(dat_5bins, method="spearman")

pdf("random100k_5bins_pearson_cor.pdf")
corrplot.mixed(pcor,number.cex=0.8, order = "alphabet", upper.col = brewer.pal(11,"RdYlBu")[11:1],
               lower.col = brewer.pal(11,"RdYlBu")[11:1], tl.cex=0.7, tl.col="black",
               main="\n\npearson correlation\n(random 100k, ±2bins)")
dev.off()

pdf("random100k_5bins_spearman_cor.pdf")
corrplot.mixed(scor,number.cex=0.8, order = "alphabet", upper.col = brewer.pal(11,"RdYlBu")[11:1],
               lower.col = brewer.pal(11,"RdYlBu")[11:1], tl.cex=0.7, tl.col="black",
               main="\n\nspearman correlation\n(random 100k, ±2bins)")
dev.off()

##################
# 7bins
##################
dat_7bins = read.table("lung0_signal_random100k_7bins.txt", header = TRUE, row.names = 1)

pcor = cor(dat_7bins, method="pearson")
scor = cor(dat_7bins, method="spearman")

pdf("random100k_7bins_pearson_cor.pdf")
corrplot.mixed(pcor,number.cex=0.8, order = "alphabet", upper.col = brewer.pal(11,"RdYlBu")[11:1],
               lower.col = brewer.pal(11,"RdYlBu")[11:1], tl.cex=0.7, tl.col="black",
               main="\n\npearson correlation\n(random 100k, ±3bins)")
dev.off()

pdf("random100k_7bins_spearman_cor.pdf")
corrplot.mixed(scor,number.cex=0.8, order = "alphabet", upper.col = brewer.pal(11,"RdYlBu")[11:1],
               lower.col = brewer.pal(11,"RdYlBu")[11:1], tl.cex=0.7, tl.col="black",
               main="\n\nspearman correlation\n(random 100k, ±3bins)")
dev.off()

##################

##################
