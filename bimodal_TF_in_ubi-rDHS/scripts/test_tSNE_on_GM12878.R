
# -- Kaili
# This script is for testing tSNE using GM12878 data.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/tSNE/")
library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)


data = read.table("hg19_ubi-rDHS_GM12878_five_dimension_matrix2.txt", header = TRUE, row.names = 1)

################
# histogram
################
jpeg("ubi-rDHS_distance_histogram.jpeg", res = 100)
hist(data$distance, breaks=200, xlab="distance2TSS (log10)", main="10,921 ubi-rDHS",ylim=c(0,100))
dev.off()
#
jpeg("ubi-rDHS_DNase_histogram.jpeg", res = 100)
hist(data$dnase, breaks=200, xlab="dnase (z-score)", main="10,921 ubi-rDHS")
dev.off()
#
jpeg("ubi-rDHS_H3K4me3_histogram.jpeg", res = 100)
hist(data$h3k4me3, breaks=200, xlab="h3k4me3 (z-score)", main="10,921 ubi-rDHS")
dev.off()
#
jpeg("ubi-rDHS_H3K27ac_histogram.jpeg", res = 100)
hist(data$h3k27ac, breaks=200, xlab="h3k27ac (z-score)", main="10,921 ubi-rDHS")
dev.off()
#
jpeg("ubi-rDHS_CTCF_histogram.jpeg", res = 100)
hist(data$ctcf, breaks=200, xlab="ctcf (z-score)", main="10,921 ubi-rDHS")
dev.off()

################
# t-SNE ss
################

set.seed(2)
tSNE <- Rtsne(data, dim = 2,perplexity=30)
plot(tSNE$Y)
info=data.frame(tSNE$Y)

pdf("hg19_ubi-rDHS_tSNE_GM12878_ss.pdf", width = 15, height = 9)
a = ggplot(info, aes(x = X1, y = X2, color = data$type)) + 
  geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by cell-type specific five-group)") + 
  scale_color_manual(values=c("#00b0f0","#06da93","#ffcd00","#e1e1e1","#ff0000"), name = 'cell-type specific five-group',
                     breaks=c("CTCF-only","DNase-only","Enhancer-like_Signatures","Inactive","Promoter-like_Signatures"))
#ggsave("ubi-rDHS_5_dimension_GM12878.jpeg")
#
b = ggplot(info, aes(x = X1, y = X2, color = data$distance)) + 
  geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by distance2TSS)") + 
  scale_color_gradient2(low = "#00b0f0", mid= "#f7f7f7", high = "#ff0000", midpoint = 0, name = "distance2TSS")
#ggsave("ubi-rDHS_5_dimension_GM12878_distance.jpeg")
#
c = ggplot(info, aes(x = X1, y = X2, color = data$dnase)) + 
  geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of DNase)")  + 
  scale_color_gradient2(low = "#00b0f0", mid= "#f7f7f7", high = "#ff0000", midpoint = 0, name = "DNase")
#ggsave("ubi-rDHS_5_dimension_GM12878_dnase.jpeg")
#
d = ggplot(info, aes(x = X1, y = X2, color = data$h3k4me3)) + 
  geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of H3K4me3)")  + 
  scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), values = rescale(c(-10,-5,0,5)), name = "H3K4me3")
#scale_color_gradient2(low = "#00b0f0", mid= "#f7f7f7", high = "#ff0000", midpoint = 0)
#ggsave("ubi-rDHS_5_dimension_GM12878_h3k4me3.jpeg")
#
e = ggplot(info, aes(x = X1, y = X2, color = data$h3k27ac)) + 
  geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of H3K27ac)")  + 
  scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), values = rescale(c(-10,-5,0,5)), name = "H3K27ac")
#ggsave("ubi-rDHS_5_dimension_GM12878_h3k27ac.jpeg")
#
f = ggplot(info, aes(x = X1, y = X2, color = data$ctcf)) + 
  geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of CTCF)")  + 
  scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), values = rescale(c(-10,-5,0,5)), name = "CTCF")
#ggsave("ubi-rDHS_5_dimension_GM12878_ctcf.jpeg")
# scale_color_gradient(low="#2166ac", high="#b2182b")
grid.arrange(a,b,c,d,e,f,nrow=2, ncol=3)
dev.off()


################
# test for more tSNE
################
pdf("hg19_ubi-rDHS_tSNE_GM12878.pdf", width = 15, height = 9)
for(i in 1:10){
  set.seed(i)
  tSNE <- Rtsne(data[,1:5], dim = 2, perplexity=30)
  info=data.frame(tSNE$Y)
  #
  a = ggplot(info, aes(x = X1, y = X2, color = data$type)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (ubi-rDHS, n=10,921)") + 
    scale_color_manual(values=c("#00b0f0","#06da93","#ffcd00","#e1e1e1","#ff0000"), name = 'cell-type specific five-group',
                       breaks=c("CTCF-only","DNase-only","Enhancer-like_Signatures","Inactive","Promoter-like_Signatures"))
  #
  b = ggplot(info, aes(x = X1, y = X2, color = data$distance)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by distance2TSS)") + 
    scale_color_gradient2(low = "#00b0f0", mid= "#f7f7f7", high = "#ff0000", midpoint = 0, name = "distance2TSS")
  #
  c = ggplot(info, aes(x = X1, y = X2, color = data$dnase)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of DNase)")  + 
    scale_color_gradient2(low = "#00b0f0", mid= "#f7f7f7", high = "#ff0000", midpoint = 0, name = "DNase")
  #
  d = ggplot(info, aes(x = X1, y = X2, color = data$h3k4me3)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of H3K4me3)")  + 
    scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), values = rescale(c(-10,-5,0,5)), name = "H3K4me3")
  #
  e = ggplot(info, aes(x = X1, y = X2, color = data$h3k27ac)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of H3K27ac)")  + 
    scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), values = rescale(c(-10,-5,0,5)), name = "H3K27ac")
  #
  f = ggplot(info, aes(x = X1, y = X2, color = data$ctcf)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of CTCF)")  + 
    scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), values = rescale(c(-10,-5,0,5)), name = "CTCF")
  #
  grid.arrange(a, b, c, d, e, f, nrow=2, ncol=3)
}
dev.off()


################
  data2 = read.table("hg19_ubi-rDHS_2_GM12878_five_dimension_matrix2.txt", header = TRUE, row.names = 1)
  
################
# histogram
################
jpeg("ubi-rDHS_2_distance_histogram.jpeg", res = 100)
hist(data2$distance, breaks=200, xlab="distance2TSS (log10)", main="29,773 ubi-rDHS",ylim=c(0,600))
dev.off()
#
jpeg("ubi-rDHS_2_DNase_histogram.jpeg", res = 100)
hist(data2$dnase, breaks=200, xlab="dnase (z-score)", main="29,773 ubi-rDHS")
dev.off()
#
jpeg("ubi-rDHS_2_H3K4me3_histogram.jpeg", res = 100)
hist(data2$h3k4me3, breaks=200, xlab="h3k4me3 (z-score)", main="29,773 ubi-rDHS")
dev.off()
#
jpeg("ubi-rDHS_2_H3K27ac_histogram.jpeg", res = 100)
hist(data2$h3k27ac, breaks=200, xlab="h3k27ac (z-score)", main="29,773 ubi-rDHS")
dev.off()
#
jpeg("ubi-rDHS_2_CTCF_histogram.jpeg", res = 100)
hist(data2$ctcf, breaks=200, xlab="ctcf (z-score)", main="29,773 ubi-rDHS")
dev.off()

################
pdf("hg19_ubi-rDHS_2_tSNE_GM12878.pdf", width = 15, height = 9)
for(i in 1:10){
  set.seed(i)
  tSNE <- Rtsne(data2[,1:5], dim = 2, perplexity=30)
  info=data.frame(tSNE$Y)
  #
  a = ggplot(info, aes(x = X1, y = X2, color = data2$type)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (ubi-rDHS, n=29,733)") + 
    scale_color_manual(values=c("#00b0f0","#06da93","#ffcd00","#e1e1e1","#ff0000"), name = 'cell-type specific five-group',
                       breaks=c("CTCF-only","DNase-only","Enhancer-like_Signatures","Inactive","Promoter-like_Signatures"))
  #
  b = ggplot(info, aes(x = X1, y = X2, color = data2$distance)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by distance2TSS)") + 
    scale_color_gradient2(low = "#00b0f0", mid= "#f7f7f7", high = "#ff0000", midpoint = 0, name = "distance2TSS")
  #
  c = ggplot(info, aes(x = X1, y = X2, color = data2$dnase)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of DNase)")  + 
    scale_color_gradient2(low = "#00b0f0", mid= "#f7f7f7", high = "#ff0000", midpoint = 0, name = "DNase")
  #
  d = ggplot(info, aes(x = X1, y = X2, color = data2$h3k4me3)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of H3K4me3)")  + 
    scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), values = rescale(c(-10,-5,0,5)), name = "H3K4me3")
  #
  e = ggplot(info, aes(x = X1, y = X2, color = data2$h3k27ac)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of H3K27ac)")  + 
    scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), values = rescale(c(-10,-5,0,5)), name = "H3K27ac")
  #
  f = ggplot(info, aes(x = X1, y = X2, color = data2$ctcf)) + 
    geom_point(alpha = 0.3) + theme_gray() + ggtitle("GM12878 (color by z-score of CTCF)")  + 
    scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), values = rescale(c(-10,-5,0,5)), name = "CTCF")
  #
  grid.arrange(a, b, c, d, e, f, nrow=2, ncol=3)
}
dev.off()


