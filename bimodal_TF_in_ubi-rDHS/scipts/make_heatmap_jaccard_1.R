
# -- Kaili
# This script is for making heatmap using given data (jaccard).
# one exp per biosample
# locally
# INPUT
# OUTPUT

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/")

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(plyr)

############
# function
############
calculate_jaccard <- function(data){
  n = ncol(data)
  out = matrix(0,n,n)
  for (i in 1:n){
    a=data[,i]
    line = c()
    for (j in 1:n){
      b=data[,j]
      c=count(a+b)
      if(length(c$x)==3){
        j=round( c[c$x==2,]$freq/(c[c$x==1,]$freq+c[c$x==2,]$freq), 2)
      }else{
        j=1
      }
      line <- c(line, j)
    }
    out[,i]=line
  }
  return(out)
}

calculate_jaccard_2 <- function(data){
  n = ncol(data)
  out = matrix(0,n,n)
  for (i in 1:n){
    a=rownames(data[data[,i]==1,])
    line = c()
    for (j in 1:n){
      b=rownames(data[data[,j]==1,])
      c=length(unlist(intersect(a,b)))
      d=length(unlist(union(a,b)))
      jaccard=c/d
      line <- c(line, jaccard)
    }
    out[,i]=line
  }
  return(out)
}

calculate_jaccard_two <- function(data1, data2){
  n = ncol(data1)
  out = matrix(0,n,n)
  for(i in 1:n){
    p = unlist(strsplit(colnames(data1)[i],"_"))[1]
    a = count(data1[,i]+data2[,p])
    a2 = a[a$x==2,]$freq
    a_name = rownames(data1[(data1[,i]+data2[,p])==2,])
    for(j in 1:n){
      q = unlist(strsplit(colnames(data1)[j],"_"))[1]
      b = count(data1[,j]+data2[,q])
      b2 = b[b$x==2,]$freq
      b_name = rownames(data1[(data1[,j]+data2[,q])==2,])
      #
      c=length(intersect(a_name, b_name))
      d=length(union(a_name, b_name))
      out[i,j]=round(c/d,2)
    }
  }
  return(out)
}

calculate_jaccard_two_2 <- function(data1, data2){
  n = ncol(data1)
  out = matrix(0,n,n)
  for(i in 1:n){
    p = colnames(data1)[i]
    a = intersect(rownames(data1[data1[,i]==1,]), rownames(data2[data2[,p]==1,]))
    for(j in 1:n){
      q = colnames(data2)[j]
      b = intersect(rownames(data1[data1[,j]==1,]), rownames(data2[data2[,q]==1,]))
      #
      c=length(intersect(a, b))
      d=length(union(a, b))
      out[i,j]=round(c/d,2)
    }
  }
  return(out)
}

############
# data
############
h3k27ac_0 <- read.table("hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
h3k4me3_0 <- read.table("hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
ctcf_1.64 <- read.table("hg19_ubi-rDHS_CTCF_2_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
ctcf_em <- read.table("hg19_ubi-rDHS_CTCF_2_zscore_em_classification_ccreid.txt", header = TRUE, row.names = 1)

pls <- as.vector(read.table("hg19_ubi-rDHS_PLS_list.txt")[,1])
no_pls <- as.vector(read.table("hg19_ubi-rDHS_non_PLS_list.txt")[,1])

over_h3k4me3 = unlist(intersect(colnames(ctcf_1.64), colnames(h3k4me3_0)))
over_h3k27ac = unlist(intersect(colnames(ctcf_1.64), colnames(h3k27ac_0)))

h3k4me3 = h3k4me3[, over_h3k4me3]
h3k27ac = h3k27ac[, over_h3k27ac]

h3k4me3_pls = h3k4me3[pls,]
h3k27ac_pls = h3k27ac[pls,]
ctcf_1.64_pls = ctcf_1.64[pls,]
ctcf_em_pls = ctcf_em[pls,]

h3k4me3_no_pls = h3k4me3[no_pls,]
h3k27ac_no_pls = h3k27ac[no_pls,]
ctcf_1.64_no_pls = ctcf_1.64[no_pls,]
ctcf_em_no_pls = ctcf_em[no_pls,]

ctcf_1.64_h3k4me3 = ctcf_1.64[, over_h3k4me3]
ctcf_1.64_h3k27ac = ctcf_1.64[, over_h3k27ac]
ctcf_em_h3k4me3 = ctcf_em[, over_h3k4me3]
ctcf_em_h3k27ac = ctcf_em[, over_h3k27ac]

ctcf_1.64_h3k4me3_pls = ctcf_1.64_h3k4me3[pls,]
ctcf_1.64_h3k27ac_pls = ctcf_1.64_h3k27ac[pls,]
ctcf_em_h3k4me3_pls = ctcf_em_h3k4me3[pls,]
ctcf_em_h3k27ac_pls = ctcf_em_h3k27ac[pls,]

ctcf_1.64_h3k4me3_no_pls = ctcf_1.64_h3k4me3[no_pls,]
ctcf_1.64_h3k27ac_no_pls = ctcf_1.64_h3k27ac[no_pls,]
ctcf_em_h3k4me3_no_pls = ctcf_em_h3k4me3[no_pls,]
ctcf_em_h3k27ac_no_pls = ctcf_em_h3k27ac[no_pls,]


###################
# H3K4me3
###################

h3k4me3_j = calculate_jaccard_2(h3k4me3)
colnames(h3k4me3_j) = rownames(h3k4me3_j) = colnames(h3k4me3)

pdf("./heatmap/ubi_rDHS_H3K4me3_Jaccard.pdf")
pheatmap(h3k4me3_j, display_numbers = F, main = "high H3K4me3 ubi-rDHS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - H3K4me3
###################

h3k4me3_pls_j = calculate_jaccard_2(h3k4me3_pls)
colnames(h3k4me3_pls_j) = rownames(h3k4me3_pls_j) = colnames(h3k4me3_pls)

pdf("./heatmap/ubi_rDHS_PLS_H3K4me3_Jaccard.pdf")
pheatmap(h3k4me3_pls_j, display_numbers = F, main = "high H3K4me3 ubi-rDHS-PLS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - H3K27ac
###################

h3k27ac_pls_j = calculate_jaccard_2(h3k27ac_pls)
colnames(h3k27ac_pls_j) = rownames(h3k27ac_pls_j) = colnames(h3k27ac_pls)

pdf("./heatmap/ubi_rDHS_PLS_H3K27ac_Jaccard.pdf")
pheatmap(h3k27ac_pls_j, display_numbers = F, main = "high H3K27ac ubi-rDHS-PLS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# H3K27ac
###################

h3k27ac_j = calculate_jaccard_2(h3k27ac)
colnames(h3k27ac_j) = rownames(h3k27ac_j) = colnames(h3k27ac)

pdf("./heatmap/ubi_rDHS_H3K27ac_Jaccard.pdf")
pheatmap(h3k27ac_j, display_numbers = F, main = "high H3K27ac ubi-rDHS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - H3K27ac
###################

h3k27ac_non_pls_j = calculate_jaccard_2(h3k27ac_no_pls)
colnames(h3k27ac_non_pls_j) = rownames(h3k27ac_non_pls_j) = colnames(h3k27ac_no_pls)

pdf("./heatmap/ubi_rDHS_nonPLS_H3K27ac_Jaccard.pdf")
pheatmap(h3k27ac_non_pls_j, display_numbers = F, main = "high H3K27ac ubi-rDHS-non-PLS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - H3K27ac
###################

h3k4me3_non_pls_j = calculate_jaccard_2(h3k4me3_no_pls)
colnames(h3k4me3_non_pls_j) = rownames(h3k4me3_non_pls_j) = colnames(h3k4me3_no_pls)

pdf("./heatmap/ubi_rDHS_nonPLS_H3K4me3_Jaccard.pdf")
pheatmap(h3k4me3_non_pls_j, display_numbers = F, main = "high H3K4me3 ubi-rDHS-non-PLS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# CTCF 1.64
###################

ctcf_1.64_j = calculate_jaccard_2(ctcf_1.64)
colnames(ctcf_1.64_j) = rownames(ctcf_1.64_j) = colnames(ctcf_1.64)

pdf("./heatmap/ubi_rDHS_CTCF_zscore_1.64_Jaccard.pdf", width = 12, height = 12)
pheatmap(ctcf_1.64_j, display_numbers = F, main = "high CTCF ubi-rDHS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# CTCF em
###################

ctcf_em_j = calculate_jaccard_2(ctcf_em)
colnames(ctcf_em_j) = rownames(ctcf_em_j) = colnames(ctcf_em)

pdf("./heatmap/ubi_rDHS_CTCF_zscore_em_Jaccard.pdf", width = 12, height = 12)
pheatmap(ctcf_em_j, display_numbers = F, main = "high CTCF ubi-rDHS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - CTCF 1.64
###################

ctcf_1.64_pls_j = calculate_jaccard_2(ctcf_1.64_pls)
colnames(ctcf_1.64_pls_j) = rownames(ctcf_1.64_pls_j) = colnames(ctcf_1.64_pls)

pdf("./heatmap/ubi-rDHS_PLS_CTCF_zscore_1.64_Jaccard.pdf", width=12, height = 12)
pheatmap(ctcf_1.64_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS-PLS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - CTCF em
###################

ctcf_em_pls_j = calculate_jaccard_2(ctcf_em_pls)
colnames(ctcf_em_pls_j) = rownames(ctcf_em_pls_j) = colnames(ctcf_em_pls)

pdf("./heatmap/ubi-rDHS_PLS_CTCF_em_Jaccard.pdf", height = 12, width = 12)
pheatmap(ctcf_em_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS-PLS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - CTCF 1.64
###################

ctcf_1.64_no_pls_j = calculate_jaccard_2(ctcf_1.64_no_pls)
colnames(ctcf_1.64_no_pls_j) = rownames(ctcf_1.64_no_pls_j) = colnames(ctcf_1.64_no_pls)

pdf("./heatmap/ubi-rDHS_nonPLS_CTCF_zscore_1.64_Jaccard.pdf", height = 12, width = 12)
pheatmap(ctcf_1.64_no_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS nonPLS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - CTCF em
###################

ctcf_em_no_pls_j = calculate_jaccard_2(ctcf_em_no_pls)
colnames(ctcf_em_no_pls_j) = rownames(ctcf_em_no_pls_j) = colnames(ctcf_em_no_pls)

pdf("./heatmap/ubi-rDHS_nonPLS_CTCF_em_Jaccard.pdf", height = 12, width = 12)
pheatmap(ctcf_em_no_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS nonPLS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# H3K27ac & CTCF 1.64
###################

ctcf_1.64_h3k27ac_j = calculate_jaccard_two_2(ctcf_1.64_h3k27ac, h3k27ac)
colnames(ctcf_1.64_h3k27ac_j) = rownames(ctcf_1.64_h3k27ac_j) = colnames(ctcf_1.64_h3k27ac)

pdf("./heatmap/ubi_rDHS_H3K27ac_CTCF_zscore_1.64_Jaccard.pdf", width = 12, height = 12)
pheatmap(ctcf_1.64_h3k27ac_j, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# H3K27ac & CTCF em
###################

ctcf_em_h3k27ac_j = calculate_jaccard_two_2(ctcf_em_h3k27ac, h3k27ac)
colnames(ctcf_em_h3k27ac_j) = rownames(ctcf_em_h3k27ac_j) = colnames(ctcf_em_h3k27ac)

pdf("./heatmap/ubi_rDHS_H3K27ac_CTCF_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_em_h3k27ac_j, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - H3K4me3 & CTCF 1.64
###################

ctcf_em_h3k4me3_pls_j = calculate_jaccard_two_2(ctcf_em_h3k4me3_pls, h3k4me3_pls)
colnames(ctcf_em_h3k4me3_pls_j) = rownames(ctcf_em_h3k4me3_pls_j) = colnames(ctcf_em_h3k4me3_pls)

pdf("./heatmap/ubi_rDHS_PLS_H3K4me3_CTCF_zscore_1.64_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_em_h3k4me3_pls_j, display_numbers = F, main = "high H3K4me3&CTCF ubi-rDHS PLS(zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - H3K4me3 & CTCF em
###################

ctcf_em_h3k4me3_pls_j = calculate_jaccard_two_2(ctcf_em_h3k4me3_pls, h3k4me3_pls)
colnames(ctcf_em_h3k4me3_pls_j) = rownames(ctcf_em_h3k4me3_pls_j) = colnames(ctcf_em_h3k4me3_pls)

pdf("./heatmap/ubi_rDHS_PLS_H3K4me3_CTCF_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_em_h3k4me3_pls_j, display_numbers = F, main = "high H3K4me3&CTCF ubi-rDHS PLS(EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - H3K27ac & CTCF 1.64
###################

ctcf_1.64_h3k27ac_no_pls_j = calculate_jaccard_two_2(ctcf_1.64_h3k27ac_no_pls, h3k27ac_no_pls)
colnames(ctcf_1.64_h3k27ac_no_pls_j) = rownames(ctcf_1.64_h3k27ac_no_pls_j) = colnames(ctcf_1.64_h3k27ac_no_pls)

pdf("./heatmap/ubi_rDHS_nonPLS_H3K27ac_CTCF_zscore_1.64_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_1.64_h3k27ac_no_pls_j, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS nonPLS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - H3K27ac & CTCF em
###################

ctcf_em_h3k27ac_no_pls_j = calculate_jaccard_two_2(ctcf_em_h3k27ac_no_pls, h3k27ac_no_pls)
colnames(ctcf_em_h3k27ac_no_pls_j) = rownames(ctcf_em_h3k27ac_no_pls_j) = colnames(ctcf_em_h3k27ac_no_pls)

pdf("./heatmap/ubi_rDHS_nonPLS_H3K27ac_CTCF_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_em_h3k27ac_no_pls_j, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS nonPLS(EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
