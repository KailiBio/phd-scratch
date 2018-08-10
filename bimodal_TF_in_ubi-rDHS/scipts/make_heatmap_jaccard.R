
# -- Kaili
# This script is for making heatmap using given data (jaccard).
# locally
# INPUT
# OUTPUT

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/")

library(pheatmap)
library(RColorBrewer)
library(dplyr)

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

############
# data
############
h3k27ac <- read.table("hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
h3k4me3 <- read.table("hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
ctcf_1.64 <- read.table("hg19_ubi-rDHS_CTCF_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
ctcf_em <- read.table("hg19_ubi-rDHS_CTCF_zscore_classification_ccreid.txt", header = TRUE, row.names = 1)

pls <- as.vector(read.table("hg19_ubi-rDHS_PLS_list.txt")[,1])
no_pls <- as.vector(read.table("hg19_ubi-rDHS_non_PLS_list.txt")[,1])

h3k4me3_pls = h3k4me3[pls,]
h3k27ac_pls = h3k27ac[pls,]
ctcf_1.64_pls = ctcf_1.64[pls,]
ctcf_em_pls = ctcf_em[pls,]

h3k4me3_no_pls = h3k4me3[no_pls,]
h3k27ac_no_pls = h3k27ac[no_pls,]
ctcf_1.64_no_pls = ctcf_1.64[no_pls,]
ctcf_em_no_pls = ctcf_em[no_pls,]

## get overlapped
c1 = c()
c2 = c()
for(i in 1:ncol(ctcf_1.64)){
  cellline = unlist(strsplit(colnames(ctcf_1.64)[i],"_"))[1]
  if(! cellline %in% colnames(h3k4me3)){
    c1 = c(c1,i)
  }
  if(! cellline %in% colnames(h3k27ac)){
    c2 = c(c2,i)
  }
}
ctcf_1.64_h3k4me3 <- ctcf_1.64[,-c1]
ctcf_1.64_h3k27ac <- ctcf_1.64[,-c2]

c3 = c()
c4 = c()
for(i in 1:ncol(ctcf_em)){
  cellline = unlist(strsplit(colnames(ctcf_em)[i],"_"))[1]
  if(! cellline %in% colnames(h3k4me3)){
    c3 = c(c3,i)
  }
  if(! cellline %in% colnames(h3k27ac)){
    c4 = c(c4,i)
  }
}
ctcf_em_h3k4me3 <- ctcf_em[,-c3]
ctcf_em_h3k27ac <- ctcf_em[,-c4]

ctcf_1.64_h3k4me3_pls <- ctcf_1.64_h3k4me3[pls,]
ctcf_1.64_h3k27ac_pls <- ctcf_1.64_h3k27ac[pls,]
ctcf_em_h3k4me3_pls <- ctcf_em_h3k4me3[pls,]
ctcf_em_h3k27ac_pls <- ctcf_em_h3k27ac[pls,]

ctcf_1.64_h3k4me3_no_pls <- ctcf_1.64_h3k4me3[no_pls,]
ctcf_1.64_h3k27ac_no_pls <- ctcf_1.64_h3k27ac[no_pls,]
ctcf_em_h3k4me3_no_pls <- ctcf_em_h3k4me3[no_pls,]
ctcf_em_h3k27ac_no_pls <- ctcf_em_h3k27ac[no_pls,]

###################
# H3K4me3
###################

h3k4me3_j = calculate_jaccard(h3k4me3)
colnames(h3k4me3_j) = rownames(h3k4me3_j) = colnames(h3k4me3)

pdf("ubi_rDHS_H3K4me3_Jaccard.pdf")
pheatmap(h3k4me3_j, display_numbers = F, main = "high H3K4me3 ubi-rDHS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - H3K4me3
###################

h3k4me3_pls_j = calculate_jaccard(h3k4me3_pls)
colnames(h3k4me3_pls_j) = rownames(h3k4me3_pls_j) = colnames(h3k4me3_pls)

pdf("ubi_rDHS_PLS_H3K4me3_Jaccard.pdf")
pheatmap(h3k4me3_pls_j, display_numbers = F, main = "high H3K4me3 ubi-rDHS-PLS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - H3K27ac
###################

h3k27ac_pls_j = calculate_jaccard(h3k27ac_pls)
colnames(h3k27ac_pls_j) = rownames(h3k27ac_pls_j) = colnames(h3k27ac_pls)

pdf("ubi_rDHS_PLS_H3K27ac_Jaccard.pdf")
pheatmap(h3k27ac_pls_j, display_numbers = F, main = "high H3K27ac ubi-rDHS-PLS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# H3K27ac
###################

h3k27ac_j = calculate_jaccard(h3k27ac)
colnames(h3k27ac_j) = rownames(h3k27ac_j) = colnames(h3k27ac)

pdf("ubi_rDHS_H3K27ac_Jaccard.pdf")
pheatmap(h3k27ac_j, display_numbers = F, main = "high H3K27ac ubi-rDHS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - H3K27ac
###################

h3k27ac_non_pls_j = calculate_jaccard(h3k27ac_no_pls)
colnames(h3k27ac_non_pls_j) = rownames(h3k27ac_non_pls_j) = colnames(h3k27ac_no_pls)

pdf("ubi_rDHS_nonPLS_H3K27ac_Jaccard.pdf")
pheatmap(h3k27ac_non_pls_j, display_numbers = F, main = "high H3K27ac ubi-rDHS-non-PLS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - H3K27ac
###################

h3k4me3_non_pls_j = calculate_jaccard(h3k4me3_no_pls)
colnames(h3k4me3_non_pls_j) = rownames(h3k4me3_non_pls_j) = colnames(h3k4me3_no_pls)

pdf("ubi_rDHS_nonPLS_H3K4me3_Jaccard.pdf")
pheatmap(h3k4me3_non_pls_j, display_numbers = F, main = "high H3K4me3 ubi-rDHS-non-PLS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# CTCF 1.64
###################

ctcf_1.64_j = calculate_jaccard(ctcf_1.64)
colnames(ctcf_1.64_j) = rownames(ctcf_1.64_j) = colnames(ctcf_1.64)

pdf("ubi_rDHS_CTCF_zscore_1.64_Jaccard.pdf", width = 12, height = 12)
pheatmap(ctcf_1.64_j, display_numbers = F, main = "high CTCF ubi-rDHS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# CTCF em
###################

ctcf_em_j = calculate_jaccard(ctcf_em)
colnames(ctcf_em_j) = rownames(ctcf_em_j) = colnames(ctcf_em)

pdf("ubi_rDHS_CTCF_zscore_em_Jaccard.pdf", width = 12, height = 12)
pheatmap(ctcf_em_j, display_numbers = F, main = "high CTCF ubi-rDHS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - CTCF 1.64
###################

ctcf_1.64_pls_j = calculate_jaccard(ctcf_1.64_pls)
colnames(ctcf_1.64_pls_j) = rownames(ctcf_1.64_pls_j) = colnames(ctcf_1.64_pls)

pdf("ubi-rDHS_PLS_CTCF_zscore_1.64_Jaccard.pdf", width=12, height = 12)
pheatmap(ctcf_1.64_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS-PLS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - CTCF em
###################

ctcf_em_pls_j = calculate_jaccard(ctcf_em_pls)
colnames(ctcf_em_pls_j) = rownames(ctcf_em_pls_j) = colnames(ctcf_em_pls)

pdf("ubi-rDHS_PLS_CTCF_em_Jaccard.pdf", height = 12, width = 12)
pheatmap(ctcf_em_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS-PLS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - CTCF 1.64
###################

ctcf_1.64_no_pls_j = calculate_jaccard(ctcf_1.64_no_pls)
colnames(ctcf_1.64_no_pls_j) = rownames(ctcf_1.64_no_pls_j) = colnames(ctcf_1.64_no_pls)

pdf("ubi-rDHS_nonPLS_CTCF_zscore_1.64_Jaccard.pdf", height = 12, width = 12)
pheatmap(ctcf_1.64_no_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS nonPLS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - CTCF em
###################

ctcf_em_no_pls_j = calculate_jaccard(ctcf_em_no_pls)
colnames(ctcf_em_no_pls_j) = rownames(ctcf_em_no_pls_j) = colnames(ctcf_em_no_pls)

pdf("ubi-rDHS_nonPLS_CTCF_em_Jaccard.pdf", height = 12, width = 12)
pheatmap(ctcf_em_no_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS nonPLS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# H3K27ac & CTCF 1.64
###################

ctcf_1.64_h3k27ac_j = calculate_jaccard_two(ctcf_1.64_h3k27ac, h3k27ac)
colnames(ctcf_1.64_h3k27ac_j) = rownames(ctcf_1.64_h3k27ac_j) = colnames(ctcf_1.64_h3k27ac)

pdf("ubi_rDHS_H3K27ac_CTCF_zscore_1.64_Jaccard.pdf", width = 12, height = 12)
pheatmap(ctcf_1.64_h3k27ac_j, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# H3K27ac & CTCF em
###################

ctcf_em_h3k27ac_j = calculate_jaccard_two(ctcf_em_h3k27ac, h3k27ac)
colnames(ctcf_em_h3k27ac_j) = rownames(ctcf_em_h3k27ac_j) = colnames(ctcf_em_h3k27ac)

pdf("ubi_rDHS_H3K27ac_CTCF_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_em_h3k27ac_j, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - H3K4me3 & CTCF 1.64
###################

ctcf_em_h3k4me3_pls_j = calculate_jaccard_two(ctcf_em_h3k4me3_pls, h3k4me3_pls)
colnames(ctcf_em_h3k4me3_pls_j) = rownames(ctcf_em_h3k4me3_pls_j) = colnames(ctcf_em_h3k4me3_pls)

pdf("ubi_rDHS_PLS_H3K4me3_CTCF_zscore_1.64_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_em_h3k4me3_pls_j, display_numbers = F, main = "high H3K4me3&CTCF ubi-rDHS PLS(zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - H3K4me3 & CTCF em
###################

ctcf_em_h3k4me3_pls_j = calculate_jaccard_two(ctcf_em_h3k4me3_pls, h3k4me3_pls)
colnames(ctcf_em_h3k4me3_pls_j) = rownames(ctcf_em_h3k4me3_pls_j) = colnames(ctcf_em_h3k4me3_pls)

pdf("ubi_rDHS_PLS_H3K4me3_CTCF_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_em_h3k4me3_pls_j, display_numbers = F, main = "high H3K4me3&CTCF ubi-rDHS PLS(EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - H3K27ac & CTCF 1.64
###################

ctcf_1.64_h3k27ac_no_pls_j = calculate_jaccard_two(ctcf_1.64_h3k27ac_no_pls, h3k27ac_no_pls)
colnames(ctcf_1.64_h3k27ac_no_pls_j) = rownames(ctcf_1.64_h3k27ac_no_pls_j) = colnames(ctcf_1.64_h3k27ac_no_pls)

pdf("ubi_rDHS_nonPLS_H3K27ac_CTCF_zscore_1.64_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_1.64_h3k27ac_no_pls_j, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS nonPLS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - H3K27ac & CTCF em
###################

ctcf_em_h3k27ac_no_pls_j = calculate_jaccard_two(ctcf_em_h3k27ac_no_pls, h3k27ac_no_pls)
colnames(ctcf_em_h3k27ac_no_pls_j) = rownames(ctcf_em_h3k27ac_no_pls_j) = colnames(ctcf_em_h3k27ac_no_pls)

pdf("ubi_rDHS_nonPLS_H3K27ac_CTCF_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_em_h3k27ac_no_pls_j, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS nonPLS(EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################