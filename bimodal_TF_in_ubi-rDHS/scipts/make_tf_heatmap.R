
# -- Kaili
# This script is for making heatmap for TF data (CTCF, SMC3, RAD21)
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

calculate_jaccard_two2 <- function(data1, data2){
  n = ncol(data1)
  out = matrix(0,n,n)
  for(i in 1:n){
    p = unlist(strsplit(colnames(data1)[i],"_"))[2]
    a = count(data1[,i]+data2[,p])
    a2 = a[a$x==2,]$freq
    a_name = rownames(data1[(data1[,i]+data2[,p])==2,])
    for(j in 1:n){
      q = unlist(strsplit(colnames(data1)[j],"_"))[2]
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

matrix = read.table("hg19_ubi-rDHS_TF_zscore_em_classification_ccreid.txt", header = TRUE, row.names = 1)
els <- read.table("hg19_cell_type_specific_ubi-rDHS_ELS_matrix_ccreid.txt", header = TRUE, row.names = 1)

h3k27ac <- read.table("hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
h3k4me3 <- read.table("hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)

pls <- as.vector(read.table("hg19_ubi-rDHS_PLS_list.txt")[,1])
no_pls <- as.vector(read.table("hg19_ubi-rDHS_non_PLS_list.txt")[,1])


###################
# ubi-rDHS
###################
ubi_j = calculate_jaccard(matrix)
colnames(ubi_j) = rownames(ubi_j) = colnames(matrix)

pdf("ubi_rDHS_TF_zscore_em_Jaccard.pdf", width = 12, height = 12)
pheatmap(ubi_j, display_numbers = F, main = "high TF ubi-rDHS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# ubi-rDHS & H3K27ac
###################
cc = c()
for(i in 1:ncol(matrix)){
  cellline = unlist(strsplit(colnames(matrix)[i],"_"))[2]
  if(! cellline %in% colnames(h3k27ac)){
    cc = c(cc,i)
  }
}
matrix_h3k27ac <- matrix[,-cc]

ubi_ac_j = calculate_jaccard_two2(matrix_h3k27ac, h3k27ac)
colnames(ubi_ac_j) = rownames(ubi_ac_j) = colnames(matrix_h3k27ac)

pdf("ubi_rDHS_H3K27ac_TF_zscore_em_Jaccard.pdf", width = 12, height = 12)
pheatmap(ubi_ac_j, display_numbers = F, main = "high H3K27ac&TF ubi-rDHS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()



###################
# ubi-rDHS PLS
###################
matrix_pls <- matrix[pls,]

ubi_pls_j = calculate_jaccard(matrix_pls)
colnames(ubi_pls_j) = rownames(ubi_pls_j) = colnames(matrix_pls)

pdf("ubi_rDHS_PLS_TF_zscore_em_Jaccard.pdf", width = 12, height = 12)
pheatmap(ubi_pls_j, display_numbers = F, main = "high TF ubi-rDHS PLS(EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# ubi-rDHS nonPLS
###################
matrix_no_pls <- matrix[no_pls,]

ubi_no_pls_j = calculate_jaccard(matrix_no_pls)
colnames(ubi_no_pls_j) = rownames(ubi_no_pls_j) = colnames(matrix_no_pls)

pdf("ubi_rDHS_non_PLS_TF_zscore_em_Jaccard.pdf", width = 12, height = 12)
pheatmap(ubi_no_pls_j, display_numbers = F, main = "high TF ubi-rDHS non_PLS(EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# cell-tytpe specific ubi-rDHS ELS
###################
cc = c()
for(i in 1:ncol(matrix)){
  cellline = unlist(strsplit(colnames(matrix)[i],"_"))[2]
  if(! cellline %in% colnames(els)){
    cc = c(cc,i)
  }
}
matrix_els <- matrix[,-cc]

tf_els_j = calculate_jaccard_two2(matrix_els, els)
colnames(tf_els_j) = rownames(tf_els_j) = colnames(matrix_els)
#
pdf("ubi_rDHS_cell_type_ELS_TF_zscore_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(tf_els_j, display_numbers = F, main = "cell-type specific ubi-rDHS ELS\n TF (EM)",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################