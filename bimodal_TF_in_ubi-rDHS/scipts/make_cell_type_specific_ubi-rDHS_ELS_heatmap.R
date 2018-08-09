
# -- Kaili
# This script is for making cell-type specific ubi-rDHS ELS heatmap
# INPUT
# OUTPUT


setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/")

library(pheatmap)
library(RColorBrewer)
library(dplyr)

############
# function
############

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
ctcf_1.64 <- read.table("hg19_ubi-rDHS_CTCF_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
ctcf_em <- read.table("hg19_ubi-rDHS_CTCF_zscore_classification_ccreid.txt", header = TRUE, row.names = 1)

els <- read.table("hg19_cell_type_specific_ubi-rDHS_ELS_matrix_ccreid.txt", header = TRUE, row.names = 1)

##
c1 = c()
for(i in 1:ncol(ctcf_1.64)){
  cellline = unlist(strsplit(colnames(ctcf_1.64)[i],"_"))[1]
  if(! cellline %in% colnames(els)){
    c1 = c(c1,i)
  }
}
ctcf_1.64_els <- ctcf_1.64[,-c1]

c2 = c()
for(i in 1:ncol(ctcf_em)){
  cellline = unlist(strsplit(colnames(ctcf_em)[i],"_"))[1]
  if(! cellline %in% colnames(els)){
    c2 = c(c2,i)
  }
}
ctcf_em_els <- ctcf_em[,-c2]

###################
# main
###################

ctcf_1.64_els_j = calculate_jaccard_two(ctcf_1.64_els, els)
colnames(ctcf_1.64_els_j) = rownames(ctcf_1.64_els_j) = colnames(ctcf_1.64_els)
#
pdf("ubi_rDHS_cell_type_ELS_CTCF_zscore_1.64_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_1.64_els_j, display_numbers = F, main = "cell-type specific ubi-rDHS ELS\n CTCF (zscore>1.64)",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()


ctcf_em_els_j = calculate_jaccard_two(ctcf_em_els, els)
colnames(ctcf_em_els_j) = rownames(ctcf_em_els_j) = colnames(ctcf_em_els)
#
pdf("ubi_rDHS_cell_type_ELS_CTCF_zscore_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(ctcf_em_els_j, display_numbers = F, main = "cell-type specific ubi-rDHS ELS\n CTCF (EM)",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()


###################
