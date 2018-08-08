
# -- Kaili
# This script is for making heatmap using given data (jaccard).
# locally
# INPUT
# OUTPUT

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/")

library(pheatmap)
library(RColorBrewer)

############
# data
############
h3k27ac <- read.table("hg19_ubi-rDHS_H3K27ac_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
h3k4me3 <- read.table("hg19_ubi-rDHS_H3K4me3_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
ctcf_1.64 <- read.table("hg19_ubi-rDHS_CTCF_zscore_1.64_classification_ccreid.txt", header = TRUE, row.names = 1)
ctcf_em <- read.table("hg19_ubi-rDHS_CTCF_zscore_classification_ccreid.txt", header = TRUE, row.names = 1)

pls <- as.vector(read.table("hg19_PLS_list.txt")[,1])
no_pls <- as.vector(read.table("hg19_non_PLS_list.txt")[,1])

h3k4me3_pls = h3k4me3[pls,]
h3k27ac_pls = h3k27ac[pls,]
ctcf_1.64_pls = ctcf_1.64[pls,]
ctcf_em_pls = ctcf_em[pls,]

h3k4me3_no_pls = h3k4me3[no_pls,]
h3k27ac_no_pls = h3k27ac[no_pls,]
ctcf_1.64_no_pls = ctcf_1.64[no_pls,]
ctcf_em_no_pls = ctcf_em[no_pls,]

###################
# H3K4me3
###################
n = ncol(h3k4me3)
h3k4me3_j = matrix(0,n,n)
colnames(h3k4me3_j) = rownames(h3k4me3_j) = colnames(h3k4me3)

for (i in 1:n){
  for (j in 1:n){
    a = nrow(h3k4me3[h3k4me3[,i]==1 & h3k4me3[,j]==1,])
    b = nrow(h3k4me3[h3k4me3[,i]==1 | h3k4me3[,j]==1,])
    jaccard = round(a/b,2)
    h3k4me3_j[i,j] = jaccard
  }
}

pdf("ubi_rDHS_H3K4me3_Jaccard.pdf")
pheatmap(h3k4me3_j, display_numbers = F, main = "high H3K4me3 ubi-rDHS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# pls - H3K4me3
###################
n = ncol(h3k4me3_pls)
h3k4me3_pls_j = matrix(0,n,n)
colnames(h3k4me3_pls_j) = rownames(h3k4me3_pls_j) = colnames(h3k4me3_pls)

for (i in 1:n){
  for (j in 1:n){
    a = nrow(h3k4me3_pls[h3k4me3_pls[,i]==1 & h3k4me3_pls[,j]==1,])
    b = nrow(h3k4me3_pls[h3k4me3_pls[,i]==1 | h3k4me3_pls[,j]==1,])
    jaccard = round(a/b,2)
    h3k4me3_pls_j[i,j] = jaccard
  }
}

pdf("ubi_rDHS_PLS_H3K4me3_Jaccard.pdf")
pheatmap(h3k4me3_pls_j, display_numbers = F, main = "high H3K4me3 ubi-rDHS-PLS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()


###################
# H3K27ac
###################
n = ncol(h3k27ac)
h3k27ac_j = matrix(0,n,n)
colnames(h3k27ac_j) = rownames(h3k27ac_j) = colnames(h3k27ac)

for (i in 1:n){
  for (j in 1:n){
    a = nrow(h3k27ac[h3k27ac[,i]==1 & h3k27ac[,j]==1,])
    b = nrow(h3k27ac[h3k27ac[,i]==1 | h3k27ac[,j]==1,])
    jaccard = round(a/b,2)
    h3k27ac_j[i,j] = jaccard
  }
}

pdf("ubi_rDHS_H3K27ac_Jaccard.pdf")
pheatmap(h3k27ac_j, display_numbers = F, main = "high H3K27ac ubi-rDHS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# no pls - H3K27ac
###################
n = ncol(h3k27ac_no_pls)
h3k27ac_non_pls_j = matrix(0,n,n)
colnames(h3k27ac_non_pls_j) = rownames(h3k27ac_non_pls_j) = colnames(h3k27ac_no_pls)

for (i in 1:n){
  for (j in 1:n){
    a = nrow(h3k27ac_no_pls[h3k27ac_no_pls[,i]==1 & h3k27ac_no_pls[,j]==1,])
    b = nrow(h3k27ac_no_pls[h3k27ac_no_pls[,i]==1 | h3k27ac_no_pls[,j]==1,])
    jaccard = round(a/b,2)
    h3k27ac_non_pls_j[i,j] = jaccard
  }
}

pdf("ubi_rDHS_nonPLS_H3K27ac_Jaccard.pdf")
pheatmap(h3k27ac_non_pls_j, display_numbers = F, main = "high H3K27ac ubi-rDHS-non-PLS\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# CTCF 1.64
###################
n = ncol(ctcf_1.64)
ctcf_1.64_j = matrix(0,n,n)
colnames(ctcf_1.64_j) = rownames(ctcf_1.64_j) = colnames(ctcf_1.64)

for (i in 1:n){
  for (j in 1:n){
    a = nrow(ctcf_1.64[ctcf_1.64[,i]==1 & ctcf_1.64[,j]==1,])
    b = nrow(ctcf_1.64[ctcf_1.64[,i]==1 | ctcf_1.64[,j]==1,])
    jaccard = round(a/b,2)
    ctcf_1.64_j[i,j] = jaccard
  }
}

pdf("ubi_rDHS_CTCF_zscore_1.64_Jaccard.pdf", width = 12, height = 12)
pheatmap(ctcf_1.64_j, display_numbers = F, main = "high CTCF ubi-rDHS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()
###################
# CTCF em
###################
n = ncol(ctcf_em)
ctcf_em_j = matrix(0,n,n)
colnames(ctcf_em_j) = rownames(ctcf_em_j) = colnames(ctcf_em)

for (i in 1:n){
  for (j in 1:n){
    a = nrow(ctcf_em[ctcf_em[,i]==1 & ctcf_em[,j]==1,])
    b = nrow(ctcf_em[ctcf_em[,i]==1 | ctcf_em[,j]==1,])
    jaccard = round(a/b,2)
    ctcf_em_j[i,j] = jaccard
  }
}

pdf("ubi_rDHS_CTCF_zscore_em_Jaccard.pdf", width = 12, height = 12)
pheatmap(ctcf_em_j, display_numbers = F, main = "high CTCF ubi-rDHS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()

###################
# H3K27ac & CTCF 1.64
###################
n = ncol(ctcf_1.64)
h3k27ac_ctcf_jaccard = matrix(0, nrow=n, ncol = n)
colnames(h3k27ac_ctcf_jaccard) = rownames(h3k27ac_ctcf_jaccard) = colnames(ctcf_1.64)

for (i in 1:n){
  for (j in 1:n){
    cell1 = unlist(strsplit(colnames(ctcf_1.64)[i],split="_"))[1]
    cell2 = unlist(strsplit(colnames(ctcf_1.64)[j],split="_"))[1]
    a=intersect(intersect(rownames(h3k27ac[,cell1]==1), rownames(ctcf_1.64[,i]==2)))
    a = nrow(ctcf_1.64[(h3k27ac[,cell1]==1 & ctcf_1.64[,i]==2) & (h3k27ac[,cell2]==1 & ctcf_1.64[,j]==2),])
    b = nrow(ctcf_1.64[(h3k27ac[,cell1]==1 & ctcf_1.64[,i]==2) | (h3k27ac[,cell2]==1 & ctcf_1.64[,j]==2),])
    jaccard = round(a/b,2)
    h3k27ac_ctcf_jaccard[i,j] = jaccard
  }
}

pdf("ubi_rDHS_H3K27ac_CTCF_zscore_1.64_Jaccard.pdf", , width = 12, height = 12)
pheatmap(h3k27ac_ctcf_jaccard, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()
###################
# H3K27ac & CTCF em
###################
n = ncol(ctcf_em)
h3k27ac_ctcf_jaccard = matrix(0, nrow=n, ncol = n)
colnames(h3k27ac_ctcf_jaccard) = rownames(h3k27ac_ctcf_jaccard) = colnames(ctcf_em)

for (i in 1:n){
  for (j in 1:n){
    cell1 = unlist(strsplit(colnames(ctcf_em)[i],split="_"))[1]
    cell2 = unlist(strsplit(colnames(ctcf_em)[j],split="_"))[1]
    a = nrow(ctcf_em[(h3k27ac[,cell1]==1 & ctcf_em[,i]==2) & (h3k27ac[,cell2]==1 & ctcf_em[,j]==2),])
    b = nrow(ctcf_em[(h3k27ac[,cell1]==1 & ctcf_em[,i]==2) | (h3k27ac[,cell2]==1 & ctcf_em[,j]==2),])
    jaccard = round(a/b,2)
    h3k27ac_ctcf_jaccard[i,j] = jaccard
  }
}

pdf("ubi_rDHS_H3K27ac_CTCF_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(h3k27ac_ctcf_jaccard, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()



###################
# pls - CTCF 1.64
###################
n = ncol(ctcf_1.64_pls)
ctcf_1.64_pls_j = matrix(0,n,n)
colnames(ctcf_1.64_pls_j) = rownames(ctcf_1.64_pls_j) = colnames(ctcf_1.64_pls)

for (i in 1:n){
  for (j in 1:n){
    a = nrow(ctcf_1.64_pls[ctcf_1.64_pls[,i]==1 & ctcf_1.64_pls[,j]==1,])
    b = nrow(ctcf_1.64_pls[ctcf_1.64_pls[,i]==1 | ctcf_1.64_pls[,j]==1,])
    jaccard = round(a/b,2)
    ctcf_1.64_pls_j[i,j] = jaccard
  }
}

pdf("ubi_ubi-rDHS_PLS_CTCF_zscore_1.64_Jaccard.pdf", width=12, height = 12)
pheatmap(ctcf_1.64_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS-PLS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()
###################
# pls - CTCF em
###################
n = ncol(ctcf_em_pls)
ctcf_em_pls_j = matrix(0,n,n)
colnames(ctcf_em_pls_j) = rownames(ctcf_em_pls_j) = colnames(ctcf_em_pls)

for (i in 1:n){
  for (j in 1:n){
    a = nrow(ctcf_em_pls[ctcf_em_pls[,i]==1 & ctcf_em_pls[,j]==1,])
    b = nrow(ctcf_em_pls[ctcf_em_pls[,i]==1 | ctcf_em_pls[,j]==1,])
    jaccard = round(a/b,2)
    ctcf_em_pls_j[i,j] = jaccard
  }
}

pdf("ubi_ubi-rDHS_PLS_CTCF_em_Jaccard.pdf")
pheatmap(ctcf_em_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS-PLS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()
###################
# pls - H3K4me3 & CTCF 1.64
###################
n = ncol(ctcf_1.64_pls)
h3k4me3_ctcf_jaccard = matrix(0, nrow=n, ncol = n)
colnames(h3k4me3_ctcf_jaccard) = rownames(h3k4me3_ctcf_jaccard) = colnames(ctcf_1.64_pls)

for (i in 1:n){
  for (j in 1:n){
    cell1 = unlist(strsplit(colnames(ctcf_1.64_pls)[i],split="_"))[1]
    cell2 = unlist(strsplit(colnames(ctcf_1.64_pls)[j],split="_"))[1]
    a = nrow(ctcf_1.64_pls[(h3k4me3[,cell1]==1 & ctcf_1.64_pls[,i]==2) & (h3k4me3[,cell2]==1 & ctcf_1.64_pls[,j]==2),])
    b = nrow(ctcf_1.64_pls[(h3k4me3[,cell1]==1 & ctcf_1.64_pls[,i]==2) | (h3k4me3[,cell2]==1 & ctcf_1.64_pls[,j]==2),])
    jaccard = round(a/b,2)
    h3k4me3_ctcf_jaccard[i,j] = jaccard
  }
}

pdf("ubi_rDHS_PLS_H3K4me3_CTCF_zscore_1.64_Jaccard.pdf", height = 10, width = 10)
pheatmap(h3k4me3_ctcf_jaccard, display_numbers = F, main = "high H3K4me3&CTCF ubi-rDHS PLS(zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()
###################
# pls - H3K4me3 & CTCF em
###################
n = ncol(ctcf_em_pls)
h3k4me3_ctcf_jaccard = matrix(0, nrow=n, ncol = n)
colnames(h3k4me3_ctcf_jaccard) = rownames(h3k4me3_ctcf_jaccard) = colnames(ctcf_em_pls)

for (i in 1:n){
  for (j in 1:n){
    cell1 = unlist(strsplit(colnames(ctcf_em_pls)[i],split="_"))[1]
    cell2 = unlist(strsplit(colnames(ctcf_em_pls)[j],split="_"))[1]
    a = nrow(ctcf_em_pls[(h3k4me3[,cell1]==1 & ctcf_em_pls[,i]==2) & (h3k4me3[,cell2]==1 & ctcf_em_pls[,j]==2),])
    b = nrow(ctcf_em_pls[(h3k4me3[,cell1]==1 & ctcf_em_pls[,i]==2) | (h3k4me3[,cell2]==1 & ctcf_em_pls[,j]==2),])
    jaccard = round(a/b,2)
    h3k4me3_ctcf_jaccard[i,j] = jaccard
  }
}

pdf("ubi_rDHS_PLS_H3K4me3_CTCF_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(h3k4me3_ctcf_jaccard, display_numbers = F, main = "high H3K4me3&CTCF ubi-rDHS PLS(EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()



###################
# no pls - CTCF 1.64
###################
n = ncol(ctcf_1.64_no_pls)
ctcf_1.64_no_pls_j = matrix(0,n,n)
colnames(ctcf_1.64_no_pls_j) = rownames(ctcf_1.64_no_pls_j) = colnames(ctcf_1.64_no_pls)

for (i in 1:n){
  for (j in 1:n){
    a = nrow(ctcf_1.64_no_pls[ctcf_1.64_no_pls[,i]==1 & ctcf_1.64_no_pls[,j]==1,])
    b = nrow(ctcf_1.64_no_pls[ctcf_1.64_no_pls[,i]==1 | ctcf_1.64_no_pls[,j]==1,])
    jaccard = round(a/b,2)
    ctcf_1.64_no_pls_j[i,j] = jaccard
  }
}

pdf("ubi_ubi-rDHS_nonPLS_CTCF_zscore_1.64_Jaccard.pdf")
pheatmap(ctcf_1.64_no_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS nonPLS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()
###################
# no pls - CTCF em
###################
n = ncol(ctcf_em_no_pls)
ctcf_em_no_pls_j = matrix(0,n,n)
colnames(ctcf_em_no_pls_j) = rownames(ctcf_em_no_pls_j) = colnames(ctcf_em_no_pls)

for (i in 1:n){
  for (j in 1:n){
    a = nrow(ctcf_em_no_pls[ctcf_em_no_pls[,i]==1 & ctcf_em_no_pls[,j]==1,])
    b = nrow(ctcf_em_no_pls[ctcf_em_no_pls[,i]==1 | ctcf_em_no_pls[,j]==1,])
    jaccard = round(a/b,2)
    ctcf_em_no_pls_j[i,j] = jaccard
  }
}

pdf("ubi_ubi-rDHS_nonPLS_CTCF_em_Jaccard.pdf")
pheatmap(ctcf_em_no_pls_j, display_numbers = F, main = "high CTCF ubi-rDHS nonPLS (EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()




###################
# no pls - H3K27ac & CTCF 1.64
###################
n = ncol(ctcf_1.64_no_pls)
h3k27ac_ctcf_jaccard = matrix(0, nrow=n, ncol = n)
colnames(h3k27ac_ctcf_jaccard) = rownames(h3k27ac_ctcf_jaccard) = colnames(ctcf_1.64_no_pls)

for (i in 1:n){
  for (j in 1:n){
    cell1 = unlist(strsplit(colnames(ctcf_1.64_no_pls)[i],split="_"))[1]
    cell2 = unlist(strsplit(colnames(ctcf_1.64_no_pls)[j],split="_"))[1]
    a = nrow(ctcf_1.64_no_pls[(h3k27ac[,cell1]==1 & ctcf_1.64_no_pls[,i]==2) & (h3k27ac[,cell2]==1 & ctcf_1.64_no_pls[,j]==2),])
    b = nrow(ctcf_1.64_no_pls[(h3k27ac[,cell1]==1 & ctcf_1.64_no_pls[,i]==2) | (h3k27ac[,cell2]==1 & ctcf_1.64_no_pls[,j]==2),])
    jaccard = round(a/b,2)
    h3k27ac_ctcf_jaccard[i,j] = jaccard
  }
}

pdf("ubi_rDHS_nonPLS_H3K27ac_CTCF_zscore_1.64_Jaccard.pdf", height = 10, width = 10)
pheatmap(h3k27ac_ctcf_jaccard, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS nonPLS (zscore>1.64)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()
###################
# no pls - H3K27ac & CTCF em
###################
n = ncol(ctcf_em_no_pls)
h3k27ac_ctcf_jaccard = matrix(0, nrow=n, ncol = n)
colnames(h3k27ac_ctcf_jaccard) = rownames(h3k27ac_ctcf_jaccard) = colnames(ctcf_em_no_pls)

for (i in 1:n){
  for (j in 1:n){
    cell1 = unlist(strsplit(colnames(ctcf_em_no_pls)[i],split="_"))[1]
    cell2 = unlist(strsplit(colnames(ctcf_em_no_pls)[j],split="_"))[1]
    a = nrow(ctcf_em_no_pls[(h3k27ac[,cell1]==1 & ctcf_em_no_pls[,i]==2) & (h3k27ac[,cell2]==1 & ctcf_em_no_pls[,j]==2),])
    b = nrow(ctcf_em_no_pls[(h3k27ac[,cell1]==1 & ctcf_em_no_pls[,i]==2) | (h3k27ac[,cell2]==1 & ctcf_em_no_pls[,j]==2),])
    jaccard = round(a/b,2)
    h3k27ac_ctcf_jaccard[i,j] = jaccard
  }
}

pdf("ubi_rDHS_nonPLS_H3K27ac_CTCF_em_Jaccard.pdf", height = 10, width = 10)
pheatmap(h3k27ac_ctcf_jaccard, display_numbers = F, main = "high H3K27ac&CTCF ubi-rDHS nonPLS(EM)\nJaccard Index",
         breaks = (0:10)/10, color = colorRampPalette(brewer.pal(9,"RdYlBu")[9:1])(10))
dev.off()