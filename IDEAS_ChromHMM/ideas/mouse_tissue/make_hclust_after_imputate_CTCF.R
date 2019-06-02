
# -- Kaili
# This script is for making state hlcust.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/new_states_after_CTCF/")

library(pheatmap)
library(RColorBrewer)

# read .para file, get only mark aveSignal
read_para_file <- function(file, filename){
  x=read.table(file, comment="!", header=T);
  k=dim(x)[2];
  l=dim(x)[1];
  p=(sqrt(9+8*(k-1))-3)/2;
  m=as.matrix(x[,1+1:p]/x[,1]);
  colnames(m) = colnames(x)[1+1:p];
  marks=colnames(m);
  rownames(m)=paste(filename, 1:l-1," (",round(x[,1]/sum(x[,1])*10000)/100,"%)",sep="");
  m_sort = m[,order(colnames(m))]
  return(m_sort)
}
################

data1 = read_para_file("all_data_44states.para0","noCTCF_")
data0 = read_para_file("all_data_impuatation.para0","withCTCF_")
data2 = data0[,c(1,3:11)]

matrix = rbind(data1, data2)
col_bar = data.frame(matrix("noCTCF",nrow=nrow(matrix)))
col_bar = transform(col_bar, CTCF_signal = 0)
colnames(col_bar) = c("has_CTCF", "CTCF_signal")
rownames(col_bar) = rownames(matrix)
col_bar$has_CTCF = as.vector(col_bar$has_CTCF)
col_bar[rownames(data2),]$has_CTCF = "withCTCF"
col_bar[rownames(data1),]$CTCF_signal = -10
col_bar[rownames(data2),]$CTCF_signal = data0[rownames(data2),2]
ann_colors=list(has_CTCF=c(noCTCF="#4d9221",withCTCF="#c51b7d"), 
                CTCF_signal = c(brewer.pal(5,"Greys")[5:1],brewer.pal(9,"Purples")))


pdf("state_after_CTCF_imputation_hlcust.pdf", height = 7, width = 15)
pheatmap(t(matrix), cluster_cols=T, cluster_rows = F, annotation_col=col_bar, 
         annotation_colors=ann_colors, main = "state_comparison_after_adding_CTCF", cellheight = 12, 
         cellwidth=10, clustering_distance_cols = "euclidean")
dev.off()






