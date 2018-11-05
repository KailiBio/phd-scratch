
# -- Kaili
# This script is for compareing state files by doing Euclidean distance/hierarchical clustering.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/")
setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/")
library(pheatmap)
library(RColorBrewer)

#####################
# function
#####################

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


#####################
# file1 = "e14.5p0_8hm_ATAC_DNAme.para0"
# file2 = "e14.5p0_8hm_ATAC_DNAme_CTCF_1.para0"
# title = "21 biosamples in e14.5&p0, 8HM+ATAC+DNAme\n(39 states vs. 37 states)"
# output = "e14.5p0_hcluster.pdf"
# type = "CTCF"
# #
# file1 = "run_IDEAS_8hm_atac_dname_pvalue.para0"
# file2 = "e14.5p0_8hm_ATAC_DNAme_CTCF_1.para0"
# title = "states comparison(38 states vs. 37 states)"
# output = "66samples_21samplesCTCF_hcluster.pdf"
# type = "CTCF"
#
file1 = "run_IDEAS_8hm_atac_dname_pvalue.para0"
file2 = "DHS_v1_100-300bp.para0"
title = "states comparison: 100-300bp(38 states vs. 36 states)"
output = "dhs_bins_v1_comparison_hcluster.pdf"
type=""
#
file1 = "run_IDEAS_8hm_atac_dname_pvalue.para0"
file2 = "DHS_v2_1-300bp.para0"
title = "states comparison: 1-300bp(38 states vs. 36 states)"
output = "dhs_bins_v2_comparison_hcluster.pdf"
type=""


args<-commandArgs(TRUE)
file1 = args[1]
file2 = args[2]
title = args[3]
output = args[4]
type = args[5]

# get matrix
matrix1 = read_para_file(file1, "normal_")
matrix2 = read_para_file(file2, "dhs_")

matrix = rbind(matrix1, matrix2)
col_bar = data.frame(matrix("normal_bins",nrow=nrow(matrix)))
colnames(col_bar) = c("bins")
rownames(col_bar) = rownames(matrix)
col_bar$bins = as.vector(col_bar$bins)
col_bar[rownames(matrix2),1] = "dhs_bins"
ann_colors=list(bins=c(normal_bins="#4d9221",dhs_bins="#c51b7d"))


pdf(output, height = 7, width = 15)
pheatmap(t(matrix), cluster_cols=T, cluster_rows = F, annotation_col=col_bar, annotation_colors=ann_colors, 
         main = title, cellheight = 12, cellwidth=10, clustering_distance_cols = "euclidean")
dev.off()


