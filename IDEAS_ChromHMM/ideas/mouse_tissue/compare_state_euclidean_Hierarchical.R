
# -- Kaili
# This script is for compareing state files by doing Euclidean distance/hierarchical clustering.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/")
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



args<-commandArgs(TRUE)
file1 = args[1]
file2 = args[2]
title = args[3]
output = args[4]
type = args[5]

# get matrix
if(is.na(type)){
  matrix1 = read_para_file(file1)
  matrix2 = read_para_file(file2)
}else{
  matrix1 = data.frame(read_para_file(file1,"10marks "))
  matrix0 = read_para_file(file2, "11marks ")
  matrix2 = data.frame(matrix0[, colnames(matrix0)!=type])
}

matrix = rbind(matrix1, matrix2)
col_bar = data.frame(matrix("without_CTCF",nrow=nrow(matrix)))
col_bar = transform(col_bar, CTCF_signal = 0)
colnames(col_bar) = c("has_CTCF", "CTCF_signal")
rownames(col_bar) = rownames(matrix)
col_bar$has_CTCF = as.vector(col_bar$has_CTCF)
col_bar[rownames(matrix2),]$has_CTCF = "with_CTCF"
col_bar[rownames(matrix1),]$CTCF_signal = -10
col_bar[rownames(matrix2),]$CTCF_signal = matrix0[rownames(matrix2),2]
ann_colors=list(has_CTCF=c(without_CTCF="#4d9221",with_CTCF="#c51b7d"), 
                CTCF_signal = c(brewer.pal(5,"Greys")[5:1],brewer.pal(9,"Purples")))


pdf(output, height = 7, width = 15)
pheatmap(t(matrix), cluster_cols=T, cluster_rows = F, annotation_col=col_bar, annotation_colors=ann_colors, 
         main = title, cellheight = 12, cellwidth=10, clustering_distance_cols = "euclidean")
dev.off()


###########
# states color
plot.new()
legend("topright",legend=c("promoter", "enhancer","CTCF", "bivalent TSS","poised Enhancer", "Quies", "Enhancer in gene/Transcription", "H3K36me3+H3K9me3"), 
       col=c("red","yellow", "purple","blue","black","gray92", "green4", "bisque"), pch=15, cex=1, bty="n")

