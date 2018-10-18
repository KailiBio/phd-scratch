
# -- Kaili
# This script is for making correlation between mean_siganl states between two IDEAS runs.
# EXP: Rscript make_state_correlation.R "e14.5p0_8hm_ATAC_DNAme.para0" "e14.5p0_8hm_ATAC_DNAme_CTCF_1.para0"
#      "21_biosamples_without_CTCF" "21_biosamples_with_CTCF" "21 biosamples in e14.5&p0, 8HM+ATAC+DNAme"
#      "e14.5p0_correlation.pdf" "CTCF"

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/")

#####################
# function
#####################

# read .para file, get only mark aveSignal
read_para_file <- function(file){
  x=read.table(file, comment="!", header=T);
  k=dim(x)[2];
  l=dim(x)[1];
  p=(sqrt(9+8*(k-1))-3)/2;
  m=as.matrix(x[,1+1:p]/x[,1]);
  colnames(m) = colnames(x)[1+1:p];
  marks=colnames(m);
  rownames(m)=paste(1:l-1," (",round(x[,1]/sum(x[,1])*10000)/100,"%)",sep="");
  m_sort = m[,order(colnames(m))]
  return(m_sort)
}

# calculate state correlation between two IDEAS runs
calculate_correlation <- function(matrix1, matrix2){
  cor_matrix = matrix(0, nrow = nrow(matrix1), ncol = nrow(matrix2))
  #
  for(i in 1:nrow(matrix1)){
    for(j in 1:nrow(matrix2)){
      cor_matrix[i,j] = cor(as.numeric(matrix1[i,]), as.numeric(matrix2[j,]), method = "pearson")
    }
  }
  rownames(cor_matrix) = rownames(matrix1)
  colnames(cor_matrix) = rownames(matrix2)
  return(cor_matrix)
}

#####################
# file1 = "e14.5p0_8hm_ATAC_DNAme.para0"
# file2 = "e14.5p0_8hm_ATAC_DNAme_CTCF_1.para0"
# name1 = "21_biosamples_without_CTCF"
# name2 = "21_biosamples_with_CTCF"
# title = "21 biosamples in e14.5&p0, 8HM+ATAC+DNAme\n(39 states vs. 37 states)"
# output = "e14.5p0_correlation.pdf"
# type = "CTCF"
# 
# file1 = "run_IDEAS_8hm_atac_dname_pvalue.para0" 
# file2 = "e14.5p0_8hm_ATAC_DNAme_CTCF_1.para0" 
# name1 = "66_biosamples_without_CTCF" 
# name2 = "21_biosamples_with_CTCF" 
# title = "states comparison(38 states vs. 37 states)" 
# output = "66samples_21samplesCTCF_correlation.pdf" 
# type = "CTCF"

args<-commandArgs(TRUE)
file1 = args[1]
file2 = args[2]
name1 = args[3]
name2 = args[4]
title = args[6]
output = args[7]
type = args[8]

# get matrix
if(is.na(type)){
  matrix1 = read_para_file(file1)
  matrix2 = read_para_file(file2)
}else{
  matrix1 = data.frame(read_para_file(file1))
  matrix0 = read_para_file(file2)
  matrix2 = data.frame(matrix0[, colnames(matrix0)!=type])
}

# calculate correlation
cor_matrix = calculate_correlation(matrix1, matrix2)

# make heatmap
pdf(output)
library(pheatmap)
library(RColorBrewer)
library(grid)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.92, height=0.92, name="vp", just=c("right","top"))), action="prepend")
pheatmap(cor_matrix, cluster_rows=T, cluster_cols=F, breaks = c(-1,seq(0,0.8,0.1),0.9,1), 
         col = c(colorRampPalette(brewer.pal(7, "Greys")[7:2])(10), "#FFFF00FF"), main = title)
setHook("grid.newpage", NULL, "replace")
grid.text(name1, y=-0.04, gp=gpar(fontsize=12))
grid.text(name2, x=-0.04, rot=90, gp=gpar(fontsize=12))
dev.off()

#####################