
# -- Kaili
# This script is for compareing state files by doing cos correlation.

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


#####################
# file1 = "e14.5p0_8hm_ATAC_DNAme.para0"
# file2 = "e14.5p0_8hm_ATAC_DNAme_CTCF_1.para0"
# name1 = "21_biosamples"
# name2 = "21_biosamples + CTCF"
# title = "21 biosamples in e14.5&p0, 8HM+ATAC+DNAme\n(39 states vs. 37 states)"
# output = "e14.5p0_correlation_spearman.pdf"
# type = "CTCF"
# #
# file1 = "run_IDEAS_8hm_atac_dname_pvalue.para0"
# file2 = "e14.5p0_8hm_ATAC_DNAme_CTCF_1.para0"
# name1 = "66_biosamples"
# name2 = "21_biosamples + CTCF"
# title = "states comparison(38 states vs. 37 states)"
# output = "66samples_21samplesCTCF_correlation_spearman.pdf"
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
