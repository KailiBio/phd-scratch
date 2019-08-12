
# -- Kaili
# This script is for ranking states by average CTCF signal.

args = commandArgs(trailingOnly=TRUE)
workDir = args[1]
filename = args[2]
outfile = args[3]

setwd(workDir)
# setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/")

order_state_CTCF <- function(filename){
  x=read.table(filename, comment="!", header=T);
  k=dim(x)[2];
  l=dim(x)[1];
  p=(sqrt(9+8*(k-1))-3)/2;
  m=data.frame(as.matrix(x[,1+1:p]/x[,1]))
  colnames(m) = colnames(x)[1+1:p];
  marks=colnames(m);
  rownames(m)=1:l-1;
  m_sort = m[order(m$CTCF,decreasing = FALSE),]
  return(rownames(m_sort))
}

####
# filename = "ctcf_9sample_impute_11sample.para0"
state_order = data.frame(order_state_CTCF(filename))
write.table(state_order, outfile, quote = FALSE, sep="\t", col.names = FALSE)




