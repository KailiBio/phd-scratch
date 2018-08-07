
# -- Kaili
# This script is for classify bimodal distribution (GMM) data by EM algorithm
# INPUT: bimodal values (signal/zscore matrix here)
# OURPUT: threshold

# EXP: Rscript classify_bimodal_EM.R /data/zusers/fankaili/ccre/tf/matrix/hg19_ubi-rDHS_CTCF_signal_matrix.txt /data/zusers/fankaili/ccre/tf/matrix/ hg19_ubi-rDHS_CTCF_signal_log10_classification.txt
# /data/zusers/fankaili/ccre/tf/figs/hg19_ubi-rDHS_CTCF_signal_log10_classification.pdf log10

args<-commandArgs(T)
inFile = args[1]
outDir = args[2]
outMatrix = args[3]
outFigure = args[4]
type = args[5] ## log2, log10, zscore

library("mixtools")
setwd(outDir)
data <- read.table(inFile, header=TRUE, row.names = 1)

m = matrix(0, nrow = nrow(data), ncol = ncol(data))
rownames(m) = rownames(data)
colnames(m) = colnames(data)
pdf(outFigure)
for(i in 1:ncol(data)){
  biosample = colnames(data)[i]

  if(type=="log2"){
    dat = log2(data[,i]+0.01)
  } else if (type=="log10"){
    dat = log10(data[,i]+0.01)
  } else if(type=="zscore"){
    dat = data[,i]
  }
  names(dat) <- rownames(data)

  myEM <- normalmixEM(dat)
  plot(myEM, whichplots=2, sub = biosample)
  p <- myEM$posterior

  m[names(dat[p[,1]>p[,2]]),i]=1
}
dev.off()

write.table(m, file = outMatrix, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
