
# -- Kaili
# This script is for classify bimodal distribution (GMM) data by EM algorithm
# INPUT: bimodal values (signal/zscore matrix here)
# OURPUT: threshold

# EXP: Rscript classify_bimodal_EM.R /data/zusers/fankaili/ccre/tf/matrix/ CTCF zscore

args<-commandArgs(T)
outDir = args[1]
tf = args[2] ## CTCF, SMC3, RAD21 et, al.
type = args[5] ## log2, log10, zscore

inFile = paste(outDir,"hg19_ubi-rDHS_",tf,"_zscore_matrix.txt", sep="")
outMatrix = paste("hg19_ubi-rDHS_",tf,"_zscore_classification.txt", sep="")
outFigure = paste("/data/zusers/fankaili/ccre/tf/figs/","hg19_ubi-rDHS_",tf,"_zscore_classification.pdf", sep="")

library("mixtools")
setwd(outDir)
data <- read.table(inFile, header=TRUE, row.names = 1)

m = matrix(0, nrow = nrow(data), ncol = ncol(data))
rownames(m) = rownames(data)
colnames(m) = colnames(data)
pdf(outFigure)
for(i in 1:ncol(data)){
  biosample = colnames(data)[i]

  # get data
  if(type=="log2"){
    dat = log2(data[,i]+0.01)
  } else if (type=="log10"){
    dat = log10(data[,i]+0.01)
  } else if(type=="zscore"){
    dat = data[data[,i]!=(-10),i]
    # remove outliner (-10) here
  }
  names(dat) <- rownames(data[data[,i]!=(-10),])

  # call EM
  if(tf=="CTCF"){
    myEM <- normalmixEM(dat, mu = c(min(dat),max(dat)), sigma = c(2,1))
  }else{
    myEM <- normalmixEM(dat, mu = c(min(dat),max(dat)), sigma = c(4,1))
  }
  
  # density plot
  hist(dat, breaks=200, freq = FALSE, main=biosample, xlab="zscore")
  curve((myEM$lambda[1]*dnorm(x, myEM$mu[1], myEM$sigma[1])), col="green", lwd=3, add=TRUE)
  curve((myEM$lambda[2]*dnorm(x, myEM$mu[2], myEM$sigma[2])), col="red", lwd=3, add=TRUE)
  lines(c(1.64,1.64), c(0,1), col="blue", lwd=3, lty=2)

  # get classification
  p <- myEM$posterior
  if(mean(dat[p[,1]>p[,2]]) > mean(dat[p[,1]<p[,2]])){
    m[names(dat[p[,1]>p[,2]]),i]=1
  }else{
    m[names(dat[p[,1]<p[,2]]),i]=1
  }
}
dev.off()

write.table(m, file = outMatrix, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
