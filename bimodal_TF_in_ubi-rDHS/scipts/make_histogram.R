
# --Kaili
# This script is for geting histogram from given data matrix.
# INPUT: data matrix (signal/zscore), giving a spcific type (log2, log10,zscore)
# OUTPUT: histogram in given Dir

# EXP: Rscript make_histogram.R /data/zusers/fankaili/ccre/tf/matrix/hg19_ubi-rDHS_CTCF_signal_matrix.txt /data/zusers/fankaili/ccre/tf/figs/ "hg19_ubi-rDHS_CTCF_signal_log10_histogram.pdf" log10

args<-commandArgs(T)
inFile = args[1]
outDir = args[2]
outFile = args[3]
type = args[4] ## log2, log10, zscore

setwd(outDir)
data <- read.table(inFile, header=TRUE, row.names = 1)


pdf(outFile)
for(i in 1:ncol(data)){
  name = colnames(data)[i]
  if(type=="log2"){
    dat = log2(data[,i]+0.01)
  } else if (type=="log10"){
    dat = log10(data[,i]+0.01)
  } else if(type=="zscore"){
    dat = data[,i]
  }
  m = paste("CTCF signal of ubi-rDHS in ",name,"\t",type ,sep="")
  hist(dat, main = name, xlab = name, col = "#bdbdbd", freq = FALSE,
       breaks=seq(min(dat)-0.05, max(dat)+0.05, by=0.05))
  if(type=="zscore"){
    lines(c(1.64,1.64), c(0,20), col="blue", lwd=2)
  }
}

dev.off()
