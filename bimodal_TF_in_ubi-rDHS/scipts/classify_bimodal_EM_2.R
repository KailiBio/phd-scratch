
# -- Kaili
# This script is for classify bimodal distribution (GMM) data by EM algorithm
# This one for H3K4me3

outDir = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/"
tf="H3K4me3"


inFile = "hg19_ubi-rDHS_H3K4me3_zscore_matrix.txt"
outMatrix = paste("hg19_ubi-rDHS_",tf,"_zscore_classification.txt", sep="")
outFigure = paste("hg19_ubi-rDHS_",tf,"_zscore_classification.pdf", sep="")

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
  # if(type=="log2"){
  #   dat = log2(data[,i]+0.01)
  # } else if (type=="log10"){
  #   dat = log10(data[,i]+0.01)
  # } else if(type=="zscore"){
  #   dat = data[data[,i]!=(-10),i]
  #   # remove outliner (-10) here
  # }
  # names(dat) <- rownames(data[data[,i]!=(-10),])
  dat_z = data[,i]
  dat = data[data[,i]!=(-10),i]
  names(dat) <- rownames(data[data[,i]!=(-10),])

  # call EM
  if(tf=="CTCF" || tf=="CTCF_2"){
    myEM <- normalmixEM(dat, mu = c(min(dat),max(dat)), sigma = c(2,1))
  }else{
    if(i==11 | i==17){
      dat0=dat
      dat=dat0[dat0<6]
      myEM <- normalmixEM(dat, mu = c(min(dat),max(dat)), sigma = c(4,1))
    }else if(i==15){
      dat0=dat
      dat=dat0[dat0<5]
      myEM <- normalmixEM(dat, mu = c(min(dat),max(dat)), sigma = c(1,10))
    }else{
      myEM <- normalmixEM(dat, mu = c(min(dat),max(dat)), sigma = c(4,1))
    }
  }

  # get classification
  p <- myEM$posterior
  if(mean(dat[p[,1]>p[,2]]) > mean(dat[p[,1]<p[,2]])){
    print("group1")
    m[names(dat[p[,1]>p[,2]]),i]=1
    t=(min(dat[p[,1]>p[,2]])+max(dat[p[,1]<p[,2]]))/2
  }else{
    print("group2")
    m[names(dat[p[,1]<p[,2]]),i]=1
    t=(min(dat[p[,1]<p[,2]])+max(dat[p[,1]>p[,2]]))/2
  }
  #t=mean(max(min(dat[p[,1]>p[,2]]),min(dat[p[,1]<p[,2]])), min(max(dat[p[,1]>p[,2]]),max(dat[p[,1]<p[,2]])))

  # get num
  n = length(dat_z)
  n_1.64 = length(dat_z[dat_z>1.64])
  if(i==11 | i==15 | i==17){
    n_em = length(dat0[dat0>=t])
  }else{
    n_em = length(dat[dat>=t])
  }
  l1 = paste("n = ", as.character(n), sep="")
  l2 = paste("n(z-score>1.64) = ", as.character(n_1.64), sep="")
  l3 = paste("n(z-score>",as.character(round(t,2)),";EM) = ", as.character(n_em), sep="")

  # density plot
  if(i==11 | i==15 | i==17){
    hist(dat0, breaks=200, freq = FALSE, main=paste("H3K4me3:",biosample, sep=" "), xlab="zscore", col="grey")
  }else{
    hist(dat, breaks=200, freq = FALSE, main=paste("H3K4me3:",biosample, sep=" "), xlab="zscore", col="grey")
  }
  curve((myEM$lambda[1]*dnorm(x, myEM$mu[1], myEM$sigma[1])), col="green", lwd=3, add=TRUE)
  curve((myEM$lambda[2]*dnorm(x, myEM$mu[2], myEM$sigma[2])), col="yellow", lwd=3, add=TRUE)
  lines(c(1.64,1.64), c(0,2), col="blue", lwd=3, lty=2)
  lines(c(t,t), c(0,2), col="red", lwd=3, lty=2)
  legend("topleft", legend=paste("ubi-rDHS:",l1,l2,l3,sep="\n"), bty = "n")
}
dev.off()

write.table(m, file = outMatrix, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
