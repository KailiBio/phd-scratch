
# -- Kaili
# This script is for making scatter plot of ubi-rDHS CTCF zscore.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/ss_zscore/")

pls <- as.vector(read.table("../hg19_ubi-rDHS_PLS_list_DHSID.txt")[,1])
no_pls <- as.vector(read.table("../hg19_ubi-rDHS_non_PLS_list_DHSID.txt")[,1])

file=list.files()

#######################
# non PLS
#######################
pdf("../ubi-rDHS_non_PLS_CTCF_scatter.pdf", height = 8, width = 8)
for (i in 1:length(file)){
  a = data.frame(read.table(file[i]))
  rownames(a) = a[,1]
  a_name = unlist(strsplit(file[i], "_"))[3]
  for (j in 1:length(file)){
    if(j>i){
      b = data.frame(read.table(file[j]))
      rownames(b) = b[,1]
      b_name = unlist(strsplit(file[j], "_"))[3]
      #
      x = a[no_pls,]
      y = b[no_pls,]
      xx = x[x$V2>1.64,]
      yy = y[y$V2>1.64,]
      xx1 = nrow(xx)
      yy1 = nrow(yy)
      m=length(intersect(rownames(xx), rownames(yy)))
      n=length(union(rownames(xx), rownames(yy)))
      jaccard = round(m/n,3)
      #
      smoothScatter(x$V2, y$V2, nrpoints = 0, main = "non-PLS (n = 1,912)\n",
                    xlab = paste(a_name, " (z-score>1.64: ", xx1, ")", sep=""),
                    ylab = paste(b_name, " (z-score>1.64: ", yy1, ")", sep=""))
      lines(c(1.64, 1.64), c(-12,10), col="red", lwd=3, lty=2)
      lines(c(-12,10), c(1.64, 1.64), col="red", lwd=3, lty=2)
      legend("topleft", legend = paste("intersect:",m,"\nunion:",n,"\nJaccardIndex=",jaccard),
             bty = "n")
    }
  }
}
dev.off()
#######################
# PLS
#######################
pdf("../ubi-rDHS_PLS_CTCF_scatter.pdf", height = 8, width = 8)
for (i in 1:length(file)){
  a = data.frame(read.table(file[i]))
  rownames(a) = a[,1]
  a_name = unlist(strsplit(file[i], "_"))[3]
  for (j in 1:length(file)){
    if(j>i){
      b = data.frame(read.table(file[j]))
      rownames(b) = b[,1]
      b_name = unlist(strsplit(file[j], "_"))[3]
      #
      x = a[pls,]
      y = b[pls,]
      xx = x[x$V2>1.64,]
      yy = y[y$V2>1.64,]
      xx1 = nrow(xx)
      yy1 = nrow(yy)
      m=length(intersect(rownames(xx), rownames(yy)))
      n=length(union(rownames(xx), rownames(yy)))
      jaccard = round(m/n,3)
      #
      smoothScatter(x$V2, y$V2, nrpoints = 0, main = "ubi-rDHS non-PLS (n = 9,009)\n",
                    xlab = paste(a_name, " (z-score>1.64: ", xx1, ")", sep=""),
                    ylab = paste(b_name, " (z-score>1.64: ", yy1, ")", sep=""))
      lines(c(1.64, 1.64), c(-12,10), col="red", lwd=3, lty=2)
      lines(c(-12,10), c(1.64, 1.64), col="red", lwd=3, lty=2)
      legend("topleft", legend = paste("intersect:",m,"\nunion:",n,"\nJaccardIndex=",jaccard),
             bty = "n")
    }
  }
}
dev.off()


#######################
# ubi-rDHS
#######################
pdf("../ubi-rDHS_CTCF_scatter.pdf", height = 8, width = 8)
for (i in 1:length(file)){
  a = data.frame(read.table(file[i]))
  rownames(a) = a[,1]
  a_name = unlist(strsplit(file[i], "_"))[3]
  for (j in 1:length(file)){
    if(j>i){
      b = data.frame(read.table(file[j]))
      rownames(b) = b[,1]
      b_name = unlist(strsplit(file[j], "_"))[3]
      #
      x = a
      y = b
      xx = x[x$V2>1.64,]
      yy = y[y$V2>1.64,]
      xx1 = nrow(xx)
      yy1 = nrow(yy)
      m=length(intersect(rownames(xx), rownames(yy)))
      n=length(union(rownames(xx), rownames(yy)))
      jaccard = round(m/n,3)
      #
      smoothScatter(x$V2, y$V2, nrpoints = 0, main = "ubi-rDHS (n = 10,921)\n",
                    xlab = paste(a_name, " (z-score>1.64: ", xx1, ")", sep=""),
                    ylab = paste(b_name, " (z-score>1.64: ", yy1, ")", sep=""))
      lines(c(1.64, 1.64), c(-12,10), col="red", lwd=3, lty=2)
      lines(c(-12,10), c(1.64, 1.64), col="red", lwd=3, lty=2)
      legend("topleft", legend = paste("intersect:",m,"\nunion:",n,"\nJaccardIndex=",jaccard),
             bty = "n")
    }
  }
}
dev.off()


#######################
