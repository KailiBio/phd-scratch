
# -- Kaili
# This script is for clustering merged bins using signal.

# setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/pool-bins/") 

library(data.table)
library(dendextend)

args <- commandArgs(trailingOnly=TRUE);
state = args[1]
mark = args[2]

# mark_list = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K9me3")
setwd("/data/zusers/fankaili/ideas/dhs_ctcf/pool-bins/")

file = paste("./state",state,"_",mark,"/state",state,"_",mark,"_signalMatrix.txt",sep="")
data0 = fread(file)
data = data0[,2:67]
dd <- dist(t(data), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
hcd <- as.dendrogram(hc)
hcd1 <- color_labels(hcd, k = 5)
#
outfile = paste("./hclust/state",state,"_",mark,"_hlcust.pdf",sep="")
name = paste("state",state," ",mark, sep="")
pdf(outfile, width = 6, height = 12)
par(mai=c(0.5,0.5,0.2,3))
plot(hcd1,  xlab = "Height",horiz = TRUE, cex=0.1, main = name)
dev.off()
    
