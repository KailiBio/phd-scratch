
# -- Kaili
# This script is for clustering states from .para0 file.

args<-commandArgs(trailingOnly=TRUE)
file = args[1]
workDir = args[2]

file="ctcf_9sample_impute.para0"
#col = "ctcf_9sample_impute.col.txt"
workDir = "/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/"

setwd(workDir)

library(pheatmap)

#####
# This function is for reading a parafile
read_para_file <- function(file){
  x=read.table(file, comment="!", header=T);
  k=dim(x)[2];
  l=dim(x)[1];
  p=(sqrt(9+8*(k-1))-3)/2;
  m=as.matrix(x[,1+1:p]/x[,1]);
  colnames(m) = colnames(x)[1+1:p];
  rownames(m) = c(0:(l-1))
  m_sort = m[,order(colnames(m))]
  return(m_sort)
}

#################

dat = read_para_file(file)

pdf("state_cluster.pdf", height = 5, width = 8)
pheatmap(t(dat), cluster_cols=T, cluster_rows = F, cellheight = 12, cellwidth=10, 
         clustering_distance_cols = "euclidean")
dev.off()

###
clusters = hclust(dist(dat))
plot(clusters)
clustercut <- cutree(clusters,2)
list2 = names(clustercut)[clustercut==1]
#
clusters2 = hclust(dist(dat[list2,]))
plot(clusters2)
clustercut2 <- cutree(clusters2,3)

###
state_type = transform(cbind(state = c(0:46), type = "other"))
rownames(state_type) = 0:46
state_type$type = as.vector(state_type$type)
state_type[names(clustercut2)[clustercut2==2],]$type = "CTCF"
state_type[names(clustercut2)[clustercut2==3],]$type = "enhancer"
state_type[names(clustercut)[clustercut==2],]$type = "promoter"
write.table(state_type, "state_type_mark.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)

#################
