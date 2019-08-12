
# -- Kaili
# This script is for doing clustering for M4.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/mm10_RNAseq/")

library(ggplot2)
library(dendextend)
library(devtools)
library(reshape2)
library(factoextra)

#############
# function
#############

make_dendrogram_euclidean <- function(data, title){
  dat = data.frame(t(data))
  dat$sample = colnames(data)
  dat$tissue = sapply(dat$sample, function(x) strsplit(as.character(x),"_")[[1]][1])
  dat$time = sapply(dat$sample, function(x) strsplit(as.character(x),"_")[[1]][2])
  
  n = ncol(dat)
  s = as.integer(n-2)
  t = as.integer(n-1)
  p = as.integer(n)
  
  dd <- dist(dat[,1:(n-3)], diag=TRUE)
  #dd <- as.dist(1-cor(t(dat[,c(1:39924)])))
  hc <- hclust(dd)
  hcd <- as.dendrogram(hc)
  
  # specific_leaf=hcd[[2]][[1]][[1]][[1]]
  # specific_leaf
  # attributes(specific_leaf)
  
  i=0
  colLab<<-function(n){
    if(is.leaf(n)){
      
      #I take the current attributes
      a=attributes(n)
      
      #I deduce the line in the original data, and so the treatment and the specie.
      ligne=match(attributes(n)$label,dat[,s])
      tissue=dat[ligne,t];
      if(tissue=="forebrain"){col="#008000"}
      if(tissue=="midbrain"){col="#06DA93"}
      if(tissue=="hindbrain"){col="#00aa00"}
      if(tissue=="neural.tube"){col="#00bfc4"}
      if(tissue=="heart"){col="#FF0000"}
      if(tissue=="facial"){col="#1262EB"}
      if(tissue=="limb"){col="#969696"}
      if(tissue=="kidney"){col="#A872E5"}
      if(tissue=="lung"){col="#6a3d9a"}
      if(tissue=="stomach"){col="#D642CA"}
      if(tissue=="intestine"){col="#b15928"}
      if(tissue=="liver"){col="#FBB30B"}
      time = dat[ligne,p];
      if(time=="10.5"){shape=0}
      if(time=="11.5"){shape=1}
      if(time=="12.5"){shape=2}
      if(time=="13.5"){shape=3}
      if(time=="14.5"){shape=15}
      if(time=="15.5"){shape=20}
      if(time=="16.5"){shape=17}
      if(time=="0"){shape=18}
      
      
      #Modification of leaf attribute
      attr(n,"nodePar")<-c(a$nodePar,list(cex=1,lab.cex=0.8,pch=shape,col=col,lab.col=col,
                                          lab.font=0.8, lab.font=1))
    }
    return(n)
  }
  
  hc_col <- dendrapply(hcd, colLab)
  
  pdf(paste("./clustering_M4/",title,".pdf",sep=""), height=15)
  par(mai=c(0.5,0.5,0.5,1.5))
  plot(hc_col , main=title, horiz=TRUE)
  dev.off()
}

make_dendrogram_pearson <- function(data, title){
  dat = data.frame(t(data))
  dat$sample = colnames(data)
  dat$tissue = sapply(dat$sample, function(x) strsplit(as.character(x),"_")[[1]][1])
  dat$time = sapply(dat$sample, function(x) strsplit(as.character(x),"_")[[1]][2])
  
  n = ncol(dat)
  s = as.integer(n-2)
  t = as.integer(n-1)
  p = as.integer(n)
  dd <- as.dist(1-cor(t(dat[,1:(n-3)])))
  hc <- hclust(dd)
  hcd <- as.dendrogram(hc)
  
  
  # specific_leaf=hcd[[2]][[1]][[1]][[1]]
  # specific_leaf
  # attributes(specific_leaf)
  
  i=0
  colLab<<-function(n){
    if(is.leaf(n)){
      
      #I take the current attributes
      a=attributes(n)
      
      #I deduce the line in the original data, and so the treatment and the specie.
      ligne=match(attributes(n)$label,dat[,s])
      tissue=dat[ligne,t];
      if(tissue=="forebrain"){col="#008000"}
      if(tissue=="midbrain"){col="#06DA93"}
      if(tissue=="hindbrain"){col="#00aa00"}
      if(tissue=="neural.tube"){col="#00bfc4"}
      if(tissue=="heart"){col="#FF0000"}
      if(tissue=="facial"){col="#1262EB"}
      if(tissue=="embryonic.facial.prominence"){col="#1262EB"}
      if(tissue=="limb"){col="#969696"}
      if(tissue=="kidney"){col="#A872E5"}
      if(tissue=="lung"){col="#6a3d9a"}
      if(tissue=="stomach"){col="#D642CA"}
      if(tissue=="intestine"){col="#b15928"}
      if(tissue=="liver"){col="#FBB30B"}
      time = dat[ligne,p];
      if(time=="10.5"){shape=0}
      if(time=="11.5"){shape=1}
      if(time=="12.5"){shape=2}
      if(time=="13.5"){shape=3}
      if(time=="14.5"){shape=15}
      if(time=="15.5"){shape=20}
      if(time=="16.5"){shape=17}
      if(time=="0"){shape=18}
      
      
      #Modification of leaf attribute
      attr(n,"nodePar")<-c(a$nodePar,list(cex=1,lab.cex=0.8,pch=shape,col=col,lab.col=col,
                                          lab.font=0.8, lab.font=1))
    }
    return(n)
  }
  
  hc_col <- dendrapply(hcd, colLab)
  
  pdf(paste("./clustering_M4/",title,".pdf",sep=""), height=15)
  par(mai=c(0.5,0.5,0.5,1.5))
  plot(hc_col , main=title, horiz=TRUE)
  #hang = -0.1
  dev.off()
}

make_dendrogram_spearman <- function(data, title){
  dat = data.frame(t(data))
  dat$sample = colnames(data)
  dat$tissue = sapply(dat$sample, function(x) strsplit(as.character(x),"_")[[1]][1])
  dat$time = sapply(dat$sample, function(x) strsplit(as.character(x),"_")[[1]][2])
  
  n = ncol(dat)
  s = as.integer(n-2)
  t = as.integer(n-1)
  p = as.integer(n)
  dd <- as.dist(1-cor(t(dat[,1:(n-3)]), method="spearman"))
  hc <- hclust(dd)
  hcd <- as.dendrogram(hc)
  
  
  # specific_leaf=hcd[[2]][[1]][[1]][[1]]
  # specific_leaf
  # attributes(specific_leaf)
  
  i=0
  colLab<<-function(n){
    if(is.leaf(n)){
      
      #I take the current attributes
      a=attributes(n)
      
      #I deduce the line in the original data, and so the treatment and the specie.
      ligne=match(attributes(n)$label,dat[,s])
      tissue=dat[ligne,t];
      if(tissue=="forebrain"){col="#008000"}
      if(tissue=="midbrain"){col="#06DA93"}
      if(tissue=="hindbrain"){col="#00aa00"}
      if(tissue=="neural.tube"){col="#00bfc4"}
      if(tissue=="heart"){col="#FF0000"}
      if(tissue=="facial"){col="#1262EB"}
      if(tissue=="limb"){col="#969696"}
      if(tissue=="kidney"){col="#A872E5"}
      if(tissue=="lung"){col="#6a3d9a"}
      if(tissue=="stomach"){col="#D642CA"}
      if(tissue=="intestine"){col="#b15928"}
      if(tissue=="liver"){col="#FBB30B"}
      time = dat[ligne,p];
      if(time=="10.5"){shape=0}
      if(time=="11.5"){shape=1}
      if(time=="12.5"){shape=2}
      if(time=="13.5"){shape=3}
      if(time=="14.5"){shape=15}
      if(time=="15.5"){shape=20}
      if(time=="16.5"){shape=17}
      if(time=="0"){shape=18}
      
      
      #Modification of leaf attribute
      attr(n,"nodePar")<-c(a$nodePar,list(cex=1,lab.cex=0.8,pch=shape,col=col,lab.col=col,
                                          lab.font=0.8, lab.font=1))
    }
    return(n)
  }
  
  hc_col <- dendrapply(hcd, colLab)
  
  pdf(paste("./clustering_M4/",title,".pdf",sep=""), height=15)
  par(mai=c(0.5,0.5,0.5,1.5))
  plot(hc_col , main=title, horiz=TRUE)
  dev.off()
}

make_PCA_scatter <- function(data, title){
  dat = t(data)
  dat.pca <- prcomp(dat)
  
  dat_pca_m = data.frame(dat.pca$x[,1:2])
  dat_pca_m$sample = substr(rownames(dat_pca_m),1,nchar(rownames(dat_pca_m))-2)
  dat_pca_m$tissue = sapply( strsplit(dat_pca_m$sample,"_"), "[", 1)
  dat_pca_m$time = sapply( strsplit(dat_pca_m$sample,"_"), "[", 2)
  
  ggplot(dat_pca_m, aes(x=PC1, y=PC2, shape=time, col=tissue)) + geom_point(size=3.5) +
    scale_shape_manual(values=c(18,0,1,2,3,15,20,17)) +
    scale_color_manual(values = c("#1262EB", "#008000", "#FF0000", "#00aa00", "#b15928", "#A872E5",
                                  "#969696", "#FBB30B", "#6a3d9a", "#06DA93", "#00bfc4", "#D642CA")) +
    theme_minimal() + labs(title = title)
  ggsave(paste(title,".pdf",sep=""), path="./PCA/")
  ggsave(paste(title,".png",sep=""), path="./PCA/")
}


library(Rtsne)
make_tsne_scatter <- function(data, title){
  set.seed(1)  
  tsne_model = Rtsne(t(data), perplexity=30, dims=2, check_duplicates = FALSE)
  
  d_tsne = as.data.frame(tsne_model$Y)
  d_tsne$sample = substr(colnames(data),1,nchar(colnames(data))-2)
  d_tsne$tissue = sapply( strsplit(d_tsne$sample,"_"), "[", 1)
  d_tsne$time = sapply( strsplit(d_tsne$sample,"_"), "[", 2)
  
  ggplot(d_tsne, aes(x=V1, y=V2, shape=time, col=tissue)) +  
    geom_point(size=3.5) +
    scale_shape_manual(values=c(18,0,1,2,3,15,20,17)) +
    scale_color_manual(values = c("#1262EB", "#008000", "#FF0000", "#00aa00", "#b15928", "#A872E5",
                                  "#969696", "#FBB30B", "#6a3d9a", "#06DA93", "#00bfc4", "#D642CA")) +
    theme_minimal() + labs(title = title)
  ggsave(paste(title,".pdf",sep=""), path="./clustering_M4/")
  ggsave(paste(title,".png",sep=""), path="./clustering_M4/")
}


#############
# read data
#############
data0 = read.table("mm10_M4_TPM_matrix.txt", row.names = 1, header = 1)
data_raw0 = data0[ , order(colnames(data0))]
dat_raw = data_raw0[unlist(apply(data_raw0,1,function(x) var(x)!=0)),]
# 35228

# log only
dat_raw_log = log10(dat_raw+0.01)

# limma
library(limma)
batch = sapply(colnames(dat_raw), function(x) strsplit(as.character(x),"_")[[1]][2])
#
dat_raw_limma <- removeBatchEffect(dat_raw, batch)
#
dat_log_limma <- removeBatchEffect(dat_raw_log, batch)
#write.table(dat_log_limma, "mm10_exp_TPM-log-limma_matrix.txt", sep="\t", quote = FALSE, 
            row.names = TRUE, col.names = TRUE)

# batch2 = sapply(colnames(dat_raw), function(x) strsplit(as.character(x),"_")[[1]][1])
# dat_log_limma2 <- removeBatchEffect(dat_raw_log, batch2)

# Combat
library(sva)
batch = as.integer(sapply(colnames(dat_raw), function(x) strsplit(as.character(x),"_")[[1]][2]))
#
dat_raw_combat <- ComBat(dat=as.matrix(dat_raw), batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
#
dat_log_combat <- ComBat(dat=as.matrix(dat_raw_log), batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
#write.table(dat_log_combat, "mm10_exp_TPM-log-combat_matrix.txt", sep="\t", quote = FALSE, 
            row.names = TRUE, col.names = TRUE)

#############
# clutering
#############

make_dendrogram_euclidean(dat_raw,"Eclustering_TPM")
make_dendrogram_euclidean(dat_raw_log,"Eclustering_TPM_log")
make_dendrogram_euclidean(dat_raw_limma,"Eclustering_TPM_limma")
make_dendrogram_euclidean(dat_log_limma,"Eclustering_TPM_log_limma")
make_dendrogram_euclidean(dat_raw_combat,"Eclustering_TPM_raw_combat")
make_dendrogram_euclidean(dat_log_combat,"Eclustering_TPM_log_combat")

make_dendrogram_pearson(dat_raw,"Pclustering_TPM")
make_dendrogram_pearson(dat_raw_log,"Pclustering_TPM_log")
make_dendrogram_pearson(dat_raw_limma,"Pclustering_TPM_limma")
make_dendrogram_pearson(dat_log_limma,"Pclustering_TPM_log_limma")
make_dendrogram_pearson(dat_raw_combat,"Pclustering_TPM_raw_combat")
make_dendrogram_pearson(dat_log_combat,"Pclustering_TPM_log_combat")


make_dendrogram_spearman(dat_raw,"Sclustering_TPM")
make_dendrogram_spearman(dat_raw_log,"Sclustering_TPM_log")
make_dendrogram_spearman(dat_raw_limma,"Sclustering_TPM_limma")
make_dendrogram_spearman(dat_log_limma,"Sclustering_TPM_log_limma")
make_dendrogram_spearman(dat_raw_combat,"Sclustering_TPM_raw_combat")
make_dendrogram_spearman(dat_log_combat,"Sclustering_TPM_log_combat")


make_tsne_scatter(dat_raw,"mm10_TPM_tSNE")
make_tsne_scatter(dat_raw_log,"mm10_TPM_log_tSNE")
make_tsne_scatter(dat_raw_limma,"mm10_TPM_limma_tSNE")
make_tsne_scatter(dat_log_limma,"mm10_TPM_log_limma_tSNE")
make_tsne_scatter(dat_raw_combat,"mm10_TPM_combat_tSNE")
make_tsne_scatter(dat_log_combat,"mm10_TPM_log_combat_tSNE")


