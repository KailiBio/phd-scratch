
# -- Kaili
# This script is for analyzing mm10 RNA-seq batch effect and make figures.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/mm10_RNAseq/")

library(ggplot2)
library(dendextend)
library(devtools)
library(reshape2)
library(factoextra)

#############
# function
#############
# # perform clustering (on transposed mtx :-/)
# clustfun <- function(x) {
#   # http://dx.doi.org/10.1016/j.csda.2005.10.006
#   #d <- as.dist(sqrt(1-cor(t(x))))
#   d <- as.dist(1-cor(t(x)))
#   p <- as.phylo(hclust(d))
#   return(p)
# }

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
  
  pdf(paste("./spikeIn96_Jul28/clustering_euclidean/",title,".pdf",sep=""), height=15)
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
  
  pdf(paste("./spikeIn96_Jul28/clustering_pearson/",title,".pdf",sep=""), height=15)
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
  
  pdf(paste("./spikeIn96_Jul28/clustering_spearman/",title,".pdf",sep=""), height=15)
  par(mai=c(0.5,0.5,0.5,1.5))
  plot(hc_col , main=title, horiz=TRUE)
  dev.off()
}

# library(ggdendro)
# dendr    <- dendro_data(hc, type="rectangle") # convert for ggplot
# 
# ggplot() + 
#   geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
#   geom_text(data=label(dendr), aes(x, y, label=label, hjust=0), 
#             size=2.5) +
#   coord_flip() + scale_y_reverse(expand=c(0.5, 0)) + 
#   theme(axis.line.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         panel.background=element_rect(fill="white"),
#         panel.grid=element_blank())
# ggsave("mm10_raw_hcluster.pdf", height=12)

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
  ggsave(paste(title,".pdf",sep=""), path="./spikeIn96_Jul28/tSNE/")
  ggsave(paste(title,".png",sep=""), path="./spikeIn96_Jul28/tSNE/")
}



#############
# read data
#############
# get sample_num and color for each tissue
tissue_num = data.frame(rbind(c(16,16,16,12,16,12,12,8,8,8,8,14),
                              c("#b3de69","#8dd3c7", "#ccebc5", "#FBB30B","#bebada", "#fb8072", "#80b1d3", 
                                "#fdb462", "#fccde5", "#d9d9d9", "#bc80bd", "#ffed6f")))
colnames(tissue_num) = c("forebrain", "midbrain", "hindbrain", "neural.tube", "heart", "facial", "limb", 
                         "kidney", "lung", "stomach", "interstine", "liver")
a = data.frame(tissue_num[,order(colnames(tissue_num)),])
colnames(a) = colnames(tissue_num)[order(colnames(tissue_num))]
rownames(a) = c("num","col")
tissue_num = a

###############
# all data
###############

# raw data matrix
data0 = read.table("mm10_exp_matrix.txt", row.names = 1, header = 1)
data_raw0 = data0[ , order(colnames(data0))]
dat_raw = data_raw0[unlist(apply(data_raw0,1,function(x) var(x)!=0)),]
# 39924
# 
# write.table(data_raw0, "mm10_exp_matrix_sorted.txt", sep="\t", quote = FALSE, 
#             row.names = TRUE, col.names = TRUE)


a = data0[,"neural.tube_14.5_2"]
b = dat0[,"neural.tube_14.5_2"]
cor(a,b, method="pearson")
cor(a,b, method="spearman")
smoothScatter(a, b)

cor(data0[,"neural.tube_14.5_2"],data0[,"hindbrain_14.5_2"], method="pearson")
cor(dat0[,"neural.tube_14.5_2"],dat0[,"hindbrain_14.5_2"], method="pearson")


cor(data0[,"hindbrain_14.5_1"],data0[,"hindbrain_14.5_2"], method="pearson")
cor(dat0[,"hindbrain_14.5_1"],dat0[,"hindbrain_14.5_2"], method="pearson")

dat0 = read.table("mm10_exp_matrix_spikeIn96_Jul28.txt", row.names = 1, header = 1)
dat_raw0 = dat0[ , order(colnames(dat0))]
dat_raw = dat_raw0[unlist(apply(dat_raw0,1,function(x) var(x)!=0)),]
# 39727
# spikeIn96_Jul28: 39906



# dat_1 = dat_raw[,seq(1,145,2)]
# dat_2 = dat_raw[,seq(2,146,2)]
# dat_avg = dat_1 + dat_2

# z-score matrix
dat1 = read.table("mm10_exp_zscore_matrix.txt", row.names = 1, header = 1)
dat_zscore0 = dat1[, order(colnames(dat1))]
dat_zscore = dat_zscore0[unlist(apply(dat_zscore0,1,function(x) var(x)>0)),]
# 39924

# log only
dat_raw_log = log10(dat_raw+0.01)

# quantile normalization
library(preprocessCore)
dat_raw_quan = normalize.quantiles(as.matrix(dat_raw))
colnames(dat_raw_quan) = colnames(dat_raw)
rownames(dat_raw_quan) = rownames(dat_raw)
#
dat_zscore_quan = normalize.quantiles(as.matrix(dat_zscore))
colnames(dat_zscore_quan) = colnames(dat_zscore)
rownames(dat_zscore_quan) = rownames(dat_zscore)
#
dat_log_quan = normalize.quantiles(as.matrix(dat_raw_log))
colnames(dat_log_quan) = colnames(dat_raw_log)
rownames(dat_log_quan) = rownames(dat_raw_log)

# limma
library(limma)
batch = sapply(colnames(dat_raw), function(x) strsplit(as.character(x),"_")[[1]][2])
#
dat_raw_limma <- removeBatchEffect(dat_raw, batch)
#
dat_log_limma <- removeBatchEffect(dat_raw_log, batch)
write.table(dat_log_limma, "mm10_exp_TPM-log-limma_matrix.txt", sep="\t", quote = FALSE, 
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
write.table(dat_log_combat, "mm10_exp_TPM-log-combat_matrix.txt", sep="\t", quote = FALSE, 
            row.names = TRUE, col.names = TRUE)

# RUVseq
library(RUVSeq)
#
count_matrix0 = read.table("mm10_exp_counts_matrix.txt", row.names=1, header=1)
count_matrix1 = count_matrix0[rownames(dat_raw),]
count_matrix = count_matrix1[, order(colnames(count_matrix1))]
#
spikes_matrix0 = read.table("mm10_spikeIn_exp_matrix.txt", row.names = 1, header = 1)
spikes_matrix = spikes_matrix0[, order(colnames(spikes_matrix0))]
spikes = rownames(spikes_matrix)
#
matrix = rbind(count_matrix, spikes_matrix)
set=newSeqExpressionSet(as.matrix(matrix))
dat_raw_RUVseq=normCounts(RUVg(set,spikes,k=1))



# edgeR
# library(edgeR)
# #
# counts <- DGEList(counts=dat_raw,group=batch)
# norFactor <- calcNormFactors(counts)
# dat_raw_edgeR <- t(t(norFactor$counts)/norFactor$samples[,"norm.factors"])
# #
# counts <- DGEList(counts=dat_raw_log, group=batch)
# norFactor <- calcNormFactors(counts)
# dat_log_edgeR <- t(t(norFactor$counts)/norFactor$samples[,"norm.factors"])


###############
# no e10.5
###############
rm_list=c("forebrain_10.5_1", "forebrain_10.5_2", "midbrain_10.5_1", "midbrain_10.5_2", "hindbrain_10.5_1",
          "hindbrain_10.5_2", "heart_10.5_1", "heart_10.5_2", "facial_10.5_1", "facial_10.5_2",
          "limb_10.5_1", "limb_10.5_2", "neural.tube_0_1", "neural.tube_0_2")
dat_noe10.5 = dat_raw[,!colnames(dat_raw)%in%rm_list]
dat2_1 = dat_noe10.5[,seq(1,132,2)]
dat2_2 = dat_noe10.5[,seq(2,132,2)]
dat2_avg = (dat2_1 + dat2_2)/2
a = unlist(lapply(colnames(dat2_avg), function(x) substr(x,start=1,stop=nchar(x)-2)))
colnames(dat2_avg) = a

# TPM
dat2_no0var_TPM = dat2_avg[unlist(apply(dat2_avg,1,function(x) var(x)>0)),]

# TPM_log
dat2_TPM_log = log10(dat2_avg+0.01)
dat2_no0var_TPM_log = dat2_TPM_log[unlist(apply(dat2_TPM_log,1,function(x) var(x)>0)),]

# limma
library(limma)
batch = sapply(colnames(dat2_avg), function(x) strsplit(as.character(x),"_")[[1]][2])
#
dat_TPM_limma <- removeBatchEffect(dat2_avg, batch)
dat2_no0var_TPM_limma = dat_TPM_limma[unlist(apply(dat_TPM_limma,1,function(x) var(x)>0)),]
#
dat_TPM_log_limma <- removeBatchEffect(dat2_TPM_log, batch)
dat2_no0var_TPM_log_limma = dat_TPM_log_limma[unlist(apply(dat_TPM_log_limma,1,function(x) var(x)>0)),]



#############
# make violin plot
#############
# raw
dat = melt(dat_raw, value.name = "value" )
dat$tissue = sapply(dat$variable, function(x) strsplit(as.character(x),"_")[[1]][1])
ggplot(dat, aes(x=variable, y=value)) + geom_violin() + facet_grid(tissue ~ .) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("violin_all_raw_signal.pdf")

# z-score
dat = melt(dat_zscore, value.name = "value" )
dat$tissue = sapply(dat$variable, function(x) strsplit(as.character(x),"_")[[1]][1])
ggplot(dat, aes(x=variable, y=value)) + geom_violin() + facet_grid(tissue ~ .) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("violin_all_zscore_signal.pdf")

#############
# distribution
#############
## z-score
pdf("histogram_zscore.pdf")
for(i in 1:146){
  dat = data.frame(dat_zscore[,1])
  colnames(dat) = c("value")
  p = ggplot(dat, aes(x=value)) + geom_histogram(binwidth = 1)
  print(p)
}
dev.off()

#############
# clutering
#############
make_dendrogram_euclidean(dat_raw,"Eclustering_TPM")
make_dendrogram_euclidean(dat_raw_log,"Eclustering_TPM_log")
make_dendrogram_euclidean(dat_zscore,"Eclustering_TPM_zscore")
make_dendrogram_euclidean(dat_raw_quan,"Eclustering_TPM_quantile")
make_dendrogram_euclidean(dat_log_quan,"Eclustering_TPM_log_quantile")
make_dendrogram_euclidean(dat_zscore_quan,"Eclustering_TPM_zscore_quantile")
make_dendrogram_euclidean(dat_raw_limma,"Eclustering_TPM_limma")
make_dendrogram_euclidean(dat_log_limma,"Eclustering_TPM_log_limma")
make_dendrogram_euclidean(dat_raw_combat,"Eclustering_TPM_raw_combat")
make_dendrogram_euclidean(dat_log_combat,"Eclustering_TPM_log_combat")
make_dendrogram_euclidean(dat_raw_RUVseq,"Eclustering_counts_ruvseq")
# make_dendrogram_euclidean(dat_raw_edgeR,"Eclustering_TPM_raw_edger")
# make_dendrogram_euclidean(dat_log_edgeR,"Eclustering_TPM_log_edger")

make_dendrogram_pearson(dat_raw,"Pclustering_TPM2")
make_dendrogram_pearson(dat_raw_log,"Pclustering_TPM_log2")
make_dendrogram_pearson(dat_raw_limma,"Pclustering_TPM_limma2")
make_dendrogram_pearson(dat_log_limma,"Pclustering_TPM_log_limma2")
make_dendrogram_pearson(dat_raw_combat,"Pclustering_TPM_raw_combat2")
make_dendrogram_pearson(dat_log_combat,"Pclustering_TPM_log_combat2")
make_dendrogram_pearson(dat_raw_RUVseq,"Pclustering_counts_ruvseq2")
# make_dendrogram_pearson(dat_raw_edgeR,"Pclustering_TPM_raw_edger")
# make_dendrogram_pearson(dat_log_edgeR,"Pclustering_TPM_log_edger")


make_dendrogram_pearson(dat_log_limma,"Pclustering_TPM_log_limma_new2")

make_dendrogram_spearman(dat_raw,"Sclustering_TPM")
make_dendrogram_spearman(dat_raw_log,"Sclustering_TPM_log")
make_dendrogram_spearman(dat_raw_limma,"Sclustering_TPM_limma")
make_dendrogram_spearman(dat_log_limma,"Sclustering_TPM_log_limma")
make_dendrogram_spearman(dat_raw_combat,"Sclustering_TPM_raw_combat")
make_dendrogram_spearman(dat_log_combat,"Sclustering_TPM_log_combat")
make_dendrogram_spearman(dat_raw_RUVseq,"Sclustering_counts_ruvseq")
# make_dendrogram_spearman(dat_raw_edgeR,"Sclustering_TPM_raw_edger")
# make_dendrogram_spearman(dat_log_edgeR,"Sclustering_TPM_log_edger")


# make_dendrogram_euclidean(dat_log_limma2,"mm10_log_signal_limma_cluster2")
# make_dendrogram_euclidean(dat_1,"mm10_raw_signal_rep1_cluster")
# make_dendrogram_euclidean(dat_avg,"mm10_raw_signal_rep_avg_cluster")
# make_dendrogram_pearson(dat_no0var,"mm10_raw_signal_rep_avg_novar_cluster")

############
# for comparison with Arjan's
############

make_dendrogram_euclidean(dat2_no0var_TPM, "Eclustering_no0var_TPM")
make_dendrogram_euclidean(dat2_no0var_TPM_log, "Eclustering_no0var_TPM_log")
make_dendrogram_euclidean(dat2_no0var_TPM_limma, "Eclustering_no0var_TPM_limma")
make_dendrogram_euclidean(dat2_no0var_TPM_log_limma, "Eclustering_no0var_TPM_log_limma")

make_dendrogram_pearson(dat2_no0var_TPM, "Pclustering_no0var_TPM")
make_dendrogram_pearson(dat2_no0var_TPM_log, "Pclustering_no0var_TPM_log")
make_dendrogram_pearson(dat2_no0var_TPM_limma, "Pclustering_no0var_TPM_limma")
make_dendrogram_pearson(dat2_no0var_TPM_log_limma, "Pclustering_no0var_TPM_log_limma")


#############
# PCA 
#############
make_PCA_scatter(dat_raw,"mm10_raw_signal_PCA")
make_PCA_scatter(dat_raw_log,"mm10_log_signal_PCA")
make_PCA_scatter(dat_zscore,"mm10_zscore_PCA")
make_PCA_scatter(dat_raw_quan,"mm10_raw_signal_quantiled_PCA")
make_PCA_scatter(dat_log_quan,"mm10_log_signal_quantiled_PCA")
make_PCA_scatter(dat_zscore_quan,"mm10_zscore_quantiled_PCA")
make_PCA_scatter(dat_raw_limma,"mm10_raw_signal_limma_PCA")
make_PCA_scatter(dat_log_limma,"mm10_log_signal_limma_PCA")

# for(i in 1:length(tissue_list)){
#   tissue = tissue_list[i]
#   #
#   d = dat_pca_m[dat_pca_m$tissue==tissue,]
#   p = ggplot(d, aes(x=PC1, y=PC2, col=time, shape=tissue)) + geom_point(size=2) +
#     theme_minimal() + 
#     labs(title = paste("raw expression - ", tissue, sep=""))
#   print(p)
# }

#############
# tSNE
#############
make_tsne_scatter(dat_raw,"mm10_TPM_tSNE")
make_tsne_scatter(dat_raw_log,"mm10_TPM_log_tSNE")
make_tsne_scatter(dat_raw_limma,"mm10_TPM_limma_tSNE")
make_tsne_scatter(dat_log_limma,"mm10_TPM_log_limma_tSNE")
make_tsne_scatter(dat_raw_combat,"mm10_TPM_combat_tSNE")
make_tsne_scatter(dat_log_combat,"mm10_TPM_log_combat_tSNE")
make_tsne_scatter(dat_raw_RUVseq,"mm10_read_counts_ruvseq_tSNE")
# make_tsne_scatter(dat_raw_edgeR,"mm10_raw_signal_edger_tSNE")
# make_tsne_scatter(dat_log_edgeR,"mm10_log_signal_edger_tSNE")

make_tsne_scatter(dat_log_limma,"mm10_TPM_log_limma_tSNE_new")
#############
# if do not remove variance=0
#############
dat_raw_log0 = log10(dat_raw0+0.01)
dat_log_limma0 <- removeBatchEffect(dat_raw_log0, batch)
make_dendrogram_pearson(dat_log_limma,"Pclustering_TPM_log_limma_withVariance0")
make_tsne_scatter(dat_log_limma0,"mm10_log_signal_limma_tSNE_withVariance0")
#############

#############
# spikeIn96_Jul28
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

#############
