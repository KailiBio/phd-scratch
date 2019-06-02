
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
make_dendrogram <- function(data, title){
  dat = data.frame(t(data))
  dat$sample = colnames(data)
  dat$tissue = sapply(dat$sample, function(x) strsplit(as.character(x),"_")[[1]][1])
  dat$time = sapply(dat$sample, function(x) strsplit(as.character(x),"_")[[1]][2])
  
  dd <- dist(dat[,c(1:47862)], diag=TRUE)
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
      ligne=match(attributes(n)$label,dat[,47863])
      tissue=dat[ligne,47864];
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
      time = dat[ligne,47865];
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
  
  pdf(paste(title,".pdf",sep=""), height=15)
  par(mai=c(0.5,0.5,0.5,1.5))
  plot(hc_col , main=title, horiz=TRUE, hang = -0.1)
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
  ggsave(paste(title,".pdf",sep=""))
  ggsave(paste(title,".png",sep=""))
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

# raw data matrix
dat0 = read.table("mm10_exp_matrix.txt", row.names = 1, header = 1)
dat_raw = dat0[ , order(colnames(dat0))]

# z-score matrix
dat1 = read.table("mm10_exp_zscore_matrix.txt", row.names = 1, header = 1)
dat_zscore = dat1[, order(colnames(dat1))]

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
write.table(dat_log_limma, "mm10_exp_log-limma_matrix.txt", sep="\t", quote = FALSE, 
            row.names = TRUE, col.names = TRUE)

batch2 = sapply(colnames(dat_raw), function(x) strsplit(as.character(x),"_")[[1]][1])
dat_log_limma2 <- removeBatchEffect(dat_raw_log, batch2)
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
# clutering for raw data
#############
make_dendrogram(dat_raw,"mm10_raw_signal_cluster")
make_dendrogram(dat_raw_log,"mm10_raw_signal_log_cluster")
make_dendrogram(dat_zscore,"mm10_zscore_cluster")
make_dendrogram(dat_raw_quan,"mm10_raw_signal_quantiled_cluster")
make_dendrogram(dat_log_quan,"mm10_log_signal_quantiled_cluster")
make_dendrogram(dat_zscore_quan,"mm10_zscore_quantiled_cluster")
make_dendrogram(dat_raw_limma,"mm10_raw_signal_limma_cluster")
make_dendrogram(dat_log_limma,"mm10_log_signal_limma_cluster")

make_dendrogram(dat_log_limma2,"mm10_log_signal_limma_cluster2")
#############
# PCA for raw exp
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
