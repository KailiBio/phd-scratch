
# -- Kaili
# This script is for calling DE by DEseq2.

setwd("/data/zusers/fankaili/ccre/mm10_rnaseq")
# setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/mm10_RNAseq/")

library(ggplot2)
library(DESeq2)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)
filename = args[1]
version = args[2]

dat0 = read.table(filename, header = TRUE, row.names = 1)
dat <- apply (dat0, c (1, 2), function (x) {x[x<0] <- 0; as.integer(x)})

##############
registerDoParallel(20)
foreach(i=1:20) %dopar% {
  a1 = 2*i-1
  a2 = 2*i
  for(j in (i+1):66){
    b1 = 2*j-1
    b2 = 2*j
    #
    d = dat[, c(a1,a2,b1,b2)]
    #
    coldata = data.frame(sapply(1:ncol(d), function(i) paste(strsplit(colnames(d)[i],"_")[[1]][1],
                                                             strsplit(colnames(d)[i],"_")[[1]][2],sep="_")))
    rownames(coldata) = colnames(d)
    colnames(coldata) = c("sample")
    #
    t = as.vector(unique(coldata$sample))
    title = paste("./de_", version, "/", as.character(t[2]),"_VS_", as.character(t[1]), ".txt" ,sep = "")
    #
    dds <- DESeqDataSetFromMatrix(countData = d, colData = coldata, design= ~ sample)
    # # normalization
    # dds <- estimateSizeFactors(dds)
    # sizeFactors(dds)
    # normalized_counts <- counts(dds, normalized=TRUE)
    #
    dds <- DESeq(dds)
    res <- results(dds)
    res=res[order(res$padj),]
    write.table(res, title, sep="\t", quote=F)
  }
}
