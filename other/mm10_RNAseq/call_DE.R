
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
dat <- apply (dat0, c (1, 2), function (x) {x[x<0] <- 0; round(x)})

# data0 = read.table("mm10_M4_counts_matrix.txt", row.names = 1, header = TRUE)
# rm_list=c("forebrain_10.5_1", "forebrain_10.5_2", "midbrain_10.5_1", "midbrain_10.5_2", "hindbrain_10.5_1",
#           "hindbrain_10.5_2", "heart_10.5_1", "heart_10.5_2", "facial_10.5_1", "facial_10.5_2",
#           "limb_10.5_1", "limb_10.5_2", "neural.tube_0_1", "neural.tube_0_2")
# dat0 = data0[,!colnames(data0)%in%rm_list]
# dat <- apply (dat0, c (1, 2), function (x) {x[x<0] <- 0; round(x)})
# write.table(dat, "mm10_M4_read_66.txt", sep="\t", quote=F)

##############
registerDoParallel(20)
foreach(i=1:73) %dopar% {
  a1 = 2*i-1
  a2 = 2*i
  for(j in (i+1):73){
    b1 = 2*j-1
    b2 = 2*j
    #
    d = dat[, c(a1,a2,b1,b2)]
    #
    coldata = data.frame(sapply(1:ncol(d), function(i) paste(strsplit(colnames(d)[i],"_")[[1]][1],
                                                             strsplit(colnames(d)[i],"_")[[1]][2],sep="_")),
                         sapply(1:ncol(d), function(i) strsplit(colnames(d)[i],"_")[[1]][2]))
    rownames(coldata) = colnames(d)
    colnames(coldata) = c("sample", "batch")
    #
    t = as.vector(unique(coldata$sample))
    title = paste("./de_", version, "/", as.character(t[2]),"_VS_", as.character(t[1]), ".txt" ,sep = "")
    #
    dds <- DESeqDataSetFromMatrix(countData = d, colData = coldata, design= ~ sample)
    # # normalization
    # dds <- estimateSizeFactors(dds)
    # sizeFactors(dds)
    # normalized_counts <- counts(dds, normalized=TRUE)
    dds <- DESeq(dds, betaPrior=TRUE)
    #dds_log <- rlog(dds)
    res <- results(dds)
    res=res[order(res$padj),]
    write.table(res, title, sep="\t", quote=F)
  }
}

####################
# 
# library("tximport")
# library("readr")
# library("tximportData")
# 
# 
# # facial_12.5 vs facial_11.5
# 
# samples <- data.frame(c("facial_11.5", "facial_11.5", "facial_12.5", "facial_12.5"))
# colnames(samples) = c("sample")
# rownames(samples) = c("facial_11.5_1", "facial_11.5_2", "facial_12.5_1", "facial_12.5_2")
# 
# 
# files <- c("/data/zusers/fankaili/ccre/mm10_rnaseq/facial_11.5_1.tsv",
#            "/data/zusers/fankaili/ccre/mm10_rnaseq/facial_11.5_2.tsv",
#            "/data/zusers/fankaili/ccre/mm10_rnaseq/facial_12.5_1.tsv",
#            "/data/zusers/fankaili/ccre/mm10_rnaseq/facial_12.5_2.tsv")
# names(files) <- rownames(samples)
# 
# txi <- tximport(files, type="rsem",txIn = FALSE, txOut = FALSE)
# sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 2)), batch=c(1,1,2,2))
# rownames(sampleTable) <- colnames(txi$counts)
# txi$length <- apply (txi$length, c (1, 2), function (x) {x[x<0] <- 1})
# 
# ddsTxi <- DESeqDataSetFromTximport(txi, sampleTable, ~condition+batch)
# dds <- DESeq(ddsTxi)
# res <- results(dds)
# res=res[order(res$padj),]
# 
# write.table(res, "/data/zusers/fankaili/ccre/mm10_rnaseq/facial_11.5_VS_facial_12.5.ss.txt", 
#             sep="\t", quote=F)
# 
# 
# library(apeglm)
# resultsNames(dds)
# resLFC <- lfcShrink(dds, coef="condition_B_vs_A", type="apeglm")
# resLFC=resLFC[order(resLFC$padj),]
# write.table(resLFC, "/data/zusers/fankaili/ccre/mm10_rnaseq/facial_11.5_VS_facial_12.5.ss2.txt", 
#             sep="\t", quote=F)
# 
# 


