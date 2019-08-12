
# -- Kaili
# This script is for making CTCF signal based curve.

library(PRROC)
library(RColorBrewer)
library(data.table)
library(ROCR)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
workDir = args[1]
pos_file_name = args[2]
neg_file_name = args[3]
title = args[4]

setwd(workDir)
#setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/recall_AUC/")

##############
pos_file = data.frame(fread(pos_file_name))
neg_file = data.frame(fread(neg_file_name))

# pos_file = data.frame(fread("liver_14.5_AUC_positive.txt"))
# neg_file = data.frame(fread("liver_14.5_AUC_negative.txt"))

# get AUC
auc <- pr.curve(scores.class0 = pos_file[,2], scores.class1 = neg_file[,2])$auc.integral

# plot PR curve
df_p = data.frame(pos_file[,2],"1")
df_n = data.frame(neg_file[,2], "0")
colnames(df_p) =  colnames(df_n) = c("value", "labels")
df = data.frame(rbind(df_p, df_n))
pred <- prediction(df$value, df$labels)
perf <- performance(pred,"prec","tpr")
#
dat = data.frame(perf@x.values, perf@y.values)
colnames(dat) = c("recall", "precision")
#
name  = paste("CTCF signal based\nAUC=",round(auc,2),"\n",title,sep="")
#
ggplot(dat, aes(x=recall, y=precision)) + geom_point(col = "grey") +
         theme_minimal() + theme(title = element_text(size=12, face="bold")) + 
         labs(title = name)
ggsave(paste(title,".png",sep=""))
