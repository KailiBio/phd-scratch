
# -- Kaili
# This script is for making PR curve (using CTCF peak as gold standard).

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

workDir="/data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/9impute11/"
pos_file_name="9impute11_liver_14.5_peak_PR_pos.txt"
neg_file_name="9impute11_liver_14.5_peak_PR_neg.txt"
title="liver_14.5_9impute11_PRcurve_peak"


setwd(workDir)


##############
pos_file = data.frame(fread(pos_file_name, sep="\t"))
neg_file = data.frame(fread(neg_file_name, sep="\t"))

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
name  = paste("CTCF peak based, states ranked\nAUC=",round(auc,2),"\n",title,sep="")
#
ggplot(dat, aes(x=recall, y=precision)) + geom_line(size=1.2) +
  theme_minimal() + theme(title = element_text(size=12, face="bold")) +
  labs(title = name) + coord_cartesian(ylim=c(0,1), xlim=c(0,1))
ggsave(paste(title,".png",sep=""))
ggsave(paste(title,".pdf",sep=""))
