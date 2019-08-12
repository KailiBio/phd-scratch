
# -- Kaili
# This script is for calculating PRAU for given file.

library(PRROC)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
workDir = args[1]
pos_file_name = args[2]
neg_file_name = args[3]

setwd(workDir)


##############
pos_file = data.frame(fread(pos_file_name, sep="\t"))
neg_file = data.frame(fread(neg_file_name, sep="\t"))

# get AUC
auc <- pr.curve(scores.class0 = pos_file[,2], scores.class1 = neg_file[,2])$auc.integral

print(auc)