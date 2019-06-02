
# -- Kaili
# This script is for calculating correlation for two given signal file.

args = commandArgs(trailingOnly=TRUE)
file1 = args[1]
file2 = args[2]
type = args[3]

# file1="/data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/heart_13.5_ATAC_normal.txt"
# file2="/data/zusers/fankaili/ideas/signal/rep2_signal_normal_bins/heart_13.5_ATAC_normal.txt"
# type="spearman"

library(data.table)
dat1 <- fread(file1)
dat2 <- fread(file2)

if(type=="spearman"){
  r = cor(dat1[,1], dat2[,1], method="spearman")[1]
}else if(type=="pearson"){
  r = cor(dat1[,1], dat2[,1], method="pearson")[1]
}

print(r)
