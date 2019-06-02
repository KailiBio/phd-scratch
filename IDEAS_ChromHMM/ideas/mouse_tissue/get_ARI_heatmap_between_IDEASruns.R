
# -- Kaili
# This script is for getting ARI heatmap for given two state files.
# EXP: Rscript get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/dhs_ctcf/
#     ctcf_8sample_impute_11sample.state ctcf_8sample_impute_11sample_2.state 11 47
#     "8 samples impute 11 samples (reproducibility)"
#     "ARI_8sample_imputation_reproducibility"

args = commandArgs(trailingOnly=TRUE)
workDir = args[1]
state_file1 = args[2]
state_file2 = args[3]
sample_num = as.integer(args[4])
state_num = as.integer(args[5])
name = args[6]
out_name = args[7]

# workDir = ""
# sample_num = 11
# state_num = 47
# state_file1 = "ctcf_8sample_impute_11sample.state"
# state_file2 = "ctcf_8sample_impute_11sample_2.state"
# name = ""
# pdf_name = ""

setwd(workDir)

library(RColorBrewer)
library(pheatmap)
library(data.table)
library(parallel)

#####################
calculate_ari<-function(x,y)
{
  t=as.matrix(table(x,y));
  TP=sum(choose(t,2));
  FN=sum(choose(apply(t,1,sum),2))-TP;
  FP=sum(choose(apply(t,2,sum),2))-TP;
  TN=sum(choose(sum(t),2))-TP-FN-FP;
  a=TP+1;
  b=TN+1;
  c=FN+1;
  d=FP+1;

  ari=(a-(a+c)*(a+d)/(a+b+c+d))/((a+c+a+d)/2-(a+c)*(a+d)/(a+b+c+d));
  return(ari);
}

switch_state <- function(vector, state){
  out = as.integer(vector==state)
  return(out)
}

#####################

# get data
data1 = data.frame(fread(state_file1))
data2 = data.frame(fread(state_file2))
print("finish read file.")

dat1 = data1[,5:(4+sample_num)]
dat2 = data2[,5:(4+sample_num)]
sample = colnames(dat1)

# get ari matrix
ari = matrix(0, nrow = state_num, ncol = sample_num)
for(i in 0:(state_num-1)){
  datM1 =apply(dat1,2,switch_state, state=i)
  datM2 =apply(dat2,2,switch_state, state=i)
  ari[i+1,] = unlist(mclapply(1:sample_num, function(j){
    calculate_ari(datM1[,j], datM2[,j])
    }, mc.cores = 8, mc.allow.recursive = TRUE))
}
colnames(ari) = sample
rownames(ari) = 0:(state_num-1)
print("get ari matrix.")

#for(i in 0:(state_num-1)){
#  datM1 =switch_state(dat1,state=i)
#  datM2 =switch_state(dat2,state=i)
#  ari[i+1,1] = calculate_ari(datM1, datM2)
#}

write.table(ari,paste(out_name,".txt", sep=""), quote = FALSE)

# make figures
pdf(paste(out_name,".pdf", sep="") ,width = 8, height=10)
pheatmap(ari, display_numbers = TRUE, main=name, breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
         col = brewer.pal(10,"RdYlBu")[10:1], fontsize = 12)
dev.off()
