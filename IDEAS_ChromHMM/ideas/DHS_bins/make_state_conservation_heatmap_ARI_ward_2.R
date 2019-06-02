
# -- Kaili
# This script is for making ARI heatmap between states.
# EXP: Rscript make_state_conservation_heatmap_ARI.R 
#       /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/ari/
#       /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_v3_100-400bp.state 0 66 dhs_bins

args = commandArgs(trailingOnly=TRUE)
workDir = args[1]
state_file = args[2]
state = as.integer(args[3])
sample_num = as.integer(args[4])
prefix = args[5]

setwd(workDir)

library(RColorBrewer)
library(pheatmap)
library(data.table)
library(parallel)
#####################
correlation<-function(x,y)
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

#switch_state <- function(vector, state){
#  out = as.integer(vector==state)
#  return(out)
#}

#####################

data = data.frame(fread(state_file))
print("finish read file.")
dat = data[,5:(4+sample_num)]
sample = colnames(dat)


#datM =apply(dat,2,switch_state, state=state)
#print("finish switch states")
# get ARI matrix
ari = matrix(0, ncol = 66, nrow = 66)
colnames(ari) = rownames(ari) = sample
for (m in 1:66){
  print(m)
  ari[m,] = unlist(mclapply(1:66, function(i){
    correlation(dat[,m], dat[,i])
  }, mc.cores = 8, mc.allow.recursive = TRUE))
}
print("finish ARI")
# make figures
pdf_name = paste("state_conservation_",prefix,"_state_",state,".pdf",sep="")
pdf(pdf_name ,width = 15, height=15)
pheatmap(ari, display_numbers = FALSE, main=paste(prefix, ": state ",as.character(state),sep=""),
         breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),  clustering_method = "ward.D2",
         col = brewer.pal(10,"RdYlBu")[10:1])
dev.off()
