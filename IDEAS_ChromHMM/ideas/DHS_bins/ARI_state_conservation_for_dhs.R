
# -- Kaili
# This script is for calculating adjusted Rand index & Matthews correlation coefficient for state conversation between samples.
# ARI & MCC

setwd("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_v3_100-400bp_result/")

library(RColorBrewer)
library(pheatmap)

args = commandArgs(trailingOnly=TRUE)
state_num=args[1]


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

#####################
# get all the data
dat = NULL
for(i in c(1:19,"X","Y")){
  file = paste("DHS_v3_100-400bp.chr",i,".state", sep="")
  data = read.table(file, header = TRUE, row.names = 1, sep=" ",comment.char = "!")
  if(i==1){
    dat = data
  }else{
    dat = rbind(dat,data)
  }
}

######------------------
# using whole genome as background
sample = c("intestine_16.5", "stomach_16.5", "limb_11.5", "kidney_15.5", "embryonic.facial.prominence_15.5",
           "midbrain_14.5", "embryonic.facial.prominence_11.5", "hindbrain_11.5", "neural.tube_14.5",
           "liver_15.5", "midbrain_16.5", "neural.tube_11.5", "embryonic.facial.prominence_14.5",
           "kidney_16.5", "neural.tube_13.5", "midbrain_13.5", "midbrain_11.5", "hindbrain_13.5",
           "hindbrain_16.5", "intestine_14.5", "hindbrain_14.5", "heart_11.5", "stomach_0", "liver_11.5",
           "heart_15.5", "hindbrain_15.5", "heart_16.5", "liver_16.5", "midbrain_0", "kidney_0",
           "liver_13.5", "heart_14.5", "lung_15.5", "lung_16.5", "liver_0", "stomach_15.5", "intestine_0",
           "stomach_14.5", "kidney_14.5", "lung_14.5", "lung_0", "heart_0", "neural.tube_15.5", "hindbrain_0",
           "limb_15.5", "embryonic.facial.prominence_12.5", "neural.tube_12.5", "intestine_15.5", "limb_13.5",
           "embryonic.facial.prominence_13.5", "limb_12.5", "hindbrain_12.5", "liver_12.5", "midbrain_12.5",
           "midbrain_15.5", "heart_12.5", "heart_13.5", "liver_14.5", "limb_14.5", "forebrain_14.5",
           "forebrain_12.5", "forebrain_13.5", "forebrain_16.5", "forebrain_15.5", "forebrain_0", "forebrain_11.5")

pdf(paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/ari/ARI_state_conservation_dhs_",state_num,".pdf",sep=""))
for(i in state_num:state_num){
  # get ARI and MCC index matrix
  ari = matrix(0, ncol = 66, nrow = 66)
  colnames(ari) = rownames(ari) = sample[]
  for (m in 1:66){
    print(i)
    sample1 = sample[m]
    x=as.integer(dat[,sample1]==i)
    # get other samples
    for(n in 1:66){
      sample2 = sample[n]
      y=as.integer(dat[,sample2]==i)
      # calculate correlation
      ari[m,n] = correlation(x,y)
    }
  }
  # make figures
  pheatmap(ari, display_numbers = TRUE,
           main=paste("ARI: state ",as.character(i),sep=""),
           breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
           col = brewer.pal(10,"RdYlBu")[10:1])
}
dev.off()

#####################
