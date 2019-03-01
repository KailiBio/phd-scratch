
# -- Kaili
# This script is for calculating adjusted Rand index & Matthews correlation coefficient for state conversation between samples.
# ARI & MCC
# ARI and MCC are the same. change to only ARI here.

setwd("/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/")
#setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/CTCF_state_file/")
library(RColorBrewer)
library(pheatmap)


# args = commandArgs(trailingOnly=TRUE)
# state_num=args[1]

#####################
ari<-function(x,y)
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
  file = paste("ctcf_samples.chr",i,".state", sep="")
  data = read.table(file, header = TRUE, row.names = 1, sep=" ",comment.char = "!")
  if(i==1){
    dat = data
  }else{
    dat = rbind(dat,data)
  }
}

######------------------
# using whole genome as background
sample = c("lung_14.5", "liver_14.5", "stomach_0", "midbrain_0", "kidney_0", "liver_0", "intestine_0",
           "lung_0", "heart_0", "hindbrain_0", "forebrain_0")
pdf(paste("/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ari_mcc/ARI_state_conservation_11sample_CTCF.pdf",sep=""))
for(i in 0:42){
  # get ARI and MCC index matrix
  ari_matrix = matrix(0, ncol = 11, nrow = 11)
  colnames(ari_matrix) = rownames(ari_matrix) = sample
  for (m in 1:11){
    sample1 = sample[m]
    x=as.integer(dat[,sample1]==i)
    # get other samples
    for(n in 1:11){
      sample2 = sample[n]
      y=as.integer(dat[,sample2]==i)
      # calculate correlation
      ari_matrix[m,n] = ari(x,y)
    }
  }
  # make figures
  pheatmap(ari_matrix, display_numbers = TRUE,
           main=paste("ARI: state ",as.character(i),sep=""),
           breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
           col = brewer.pal(10,"RdYlBu")[10:1])
}
dev.off()

######------------------
# # use union bins as background
# #
# sample = c("lung_14.5", "liver_14.5", "stomach_0", "midbrain_0", "kidney_0", "liver_0", "intestine_0",
#            "lung_0", "heart_0", "hindbrain_0", "forebrain_0")
# for(i in start:end){
#   # get ARI and MCC index matrix
#   ari = matrix(0, ncol = 11, nrow = 11)
#   mcc = matrix(0, ncol = 11, nrow = 11)
#   colnames(ari) = rownames(ari) = sample
#   colnames(mcc) = rownames(mcc) = sample
#   for (m in 1:11){
#     sample1 = sample[m]
#     xx = dat[,sample1]
#     # get other samples
#     for(n in 1:11){
#       sample2 = sample[n]
#       yy = dat[,sample2]
#       matrix0 = cbind(xx,yy)
#       matrix = matrix0[matrix0[,1]==i | matrix0[,2]==i,]
#       x0 = matrix[,1]
#       y0 = matrix[,2]
#       #
#       if(i==0){
#         x1 = replace(x0,x0!=0,100)
#         x2 = replace(x1,x1==0,1)
#         x = replace(x2,x2==100,0)
#       }else if(i==1){
#         x = replace(x0,x0!=1,0)
#       }else{
#         x1 = replace(x0,x0==i,1)
#         x = replace(x1,x1!=1,0)
#       }
#       #
#       if(i==0){
#         y1 = replace(y0,y0!=0,100)
#         y2 = replace(y1,y1==0,1)
#         y = replace(y2,y2==100,0)
#       }else if(i==1){
#         y = replace(y0,y0!=1,0)
#       }else{
#         y1 = replace(y0,y0==i,1)
#         y = replace(y1,y1!=1,0)
#       }
#       # calculate correlation
#       r = correlation(x,y)
#       ari[m,n] = r[1]
#       mcc[m,n] = r[2]
#     }
#   }
#   # make figures
#   pdf(paste("/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ari_mcc/ARI_state_conservation_",i,".pdf",sep=""))
#   pheatmap(ari, display_numbers = TRUE, 
#            main=paste("ARI: state ",as.character(i),sep=""), 
#            breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
#            col = brewer.pal(10,"RdYlBu")[10:1])
#   dev.off()
#   #
#   pdf(paste("/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ari_mcc/MCC_state_conservation_",i,".pdf",sep=""))
#   pheatmap(mcc, display_numbers = TRUE, 
#            main=paste("MCC: state ",as.character(i),sep=""), 
#            breaks = c(-1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
#            col = c("#000000",brewer.pal(10,"RdYlBu")[10:1]))
#   dev.off()
# }


#####################