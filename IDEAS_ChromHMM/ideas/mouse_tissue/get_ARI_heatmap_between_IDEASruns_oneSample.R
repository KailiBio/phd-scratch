
# -- Kaili
# # This script is for getting ARI heatmap for given two state files. (one sample).

args = commandArgs(trailingOnly=TRUE)
workDir = args[1]
state_file1 = args[2]
state_file2 = args[3]
sample_num = as.integer(args[4])
state_num = as.integer(args[5])
name = args[6]
out_name = args[7]

# workDir = "/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/rep2/"
# state_file1 = "ctcf_samples_lung14.5.state"
# state_file2 = "ctcf_samples_rep1-impute-rep2_oneSample.state"
# sample_num = 1
# state_num = 42
# name = "rep1 vs. rep1-impute-rep2\n(lung14.5)"
# out_name = "Compare_rep1_rep1-impute-rep2_ARI_lung_14.5"

setwd(workDir)

library(RColorBrewer)
library(ggplot2)
library(data.table)

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

#############
data1 = data.frame(fread(state_file1))
data2 = data.frame(fread(state_file2))
print("finish read file.")

dat1 = data1[,5:(4+sample_num)]
dat2 = data2[,5:(4+sample_num)]
sample = colnames(data1)[5]

# get ari matrix
ari = data.frame(matrix(0, nrow = state_num, ncol = 2))
ari[,1] = 0:(state_num-1)
for(i in 0:(state_num-1)){
 datM1 =switch_state(dat1,state=i)
 datM2 =switch_state(dat2,state=i)
 ari[i+1,2] = calculate_ari(datM1, datM2)
}
colnames(ari) = c("state","ari")
print("get ari matrix.")

write.table(ari,paste(out_name,".txt", sep=""), quote = FALSE)

# make figures
ggplot(ari, aes(x = 1, y = state, fill = ari)) + 
  geom_tile() + xlab("") + scale_fill_continuous(low="#f7f7f7",high="#ef3b2c") + 
  labs(title = name) + theme_minimal() + 
  scale_y_continuous(breaks=0:(state_num-1), labels = 0:(state_num-1)) +
  scale_x_discrete(label="")

ggsave(paste(out_name,".pdf", sep=""), width=3)
ggsave(paste(out_name,".png", sep=""), width=3)

#############