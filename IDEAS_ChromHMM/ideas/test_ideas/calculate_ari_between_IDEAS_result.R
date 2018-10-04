
# -- Kaili
# This script is for calculating adjusted Rand index of two IDEAS results.
# EXP: Rscript calculate_ari_between_IDEAS_result.R "/data/zusers/fankaili/ideas/yu_input/IDEAS_yu_input_2_result/" \
#      "IDEAS_yu_input_2.chr" "/data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/" \
#      "run_IDEAS_8hm_atac_dname_pvalue.chr" "yuInput_pvalue"

setwd("/data/zusers/fankaili/ideas/ari/")

args = commandArgs(trailingOnly=TRUE)
dir1 = args[1]
prefix1 = args[2]
dir2 = args[3]
prefix2 = args[4]
out = args[5]

# setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/")
# dir1 = "/data/zusers/fankaili/ideas/yu_input/IDEAS_yu_input_2_result/"
# prefix1 = "IDEAS_yu_input_2.chr"
# dir2 = "/data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/"
# prefix2 = "run_IDEAS_8hm_atac_dname_pvalue.chr"
# out = "yuInput_pvalue"
# 
# file1 = "IDEAS_yu_input_2.chr10.state"
# file2 = "run_IDEAS_8hm_atac_dname_pvalue.chr10.state"
###############
# function
###############
calculate_ari<-function(line)
{
  
  t=as.matrix(table(line[1:66],line[67:132]));
  a=sum(choose(t,2));#TP
  c=sum(choose(apply(t,1,sum),2)) - a; #FN
  d=sum(choose(apply(t,2,sum),2)) -  a; #FP
  b=sum(choose(sum(t),2))-a-c-d;#TN
  a=a+1;
  b=b+1;
  c=c+1;
  d=d+1;
  
  ri=(a+b)/(a+b+c+d);
  ari=(a-(a+c)*(a+d)/(a+b+c+d))/((a+c+a+d)/2-(a+c)*(a+d)/(a+b+c+d));
  jac=a/(a+c+d);
  fm=sqrt(a/(a+d)*a/(a+c));
  return(c(ri,ari,jac,fm));#,a-1,b-1,c-1,d-1));
}

myrand<-function(x,y)
{
  t=as.matrix(table(x,y));
  a=sum(choose(t,2));#TP
  c=sum(choose(apply(t,1,sum),2)) - a; #FN
  d=sum(choose(apply(t,2,sum),2)) -  a; #FP
  b=sum(choose(sum(t),2))-a-c-d;#TN
  a=a+1;
  b=b+1;
  c=c+1;
  d=d+1;
  
  ri=(a+b)/(a+b+c+d);
  ari=(a-(a+c)*(a+d)/(a+b+c+d))/((a+c+a+d)/2-(a+c)*(a+d)/(a+b+c+d));
  jac=a/(a+c+d);
  fm=sqrt(a/(a+d)*a/(a+c));
  return(c(ri,ari,jac,fm));#,a-1,b-1,c-1,d-1));
}

###############
# data processing
###############
rand_bin = matrix(0, 21,2)
rand_pos = matrix(0, 21,2)

n=0
#for(i in c(1:19,"X","Y")){
for(i in c(1:2)){
  n=n+1
  
  file1 = paste(dir1,prefix1,i,".state", sep="")
  file2 = paste(dir2,prefix2,i,".state", sep="")
  
  data1 = read.table(file1, header = TRUE, row.names = 1, sep=" ",comment.char = "!")
  data2 = read.table(file2, header = TRUE, row.names = 1, sep=" ",comment.char = "!")
  
  dat1 = data1[,order(colnames(data1)[-c(1:3,70)])+3]
  dat2 = data2[,order(colnames(data2)[-c(1:3,70)])+3]
  
  d1 = dat1[order(rownames(dat1)),]
  d2 = dat2[order(rownames(dat2)),]
  
  matrix = cbind(d1,d2)
  
  index = apply(matrix, 1, calculate_ari)
  rand_bin[n,1] = mean(index[1,])
  rand_bin[n,2] = mean(index[2,])
  rand_pos[n,1:2] = myrand(data1$PosClass, data2[rownames(data1),]$PosClass)[1:2]
}

colnames(rand_bin) = colnames(rand_pos) = c("ri", "ari")
#rownames(rand_bin) = rownames(rand_pos) = paste("chr",c(1:19,"X","Y"), sep="")
rownames(rand_bin) = rownames(rand_pos) = paste("chr",c(1:2), sep="")

write.table(rand_bin, paste(out,"_rand_bin.txt",sep=""),sep="\t", row.names = TRUE, 
            col.names = TRUE, quote = FALSE)
write.table(rand_pos, paste(out,"_rand_posClass.txt",sep=""),sep="\t", row.names = TRUE, 
            col.names = TRUE, quote = FALSE)
