
# -- Kaili
# This script is for doing gene regression for reps in DHS-bins.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/rep2/dhs_reps_state_proportion/")

library(ggplot2)
library(data.table)
library(reshape2)

#################
# get matrix
#################
exp = read.table("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/regression/mm10_RNA_protein-coding_tpm_matrix_matched.txt", 
                 row.names = 1,header=TRUE)
#
R1_matrix=matrix(0, ncol=ncol(exp), nrow=20)
colnames(R1_matrix) = rep("a",ncol(exp))
for(i in 1:ncol(exp)){
  print(i)
  sample = colnames(exp)[i]
  sample = gsub("embryonic.facial.prominence_","embryonic-facial-prominence_",sample)
  sample = gsub("neural.tube_","neural-tube_",sample)
  
  R1=NULL
  for(j in 1:20){
    data1 = read.table(paste("mm10_gene_",sample,"_DHS_reps_state_count_window_",
                             as.character(j),".txt",sep=""), row.names=1)
    #
    r1 = NULL
    model1 = lm(log(exp[rownames(data1),11]+1e-5) ~ log(data1[,2]+1e-5)+log(data1[,3]+1e-5)+
                  log(data1[,4]+1e-5)+log(data1[,5]+1e-5)+log(data1[,6]+1e-5)+
                  log(data1[,7]+1e-5)+log(data1[,8]+1e-5)+log(data1[,9]+1e-5)+
                  log(data1[,10]+1e-5)+log(data1[,11]+1e-5)+log(data1[,12]+1e-5)+
                  log(data1[,13]+1e-5)+log(data1[,14]+1e-5)+log(data1[,15]+1e-5)+
                  log(data1[,16]+1e-5)+log(data1[,17]+1e-5)+log(data1[,18]+1e-5)+
                  log(data1[,19]+1e-5)+log(data1[,20]+1e-5)+log(data1[,21]+1e-5)+
                  log(data1[,22]+1e-5)+log(data1[,23]+1e-5)+log(data1[,24]+1e-5)+
                  log(data1[,25]+1e-5)+log(data1[,26]+1e-5)+log(data1[,27]+1e-5)+
                  log(data1[,28]+1e-5)+log(data1[,29]+1e-5)+log(data1[,30]+1e-5)+
                  log(data1[,31]+1e-5)+log(data1[,32]+1e-5)+log(data1[,33]+1e-5)+
                  log(data1[,34]+1e-5)+log(data1[,35]+1e-5)+log(data1[,36]+1e-5)+
                  log(data1[,37]+1e-5)+log(data1[,38]+1e-5)+log(data1[,39]+1e-5)+
                  log(data1[,40]+1e-5)+log(data1[,41]+1e-5)+log(data1[,42]+1e-5)+
                  log(data1[,43]+1e-5)+log(data1[,44]+1e-5)+log(data1[,45]+1e-5))
    r1 = summary(model1)$adj.r.squared
    R1[j] = round(r1,4)
  }
  R1_matrix[,i] = R1
  colnames(R1_matrix)[i] = sample
}
rownames(R1_matrix) = c(1:20)

write.table(R1_matrix, "adj_r_square_dhs_reps.txt", sep="\t", quote=FALSE,
            row.names = TRUE, col.names = TRUE)


R1_matrix = read.table("adj_r_square_dhs_reps.txt",  sep="\t", header = TRUE)
R2_matrix = read.table("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/regression/adj_r_square_dhs_bins.txt",  
                       sep="\t", header = TRUE)
R3_matrix = read.table("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/regression/adj_r_square_normal_bins.txt",  
                       sep="\t", header = TRUE)

#################
# make regression figure
#################
R1_matrix_melt = melt(R1_matrix, id = NULL)
R2_matrix_melt = melt(R2_matrix, id = NULL)
R3_matrix_melt = melt(R3_matrix, id = NULL)
matrix = data.frame(rbind(cbind(R1_matrix_melt, type="reps"),
                          cbind(R2_matrix_melt, type="rep1"),
                          cbind(R3_matrix_melt, type="200bp")))

matrix$variable = factor(rownames(R1_matrix), levels = 1:20)

ggplot(matrix, aes(x=variable, y=value, fill = type)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  theme_minimal() +
  scale_fill_manual(values = c("#B79F00","#f8766d","#00bfc4")) +
  ylab("adjusted R-squared") + xlab("TSS±2kb") +
  theme(title = element_text(face="bold", size=12),
        axis.text.x = element_text(size=12, face="bold", angle=45),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels=c("-2000,-1800","-1800,-1600","-1600,-1400","-1400,-1200","-1200,-1000",
                            "-1000,-800","-800,-600","-600,-400","-400,-200","-200,0","0,200","200,400",
                            "400,600","600,800","800,1000","1000,1200","1200,1400","1400,1600","1600,1800",
                            "1800,2000")) +
  theme(legend.position = c(0.85,0.2),legend.title = element_blank(),
        legend.text = element_text(face="bold",size=10))

ggsave("regression_adj_rsquare_reps__boxplot.pdf", width=8)
ggsave("regression_adj_rsquare_reps__boxplot.png", width=8)
#################
