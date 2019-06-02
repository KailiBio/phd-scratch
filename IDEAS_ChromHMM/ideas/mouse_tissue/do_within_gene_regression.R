
# -- Kaili
# This script is for doing gene expression regression after add CTCF.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/gene_expression_regression/")

library(reshape2)
library(ggplot2)

exp = read.table("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/regression/mm10_RNA_protein-coding_tpm_matrix_matched.txt", row.names = 1,
                 header=TRUE)

###################
# get r-square for all
###################
R_matrix=matrix(0, ncol=ncol(exp), nrow=20)
colnames(R_matrix) = rep("a",ncol(exp))
for(i in 1:ncol(exp)){
  print(i)
  sample = colnames(exp)[i]
  sample = gsub("embryonic.facial.prominence_","embryonic-facial-prominence_",sample)
  sample = gsub("neural.tube_","neural-tube_",sample)
  
  R=NULL
  for(j in 1:20){
    data1 = read.table(paste("./state_proportion/mm10_gene_",sample,"_dhs-ctcf_state_count_window_",
                             as.character(j),".txt",sep=""), row.names=1)
    #
    r = NULL
    model = lm(log(exp[rownames(data1),11]+1e-5) ~ log(data1[,2]+1e-5)+log(data1[,3]+1e-5)+
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
                  log(data1[,43]+1e-5)+log(data1[,44]+1e-5)+log(data1[,45]+1e-5)+
                 log(data1[,46]+1e-5)+log(data1[,47]+1e-5)+log(data1[,48]+1e-5))
    r = summary(model)$adj.r.squared
    R[j] = round(r,4)
  }
  R_matrix[,i] = R
  colnames(R_matrix)[i] = colnames(R_matrix)[i] = sample
}
rownames(R_matrix) = rownames(R_matrix) = c(1:20)

write.table(R_matrix, "adj_r_square_addCTCF.txt", sep="\t", quote=FALSE,
            row.names = TRUE, col.names = TRUE)


R_matrix = read.table("adj_r_square_addCTCF.txt",  sep="\t", row.names = 1, header = TRUE)
R2_matrix = read.table("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/regression/adj_r_square_dhs_bins.txt",  
                       sep="\t", row.names = 1, header = TRUE)

###################
# for state proportion
###################
R_matrix_melt = melt(R_matrix, measure.vars = 1:66)
R2_matrix_melt = melt(R2_matrix, measure.vars = 1:66)
matrix = data.frame(rbind(cbind(R_matrix_melt, type="addCTCF"),
                          cbind(R2_matrix_melt, type="noCTCT")))

matrix$variable = factor(rownames(R2_matrix), levels = 1:20)

ggplot(matrix, aes(x=variable, y=value, fill = type)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  theme_minimal() +
  scale_fill_manual(values = c("#FBB30B","#397AF2")) +
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

ggsave("adj_rsquare_addCTCF_20windows_boxplot.pdf", width=10)
ggsave("adj_rsquare_addCTCF_20windows_boxplot.png", width=10)

p=NULL
for(i in 1:20){
  p[i] = t.test(as.numeric(R_matrix[i,]),as.numeric(R2_matrix[i,]), paired = TRUE)$p.value
  print(i)
  print(p[i])
}
q = p.adjust(p, method = "bonferroni", n=20)

format(q, scientific = F)

a = apply(R1_matrix,1,median)
b = apply(R2_matrix,1,median)
a<b
###################