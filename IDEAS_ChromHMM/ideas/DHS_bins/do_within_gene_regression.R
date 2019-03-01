
# -- Kaili
# This script is for making aggregation within cell type.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/regression/")
library(reshape2)

exp = read.table("mm10_RNA_protein-coding_tpm_matrix_matched.txt", row.names = 1,
                 header=TRUE)

###################
# try with H3K4me3
###################
ocr = read.table("try_regression_H3K4me3_heart12.5_ocr.txt")
rownames(ocr) = ocr[,1]

exp1 = read.table("gene_exp_tmp_ocr_heart12.5.txt")
rownames(exp1) = exp1[,1]

names = intersect(rownames(exp), rownames(ocr))
y = exp1[names,]
x = ocr[names,]

cor(y[,2], x[,2], method = "spearman")
model1 = lm(log10(y[,2]+1e-5)~log(x[,2]+1e-5))
summary(model1)$adj.r.squared


###################
# get r-square for all
###################
R1_matrix=matrix(0, ncol=ncol(exp), nrow=20)
R2_matrix=matrix(0, ncol=ncol(exp), nrow=20)
colnames(R1_matrix) = colnames(R2_matrix) = rep("a",ncol(exp))
for(i in 1:ncol(exp)){
  print(i)
  sample = colnames(exp)[i]
  sample = gsub("embryonic.facial.prominence_","embryonic-facial-prominence_",sample)
  sample = gsub("neural.tube_","neural-tube_",sample)
  
  R1=NULL
  R2=NULL
  for(j in 1:20){
    data1 = read.table(paste("./normal_window/mm10_gene_",sample,"_normal_state_count_window_",
                             as.character(j),".txt",sep=""), row.names=1)
    data2 = read.table(paste("./dhs_window/mm10_gene_",sample,"_dhs_state_count_window_",
                             as.character(j),".txt",sep=""), row.names=1)
    #
    r1 = NULL
    r2 = NULL
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
    model2 = lm(log(exp[rownames(data2),11]+1e-5) ~ log(data2[,2]+1e-5)+log(data2[,3]+1e-5)+
                  log(data2[,4]+1e-5)+log(data2[,5]+1e-5)+log(data2[,6]+1e-5)+
                  log(data2[,7]+1e-5)+log(data2[,8]+1e-5)+log(data2[,9]+1e-5)+
                  log(data2[,10]+1e-5)+log(data2[,11]+1e-5)+log(data2[,12]+1e-5)+
                  log(data2[,13]+1e-5)+log(data2[,14]+1e-5)+log(data2[,15]+1e-5)+
                  log(data2[,16]+1e-5)+log(data2[,17]+1e-5)+log(data2[,18]+1e-5)+
                  log(data2[,19]+1e-5)+log(data2[,20]+1e-5)+log(data2[,21]+1e-5)+
                  log(data2[,22]+1e-5)+log(data2[,23]+1e-5)+log(data2[,24]+1e-5)+
                  log(data2[,25]+1e-5)+log(data2[,26]+1e-5)+log(data2[,27]+1e-5)+
                  log(data2[,28]+1e-5)+log(data2[,29]+1e-5)+log(data2[,30]+1e-5)+
                  log(data2[,31]+1e-5)+log(data2[,32]+1e-5)+log(data2[,33]+1e-5)+
                  log(data2[,34]+1e-5)+log(data2[,35]+1e-5)+log(data2[,36]+1e-5)+
                  log(data2[,37]+1e-5)+log(data2[,38]+1e-5)+log(data2[,39]+1e-5)+
                  log(data2[,40]+1e-5)+log(data2[,41]+1e-5)+log(data2[,42]+1e-5)+
                  log(data2[,43]+1e-5)+log(data2[,44]+1e-5))
    r1 = summary(model1)$adj.r.squared
    r2 = summary(model2)$adj.r.squared
    R1[j] = round(r1,4)
    R2[j] = round(r2,4)
  }
  R1_matrix[,i] = R1
  R2_matrix[,i] = R2
  colnames(R1_matrix)[i] = colnames(R2_matrix)[i] = sample
}
rownames(R1_matrix) = rownames(R2_matrix) = c(1:20)

write.table(R1_matrix, "adj_r_square_normal_bins.txt", sep="\t", quote=FALSE,
            row.names = TRUE, col.names = TRUE)
write.table(R2_matrix, "adj_r_square_dhs_bins.txt", sep="\t", quote=FALSE,
            row.names = TRUE, col.names = TRUE)


R1_matrix = read.table("adj_r_square_normal_bins.txt",  sep="\t", row.names = 1, header = TRUE)
R2_matrix = read.table("adj_r_square_dhs_bins.txt",  sep="\t", row.names = 1, header = TRUE)
###################
# for state proportion
###################
pdf("within_gene_expression_regression.pdf")
for(i in 1:66){
  sample = colnames(exp)[i]
  sample = gsub("embryonic.facial.prominence_","embryonic-facial-prominence_",sample)
  sample = gsub("neural.tube_","neural-tube_",sample)
  #
  matrix = data.frame(rbind(cbind(1:20,R1_matrix[,i],"normal_bins"),
                            cbind(1:20,R2_matrix[,i],"dhs_bins")))
  colnames(matrix) = c("window","r", "type")
  matrix$window = factor(matrix$window, levels = 1:20)
  
  p = ggplot(matrix, aes(x=window, y=r, group=type)) + 
    geom_line(aes(col=type),size=1.5) +
    theme_minimal() + geom_point(size=1) +
    xlab("") + ylab("adjusted R-squared") +
    theme(title=element_text(face="bold", size=12),
          axis.text.x = element_blank()) + 
    labs(title=sample) +
    scale_x_discrete(name="TSS±2kb", breaks=c(1,5,10,15,20)) +
    theme(legend.position = c(0.15,0.85), 
          legend.text = element_text(face="bold", size=12))
  print(p)
}
dev.off()

# pdf("within_gene_expression_regression.pdf")
# for(i in 1:ncol(exp)){
#   sample = colnames(exp)[i]
#   sample = gsub("embryonic.facial.prominence_","embryonic-facial-prominence_",sample)
#   sample = gsub("neural.tube_","neural-tube_",sample)
#   
#   R1=NULL
#   R2=NULL
#   for(j in 1:20){
#     data1 = read.table(paste("./normal_window/mm10_gene_",sample,"_normal_state_count_window_",
#                              as.character(j),".txt",sep=""), row.names=1)
#     data2 = read.table(paste("./dhs_window/mm10_gene_",sample,"_dhs_state_count_window_",
#                              as.character(j),".txt",sep=""), row.names=1)
#     #
#     r1 = NULL
#     r2 = NULL
#     model1 = lm(log(exp[rownames(data1),11]+1e-5) ~ log(data1[,2]+1e-5)+log(data1[,3]+1e-5)+
#                   log(data1[,4]+1e-5)+log(data1[,5]+1e-5)+log(data1[,6]+1e-5)+
#                   log(data1[,7]+1e-5)+log(data1[,8]+1e-5)+log(data1[,9]+1e-5)+
#                   log(data1[,10]+1e-5)+log(data1[,11]+1e-5)+log(data1[,12]+1e-5)+
#                   log(data1[,13]+1e-5)+log(data1[,14]+1e-5)+log(data1[,15]+1e-5)+
#                   log(data1[,16]+1e-5)+log(data1[,17]+1e-5)+log(data1[,18]+1e-5)+
#                   log(data1[,19]+1e-5)+log(data1[,20]+1e-5)+log(data1[,21]+1e-5)+
#                   log(data1[,22]+1e-5)+log(data1[,23]+1e-5)+log(data1[,24]+1e-5)+
#                   log(data1[,25]+1e-5)+log(data1[,26]+1e-5)+log(data1[,27]+1e-5)+
#                   log(data1[,28]+1e-5)+log(data1[,29]+1e-5)+log(data1[,30]+1e-5)+
#                   log(data1[,31]+1e-5)+log(data1[,32]+1e-5)+log(data1[,33]+1e-5)+
#                   log(data1[,34]+1e-5)+log(data1[,35]+1e-5)+log(data1[,36]+1e-5)+
#                   log(data1[,37]+1e-5)+log(data1[,38]+1e-5)+log(data1[,39]+1e-5)+
#                   log(data1[,40]+1e-5)+log(data1[,41]+1e-5)+log(data1[,42]+1e-5)+
#                   log(data1[,43]+1e-5))
#     model2 = lm(log(exp[rownames(data2),11]+1e-5) ~ log(data2[,2]+1e-5)+log(data2[,3]+1e-5)+
#                   log(data2[,4]+1e-5)+log(data2[,5]+1e-5)+log(data2[,6]+1e-5)+
#                   log(data2[,7]+1e-5)+log(data2[,8]+1e-5)+log(data2[,9]+1e-5)+
#                   log(data2[,10]+1e-5)+log(data2[,11]+1e-5)+log(data2[,12]+1e-5)+
#                   log(data2[,13]+1e-5)+log(data2[,14]+1e-5)+log(data2[,15]+1e-5)+
#                   log(data2[,16]+1e-5)+log(data2[,17]+1e-5)+log(data2[,18]+1e-5)+
#                   log(data2[,19]+1e-5)+log(data2[,20]+1e-5)+log(data2[,21]+1e-5)+
#                   log(data2[,22]+1e-5)+log(data2[,23]+1e-5)+log(data2[,24]+1e-5)+
#                   log(data2[,25]+1e-5)+log(data2[,26]+1e-5)+log(data2[,27]+1e-5)+
#                   log(data2[,28]+1e-5)+log(data2[,29]+1e-5)+log(data2[,30]+1e-5)+
#                   log(data2[,31]+1e-5)+log(data2[,32]+1e-5)+log(data2[,33]+1e-5)+
#                   log(data2[,34]+1e-5)+log(data2[,35]+1e-5)+log(data2[,36]+1e-5)+
#                   log(data2[,37]+1e-5)+log(data2[,38]+1e-5)+log(data2[,39]+1e-5)+
#                   log(data2[,40]+1e-5)+log(data2[,41]+1e-5)+log(data2[,42]+1e-5))
#     r1 = summary(model1)$adj.r.squared
#     r2 = summary(model2)$adj.r.squared
#     R1[j] = round(r1,4)
#     R2[j] = round(r2,4)
#   }
#   matrix = data.frame(rbind(cbind(1:20,R1,"normal_bins"),
#                             cbind(1:20,R2,"dhs_bins")))
#   colnames(matrix) = c("window","r", "type")
#   matrix$window = factor(matrix$window, levels = 1:20)
#   
#   p = ggplot(matrix, aes(x=window, y=r, group=type)) + 
#     geom_line(aes(col=type),size=1.5) +
#     theme_minimal() + geom_point(size=1) +
#     xlab("") + ylab("adjusted R-squared") +
#     theme(title=element_text(face="bold", size=12),
#           axis.text.x = element_blank()) + 
#     labs(title=sample) +
#     scale_x_discrete(name="TSS±2kb", breaks=c(1,5,10,15,20)) +
#     theme(legend.position = c(0.15,0.85), 
#           legend.text = element_text(face="bold", size=12))
#   print(p)
# }
# dev.off()

###################
# for each window
###################

R1_matrix_melt = melt(R1_matrix, id="group")
R2_matrix_melt = melt(R2_matrix, id="group")
matrix = data.frame(rbind(cbind(R1_matrix_melt, type="200bp_bins"),
                    cbind(R2_matrix_melt, type="dhs_bins")))

matrix$Var1 = factor(matrix$Var1, levels = 1:20)

ggplot(matrix, aes(x=Var1, y=value, fill = type)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  theme_minimal() +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
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

ggsave("adj_rsquare_66samples_20windows_boxplot.pdf")


p=NULL
for(i in 1:20){
  p[i] = t.test(as.numeric(R1_matrix[i,]),as.numeric(R2_matrix[i,]), paired = TRUE)$p.value
  print(i)
  print(p[i])
}
q = p.adjust(p, method = "bonferroni", n=20)

format(q, scientific = F)

a = apply(R1_matrix,1,median)
b = apply(R2_matrix,1,median)
a<b

###################
# for state 10
###################

matrix2 = matrix[matrix$Var1==10,]
p = formatC(t.test(matrix2[matrix2$type=="200bp_bins",]$value,matrix2[matrix2$type=="dhs_bins",]$value
                   ,paired = TRUE)$p.value, format = "e", digits = 2)



ggplot(matrix2, aes(x=Var1, y=value, fill = type)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  theme_minimal() +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  xlab("TSS upstream 200bp") + ylab("adjusted R-squared") +
  theme(title=element_text(face="bold", size=12),
        axis.text.x = element_blank()) + 
  labs(subtitle=paste("paired t-test p-value=",p,sep=""))

ggsave("state10_TSSup200bp_adjRsqaure_boxplot.pdf")
###################
# TSS±200bp, 400bp, 500bp, 2kb
###################  

get_within_celltype_rsquare <- function(length){
  Rmatrix=matrix(0,ncol=66, nrow=2)
  colnames(Rmatrix) = rep("a",66)
  rownames(Rmatrix) = c("normal", "dhs")
  #
  for(i in 1:66){
    print(i)
    sample = colnames(exp)[i]
    sample = gsub("embryonic.facial.prominence_","embryonic-facial-prominence_",sample)
    sample = gsub("neural.tube_","neural-tube_",sample)
    #
    data_normal = read.table(paste("./normal_window/mm10_gene_",sample,"_normal_state_TSS",length,
                                   "_count.txt",sep=""), row.names=1)
    data_dhs = read.table(paste("./dhs_window/mm10_gene_",sample,"_dhs_state_TSS",length,
                                 "_count.txt",sep=""), row.names=1)
    #
    model_normal = lm(log(exp[rownames(data_normal),11]+1e-5) ~ log(data_normal[,1]+1e-5)+
                        log(data_normal[,2]+1e-5)+log(data_normal[,3]+1e-5)+
                      log(data_normal[,4]+1e-5)+log(data_normal[,5]+1e-5)+log(data_normal[,6]+1e-5)+
                      log(data_normal[,7]+1e-5)+log(data_normal[,8]+1e-5)+log(data_normal[,9]+1e-5)+
                      log(data_normal[,10]+1e-5)+log(data_normal[,11]+1e-5)+log(data_normal[,12]+1e-5)+
                      log(data_normal[,13]+1e-5)+log(data_normal[,14]+1e-5)+log(data_normal[,15]+1e-5)+
                      log(data_normal[,16]+1e-5)+log(data_normal[,17]+1e-5)+log(data_normal[,18]+1e-5)+
                      log(data_normal[,19]+1e-5)+log(data_normal[,20]+1e-5)+log(data_normal[,21]+1e-5)+
                      log(data_normal[,22]+1e-5)+log(data_normal[,23]+1e-5)+log(data_normal[,24]+1e-5)+
                      log(data_normal[,25]+1e-5)+log(data_normal[,26]+1e-5)+log(data_normal[,27]+1e-5)+
                      log(data_normal[,28]+1e-5)+log(data_normal[,29]+1e-5)+log(data_normal[,30]+1e-5)+
                      log(data_normal[,31]+1e-5)+log(data_normal[,32]+1e-5)+log(data_normal[,33]+1e-5)+
                      log(data_normal[,34]+1e-5)+log(data_normal[,35]+1e-5)+log(data_normal[,36]+1e-5)+
                      log(data_normal[,37]+1e-5)+log(data_normal[,38]+1e-5)+log(data_normal[,39]+1e-5)+
                      log(data_normal[,40]+1e-5)+log(data_normal[,41]+1e-5)+log(data_normal[,42]+1e-5)+
                      log(data_normal[,43]+1e-5)+log(data_normal[,44]+1e-5))
    model_dhs = lm(log(exp[rownames(data_dhs),11]+1e-5) ~ log(data_dhs[,1]+1e-5)+
                     log(data_dhs[,2]+1e-5)+log(data_dhs[,3]+1e-5)+
                      log(data_dhs[,4]+1e-5)+log(data_dhs[,5]+1e-5)+log(data_dhs[,6]+1e-5)+
                      log(data_dhs[,7]+1e-5)+log(data_dhs[,8]+1e-5)+log(data_dhs[,9]+1e-5)+
                      log(data_dhs[,10]+1e-5)+log(data_dhs[,11]+1e-5)+log(data_dhs[,12]+1e-5)+
                      log(data_dhs[,13]+1e-5)+log(data_dhs[,14]+1e-5)+log(data_dhs[,15]+1e-5)+
                      log(data_dhs[,16]+1e-5)+log(data_dhs[,17]+1e-5)+log(data_dhs[,18]+1e-5)+
                      log(data_dhs[,19]+1e-5)+log(data_dhs[,20]+1e-5)+log(data_dhs[,21]+1e-5)+
                      log(data_dhs[,22]+1e-5)+log(data_dhs[,23]+1e-5)+log(data_dhs[,24]+1e-5)+
                      log(data_dhs[,25]+1e-5)+log(data_dhs[,26]+1e-5)+log(data_dhs[,27]+1e-5)+
                      log(data_dhs[,28]+1e-5)+log(data_dhs[,29]+1e-5)+log(data_dhs[,30]+1e-5)+
                      log(data_dhs[,31]+1e-5)+log(data_dhs[,32]+1e-5)+log(data_dhs[,33]+1e-5)+
                      log(data_dhs[,34]+1e-5)+log(data_dhs[,35]+1e-5)+log(data_dhs[,36]+1e-5)+
                      log(data_dhs[,37]+1e-5)+log(data_dhs[,38]+1e-5)+log(data_dhs[,39]+1e-5)+
                      log(data_dhs[,40]+1e-5)+log(data_dhs[,41]+1e-5)+log(data_dhs[,42]+1e-5)+
                      log(data_dhs[,43]+1e-5))
    #
    Rmatrix[1,i] = summary(model_normal)$adj.r.squared
    Rmatrix[2,i] = summary(model_dhs)$adj.r.squared
    colnames(Rmatrix)[i] = sample
  }
  return(Rmatrix)
}

## TSS up 200bp
# Rmatrix=get_within_celltype_rsquare(200)
# dat = t(Rmatrix)
# matrix = melt(dat)
# p = formatC(t.test(Rmatrix[1,],Rmatrix[2,], paired = TRUE)$p.value, format = "e", digits = 2)
# ggplot(matrix, aes(x=Var2,y=value,fill=Var2)) +
#   geom_boxplot(width=0.5) + theme_minimal() +
#   scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
#   ylab("adjusted R-squared") + xlab("TSup200bp") +
#   labs(subtitle = paste("paired t-test p-value=",p,sep="")) +
#   theme(title = element_text(face="bold", size=12),legend.title=element_blank())
# ggsave("TSSup200bp_adjRsqaure_boxplot.pdf")

## TSS±200
Rmatrix=get_within_celltype_rsquare(200)
dat = t(Rmatrix)
matrix = melt(dat)
p = formatC(t.test(Rmatrix[1,],Rmatrix[2,], paired = TRUE)$p.value, format = "e", digits = 2)
ggplot(matrix, aes(x=Var2,y=value,fill=Var2)) +
  geom_boxplot(width=0.5) + theme_minimal() +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  ylab("adjusted R-squared") + xlab("TSS±200bp") +
  labs(subtitle = paste("paired t-test p-value=",p,sep="")) +
  theme(title = element_text(face="bold", size=12),legend.title=element_blank())
ggsave("TSS200bp_adjRsqaure_boxplot.pdf")

## TSS±400
Rmatrix=get_within_celltype_rsquare(400)
dat = t(Rmatrix)
matrix = melt(dat)
p = formatC(t.test(Rmatrix[1,],Rmatrix[2,], paired = TRUE)$p.value, format = "e", digits = 2)
ggplot(matrix, aes(x=Var2,y=value,fill=Var2)) +
  geom_boxplot(width=0.5) + theme_minimal() +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  ylab("adjusted R-squared") + xlab("TSS±400bp") +
  labs(subtitle = paste("paired t-test p-value=",p,sep="")) +
  theme(title = element_text(face="bold", size=12),legend.title=element_blank())
ggsave("TSS400bp_adjRsqaure_boxplot.pdf")

## TSS±500
Rmatrix=get_within_celltype_rsquare(500)
dat = t(Rmatrix)
matrix = melt(dat)
p = formatC(t.test(Rmatrix[1,],Rmatrix[2,], paired = TRUE)$p.value, format = "e", digits = 2)
ggplot(matrix, aes(x=Var2,y=value,fill=Var2)) +
  geom_boxplot(width=0.5) + theme_minimal() +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  ylab("adjusted R-squared") + xlab("TSS±500bp") +
  labs(subtitle = paste("paired t-test p-value=",p,sep="")) +
  theme(title = element_text(face="bold", size=12),legend.title=element_blank())
ggsave("TSS500bp_adjRsqaure_boxplot.pdf")


## TSS±2kb
Rmatrix=get_within_celltype_rsquare(2000)
dat = t(Rmatrix)
matrix = melt(dat)
p = formatC(t.test(Rmatrix[1,],Rmatrix[2,], paired = TRUE)$p.value, format = "e", digits = 2)
ggplot(matrix, aes(x=Var2,y=value,fill=Var2)) +
  geom_boxplot(width=0.5) + theme_minimal() +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  ylab("adjusted R-squared") + xlab("TSS±2kb") +
  labs(subtitle = paste("paired t-test p-value=",p,sep="")) +
  theme(title = element_text(face="bold", size=12),legend.title=element_blank())
ggsave("TSS2kb_adjRsqaure_boxplot.pdf")


###################
# TSS±400bp, given 4*42 features
###################
Rmatrix_400 = matrix(0,ncol=66, nrow=2)
colnames(Rmatrix_400) = rep("a",66)
rownames(Rmatrix_400) = c("normal", "dhs")

for(i in 1:ncol(exp)){
  print(i)
  sample = colnames(exp)[i]
  sample = gsub("embryonic.facial.prominence_","embryonic-facial-prominence_",sample)
  sample = gsub("neural.tube_","neural-tube_",sample)
  #
  data_200bp_9 = read.table(paste("./normal_window/mm10_gene_",sample,
                                  "_normal_state_count_window_9.txt",sep=""), row.names=1)
  data_200bp_10 = read.table(paste("./normal_window/mm10_gene_",sample,
                                  "_normal_state_count_window_10.txt",sep=""), row.names=1)
  data_200bp_11 = read.table(paste("./normal_window/mm10_gene_",sample,
                                  "_normal_state_count_window_11.txt",sep=""), row.names=1)
  data_200bp_12 = read.table(paste("./normal_window/mm10_gene_",sample,
                                  "_normal_state_count_window_12.txt",sep=""), row.names=1)
  names_200bp = as.vector(intersect(intersect(intersect(rownames(data_200bp_9),rownames(data_200bp_10)),
                                              rownames(data_200bp_11)),rownames(data_200bp_12)))
  #
  data_dhs_9 = read.table(paste("./dhs_window/mm10_gene_",sample,
                                "_dhs_state_count_window_9.txt",sep=""), row.names=1)
  data_dhs_10 = read.table(paste("./dhs_window/mm10_gene_",sample,
                                "_dhs_state_count_window_10.txt",sep=""), row.names=1)
  data_dhs_11 = read.table(paste("./dhs_window/mm10_gene_",sample,
                                "_dhs_state_count_window_11.txt",sep=""), row.names=1)
  data_dhs_12 = read.table(paste("./dhs_window/mm10_gene_",sample,
                                "_dhs_state_count_window_12.txt",sep=""), row.names=1)
  names_dhs = as.vector(intersect(intersect(intersect(rownames(data_dhs_9),rownames(data_dhs_10)),
                                              rownames(data_dhs_11)),rownames(data_dhs_12)))
  #
  model_200bp = lm(log(exp[names_200bp,11]+1e-5) ~ log(data_200bp_9[names_200bp,2]+1e-5)+log(data_200bp_9[names_200bp,3]+1e-5)+
                      log(data_200bp_9[names_200bp,4]+1e-5)+log(data_200bp_9[names_200bp,5]+1e-5)+log(data_200bp_9[names_200bp,6]+1e-5)+
                      log(data_200bp_9[names_200bp,7]+1e-5)+log(data_200bp_9[names_200bp,8]+1e-5)+log(data_200bp_9[names_200bp,9]+1e-5)+
                      log(data_200bp_9[names_200bp,10]+1e-5)+log(data_200bp_9[names_200bp,11]+1e-5)+log(data_200bp_9[names_200bp,12]+1e-5)+
                      log(data_200bp_9[names_200bp,13]+1e-5)+log(data_200bp_9[names_200bp,14]+1e-5)+log(data_200bp_9[names_200bp,15]+1e-5)+
                      log(data_200bp_9[names_200bp,16]+1e-5)+log(data_200bp_9[names_200bp,17]+1e-5)+log(data_200bp_9[names_200bp,18]+1e-5)+
                      log(data_200bp_9[names_200bp,19]+1e-5)+log(data_200bp_9[names_200bp,20]+1e-5)+log(data_200bp_9[names_200bp,21]+1e-5)+
                      log(data_200bp_9[names_200bp,22]+1e-5)+log(data_200bp_9[names_200bp,23]+1e-5)+log(data_200bp_9[names_200bp,24]+1e-5)+
                      log(data_200bp_9[names_200bp,25]+1e-5)+log(data_200bp_9[names_200bp,26]+1e-5)+log(data_200bp_9[names_200bp,27]+1e-5)+
                      log(data_200bp_9[names_200bp,28]+1e-5)+log(data_200bp_9[names_200bp,29]+1e-5)+log(data_200bp_9[names_200bp,30]+1e-5)+
                      log(data_200bp_9[names_200bp,31]+1e-5)+log(data_200bp_9[names_200bp,32]+1e-5)+log(data_200bp_9[names_200bp,33]+1e-5)+
                      log(data_200bp_9[names_200bp,34]+1e-5)+log(data_200bp_9[names_200bp,35]+1e-5)+log(data_200bp_9[names_200bp,36]+1e-5)+
                      log(data_200bp_9[names_200bp,37]+1e-5)+log(data_200bp_9[names_200bp,38]+1e-5)+log(data_200bp_9[names_200bp,39]+1e-5)+
                      log(data_200bp_9[names_200bp,40]+1e-5)+log(data_200bp_9[names_200bp,41]+1e-5)+log(data_200bp_9[names_200bp,42]+1e-5)+
                      log(data_200bp_9[names_200bp,43]+1e-5)+log(data_200bp_9[names_200bp,44]+1e-5)+log(data_200bp_9[names_200bp,45]+1e-5)+
                     log(data_200bp_10[names_200bp,2]+1e-5)+log(data_200bp_10[names_200bp,3]+1e-5)+
                     log(data_200bp_10[names_200bp,4]+1e-5)+log(data_200bp_10[names_200bp,5]+1e-5)+log(data_200bp_10[names_200bp,6]+1e-5)+
                     log(data_200bp_10[names_200bp,7]+1e-5)+log(data_200bp_10[names_200bp,8]+1e-5)+log(data_200bp_10[names_200bp,9]+1e-5)+
                     log(data_200bp_10[names_200bp,10]+1e-5)+log(data_200bp_10[names_200bp,11]+1e-5)+log(data_200bp_10[names_200bp,12]+1e-5)+
                     log(data_200bp_10[names_200bp,13]+1e-5)+log(data_200bp_10[names_200bp,14]+1e-5)+log(data_200bp_10[names_200bp,15]+1e-5)+
                     log(data_200bp_10[names_200bp,16]+1e-5)+log(data_200bp_10[names_200bp,17]+1e-5)+log(data_200bp_10[names_200bp,18]+1e-5)+
                     log(data_200bp_10[names_200bp,19]+1e-5)+log(data_200bp_10[names_200bp,20]+1e-5)+log(data_200bp_10[names_200bp,21]+1e-5)+
                     log(data_200bp_10[names_200bp,22]+1e-5)+log(data_200bp_10[names_200bp,23]+1e-5)+log(data_200bp_10[names_200bp,24]+1e-5)+
                     log(data_200bp_10[names_200bp,25]+1e-5)+log(data_200bp_10[names_200bp,26]+1e-5)+log(data_200bp_10[names_200bp,27]+1e-5)+
                     log(data_200bp_10[names_200bp,28]+1e-5)+log(data_200bp_10[names_200bp,29]+1e-5)+log(data_200bp_10[names_200bp,30]+1e-5)+
                     log(data_200bp_10[names_200bp,31]+1e-5)+log(data_200bp_10[names_200bp,32]+1e-5)+log(data_200bp_10[names_200bp,33]+1e-5)+
                     log(data_200bp_10[names_200bp,34]+1e-5)+log(data_200bp_10[names_200bp,35]+1e-5)+log(data_200bp_10[names_200bp,36]+1e-5)+
                     log(data_200bp_10[names_200bp,37]+1e-5)+log(data_200bp_10[names_200bp,38]+1e-5)+log(data_200bp_10[names_200bp,39]+1e-5)+
                     log(data_200bp_10[names_200bp,40]+1e-5)+log(data_200bp_10[names_200bp,41]+1e-5)+log(data_200bp_10[names_200bp,42]+1e-5)+
                     log(data_200bp_10[names_200bp,43]+1e-5)+log(data_200bp_10[names_200bp,44]+1e-5)+log(data_200bp_10[names_200bp,45]+1e-5)+
                     log(data_200bp_11[names_200bp,2]+1e-5)+log(data_200bp_11[names_200bp,3]+1e-5)+
                     log(data_200bp_11[names_200bp,4]+1e-5)+log(data_200bp_11[names_200bp,5]+1e-5)+log(data_200bp_11[names_200bp,6]+1e-5)+
                     log(data_200bp_11[names_200bp,7]+1e-5)+log(data_200bp_11[names_200bp,8]+1e-5)+log(data_200bp_11[names_200bp,9]+1e-5)+
                     log(data_200bp_11[names_200bp,10]+1e-5)+log(data_200bp_11[names_200bp,11]+1e-5)+log(data_200bp_11[names_200bp,12]+1e-5)+
                     log(data_200bp_11[names_200bp,13]+1e-5)+log(data_200bp_11[names_200bp,14]+1e-5)+log(data_200bp_11[names_200bp,15]+1e-5)+
                     log(data_200bp_11[names_200bp,16]+1e-5)+log(data_200bp_11[names_200bp,17]+1e-5)+log(data_200bp_11[names_200bp,18]+1e-5)+
                     log(data_200bp_11[names_200bp,19]+1e-5)+log(data_200bp_11[names_200bp,20]+1e-5)+log(data_200bp_11[names_200bp,21]+1e-5)+
                     log(data_200bp_11[names_200bp,22]+1e-5)+log(data_200bp_11[names_200bp,23]+1e-5)+log(data_200bp_11[names_200bp,24]+1e-5)+
                     log(data_200bp_11[names_200bp,25]+1e-5)+log(data_200bp_11[names_200bp,26]+1e-5)+log(data_200bp_11[names_200bp,27]+1e-5)+
                     log(data_200bp_11[names_200bp,28]+1e-5)+log(data_200bp_11[names_200bp,29]+1e-5)+log(data_200bp_11[names_200bp,30]+1e-5)+
                     log(data_200bp_11[names_200bp,31]+1e-5)+log(data_200bp_11[names_200bp,32]+1e-5)+log(data_200bp_11[names_200bp,33]+1e-5)+
                     log(data_200bp_11[names_200bp,34]+1e-5)+log(data_200bp_11[names_200bp,35]+1e-5)+log(data_200bp_11[names_200bp,36]+1e-5)+
                     log(data_200bp_11[names_200bp,37]+1e-5)+log(data_200bp_11[names_200bp,38]+1e-5)+log(data_200bp_11[names_200bp,39]+1e-5)+
                     log(data_200bp_11[names_200bp,40]+1e-5)+log(data_200bp_11[names_200bp,41]+1e-5)+log(data_200bp_11[names_200bp,42]+1e-5)+
                     log(data_200bp_11[names_200bp,43]+1e-5)+log(data_200bp_11[names_200bp,44]+1e-5)+log(data_200bp_11[names_200bp,45]+1e-5)+
                     log(data_200bp_12[names_200bp,2]+1e-5)+log(data_200bp_12[names_200bp,3]+1e-5)+
                     log(data_200bp_12[names_200bp,4]+1e-5)+log(data_200bp_12[names_200bp,5]+1e-5)+log(data_200bp_12[names_200bp,6]+1e-5)+
                     log(data_200bp_12[names_200bp,7]+1e-5)+log(data_200bp_12[names_200bp,8]+1e-5)+log(data_200bp_12[names_200bp,9]+1e-5)+
                     log(data_200bp_12[names_200bp,10]+1e-5)+log(data_200bp_12[names_200bp,11]+1e-5)+log(data_200bp_12[names_200bp,12]+1e-5)+
                     log(data_200bp_12[names_200bp,13]+1e-5)+log(data_200bp_12[names_200bp,14]+1e-5)+log(data_200bp_12[names_200bp,15]+1e-5)+
                     log(data_200bp_12[names_200bp,16]+1e-5)+log(data_200bp_12[names_200bp,17]+1e-5)+log(data_200bp_12[names_200bp,18]+1e-5)+
                     log(data_200bp_12[names_200bp,19]+1e-5)+log(data_200bp_12[names_200bp,20]+1e-5)+log(data_200bp_12[names_200bp,21]+1e-5)+
                     log(data_200bp_12[names_200bp,22]+1e-5)+log(data_200bp_12[names_200bp,23]+1e-5)+log(data_200bp_12[names_200bp,24]+1e-5)+
                     log(data_200bp_12[names_200bp,25]+1e-5)+log(data_200bp_12[names_200bp,26]+1e-5)+log(data_200bp_12[names_200bp,27]+1e-5)+
                     log(data_200bp_12[names_200bp,28]+1e-5)+log(data_200bp_12[names_200bp,29]+1e-5)+log(data_200bp_12[names_200bp,30]+1e-5)+
                     log(data_200bp_12[names_200bp,31]+1e-5)+log(data_200bp_12[names_200bp,32]+1e-5)+log(data_200bp_12[names_200bp,33]+1e-5)+
                     log(data_200bp_12[names_200bp,34]+1e-5)+log(data_200bp_12[names_200bp,35]+1e-5)+log(data_200bp_12[names_200bp,36]+1e-5)+
                     log(data_200bp_12[names_200bp,37]+1e-5)+log(data_200bp_12[names_200bp,38]+1e-5)+log(data_200bp_12[names_200bp,39]+1e-5)+
                     log(data_200bp_12[names_200bp,40]+1e-5)+log(data_200bp_12[names_200bp,41]+1e-5)+log(data_200bp_12[names_200bp,42]+1e-5)+
                     log(data_200bp_12[names_200bp,43]+1e-5)+log(data_200bp_12[names_200bp,44]+1e-5)+log(data_200bp_12[names_200bp,45]+1e-5))
  
  model_dhs = lm(log(exp[names_dhs,11]+1e-5) ~ log(data_dhs_9[names_dhs,2]+1e-5)+log(data_dhs_9[names_dhs,3]+1e-5)+
                   log(data_dhs_9[names_dhs,4]+1e-5)+log(data_dhs_9[names_dhs,5]+1e-5)+log(data_dhs_9[names_dhs,6]+1e-5)+
                   log(data_dhs_9[names_dhs,7]+1e-5)+log(data_dhs_9[names_dhs,8]+1e-5)+log(data_dhs_9[names_dhs,9]+1e-5)+
                   log(data_dhs_9[names_dhs,10]+1e-5)+log(data_dhs_9[names_dhs,11]+1e-5)+log(data_dhs_9[names_dhs,12]+1e-5)+
                   log(data_dhs_9[names_dhs,13]+1e-5)+log(data_dhs_9[names_dhs,14]+1e-5)+log(data_dhs_9[names_dhs,15]+1e-5)+
                   log(data_dhs_9[names_dhs,16]+1e-5)+log(data_dhs_9[names_dhs,17]+1e-5)+log(data_dhs_9[names_dhs,18]+1e-5)+
                   log(data_dhs_9[names_dhs,19]+1e-5)+log(data_dhs_9[names_dhs,20]+1e-5)+log(data_dhs_9[names_dhs,21]+1e-5)+
                   log(data_dhs_9[names_dhs,22]+1e-5)+log(data_dhs_9[names_dhs,23]+1e-5)+log(data_dhs_9[names_dhs,24]+1e-5)+
                   log(data_dhs_9[names_dhs,25]+1e-5)+log(data_dhs_9[names_dhs,26]+1e-5)+log(data_dhs_9[names_dhs,27]+1e-5)+
                   log(data_dhs_9[names_dhs,28]+1e-5)+log(data_dhs_9[names_dhs,29]+1e-5)+log(data_dhs_9[names_dhs,30]+1e-5)+
                   log(data_dhs_9[names_dhs,31]+1e-5)+log(data_dhs_9[names_dhs,32]+1e-5)+log(data_dhs_9[names_dhs,33]+1e-5)+
                   log(data_dhs_9[names_dhs,34]+1e-5)+log(data_dhs_9[names_dhs,35]+1e-5)+log(data_dhs_9[names_dhs,36]+1e-5)+
                   log(data_dhs_9[names_dhs,37]+1e-5)+log(data_dhs_9[names_dhs,38]+1e-5)+log(data_dhs_9[names_dhs,39]+1e-5)+
                   log(data_dhs_9[names_dhs,40]+1e-5)+log(data_dhs_9[names_dhs,41]+1e-5)+log(data_dhs_9[names_dhs,42]+1e-5)+
                   log(data_dhs_9[names_dhs,43]+1e-5)+log(data_dhs_9[names_dhs,44]+1e-5)+
                   log(data_dhs_10[names_dhs,2]+1e-5)+log(data_dhs_10[names_dhs,3]+1e-5)+
                   log(data_dhs_10[names_dhs,4]+1e-5)+log(data_dhs_10[names_dhs,5]+1e-5)+log(data_dhs_10[names_dhs,6]+1e-5)+
                   log(data_dhs_10[names_dhs,7]+1e-5)+log(data_dhs_10[names_dhs,8]+1e-5)+log(data_dhs_10[names_dhs,9]+1e-5)+
                   log(data_dhs_10[names_dhs,10]+1e-5)+log(data_dhs_10[names_dhs,11]+1e-5)+log(data_dhs_10[names_dhs,12]+1e-5)+
                   log(data_dhs_10[names_dhs,13]+1e-5)+log(data_dhs_10[names_dhs,14]+1e-5)+log(data_dhs_10[names_dhs,15]+1e-5)+
                   log(data_dhs_10[names_dhs,16]+1e-5)+log(data_dhs_10[names_dhs,17]+1e-5)+log(data_dhs_10[names_dhs,18]+1e-5)+
                   log(data_dhs_10[names_dhs,19]+1e-5)+log(data_dhs_10[names_dhs,20]+1e-5)+log(data_dhs_10[names_dhs,21]+1e-5)+
                   log(data_dhs_10[names_dhs,22]+1e-5)+log(data_dhs_10[names_dhs,23]+1e-5)+log(data_dhs_10[names_dhs,24]+1e-5)+
                   log(data_dhs_10[names_dhs,25]+1e-5)+log(data_dhs_10[names_dhs,26]+1e-5)+log(data_dhs_10[names_dhs,27]+1e-5)+
                   log(data_dhs_10[names_dhs,28]+1e-5)+log(data_dhs_10[names_dhs,29]+1e-5)+log(data_dhs_10[names_dhs,30]+1e-5)+
                   log(data_dhs_10[names_dhs,31]+1e-5)+log(data_dhs_10[names_dhs,32]+1e-5)+log(data_dhs_10[names_dhs,33]+1e-5)+
                   log(data_dhs_10[names_dhs,34]+1e-5)+log(data_dhs_10[names_dhs,35]+1e-5)+log(data_dhs_10[names_dhs,36]+1e-5)+
                   log(data_dhs_10[names_dhs,37]+1e-5)+log(data_dhs_10[names_dhs,38]+1e-5)+log(data_dhs_10[names_dhs,39]+1e-5)+
                   log(data_dhs_10[names_dhs,40]+1e-5)+log(data_dhs_10[names_dhs,41]+1e-5)+log(data_dhs_10[names_dhs,42]+1e-5)+
                   log(data_dhs_10[names_dhs,43]+1e-5)+log(data_dhs_10[names_dhs,44]+1e-5)+
                   log(data_dhs_11[names_dhs,2]+1e-5)+log(data_dhs_11[names_dhs,3]+1e-5)+
                   log(data_dhs_11[names_dhs,4]+1e-5)+log(data_dhs_11[names_dhs,5]+1e-5)+log(data_dhs_11[names_dhs,6]+1e-5)+
                   log(data_dhs_11[names_dhs,7]+1e-5)+log(data_dhs_11[names_dhs,8]+1e-5)+log(data_dhs_11[names_dhs,9]+1e-5)+
                   log(data_dhs_11[names_dhs,10]+1e-5)+log(data_dhs_11[names_dhs,11]+1e-5)+log(data_dhs_11[names_dhs,12]+1e-5)+
                   log(data_dhs_11[names_dhs,13]+1e-5)+log(data_dhs_11[names_dhs,14]+1e-5)+log(data_dhs_11[names_dhs,15]+1e-5)+
                   log(data_dhs_11[names_dhs,16]+1e-5)+log(data_dhs_11[names_dhs,17]+1e-5)+log(data_dhs_11[names_dhs,18]+1e-5)+
                   log(data_dhs_11[names_dhs,19]+1e-5)+log(data_dhs_11[names_dhs,20]+1e-5)+log(data_dhs_11[names_dhs,21]+1e-5)+
                   log(data_dhs_11[names_dhs,22]+1e-5)+log(data_dhs_11[names_dhs,23]+1e-5)+log(data_dhs_11[names_dhs,24]+1e-5)+
                   log(data_dhs_11[names_dhs,25]+1e-5)+log(data_dhs_11[names_dhs,26]+1e-5)+log(data_dhs_11[names_dhs,27]+1e-5)+
                   log(data_dhs_11[names_dhs,28]+1e-5)+log(data_dhs_11[names_dhs,29]+1e-5)+log(data_dhs_11[names_dhs,30]+1e-5)+
                   log(data_dhs_11[names_dhs,31]+1e-5)+log(data_dhs_11[names_dhs,32]+1e-5)+log(data_dhs_11[names_dhs,33]+1e-5)+
                   log(data_dhs_11[names_dhs,34]+1e-5)+log(data_dhs_11[names_dhs,35]+1e-5)+log(data_dhs_11[names_dhs,36]+1e-5)+
                   log(data_dhs_11[names_dhs,37]+1e-5)+log(data_dhs_11[names_dhs,38]+1e-5)+log(data_dhs_11[names_dhs,39]+1e-5)+
                   log(data_dhs_11[names_dhs,40]+1e-5)+log(data_dhs_11[names_dhs,41]+1e-5)+log(data_dhs_11[names_dhs,42]+1e-5)+
                   log(data_dhs_11[names_dhs,43]+1e-5)+log(data_dhs_11[names_dhs,44]+1e-5)+
                   log(data_dhs_12[names_dhs,2]+1e-5)+log(data_dhs_12[names_dhs,3]+1e-5)+
                   log(data_dhs_12[names_dhs,4]+1e-5)+log(data_dhs_12[names_dhs,5]+1e-5)+log(data_dhs_12[names_dhs,6]+1e-5)+
                   log(data_dhs_12[names_dhs,7]+1e-5)+log(data_dhs_12[names_dhs,8]+1e-5)+log(data_dhs_12[names_dhs,9]+1e-5)+
                   log(data_dhs_12[names_dhs,10]+1e-5)+log(data_dhs_12[names_dhs,11]+1e-5)+log(data_dhs_12[names_dhs,12]+1e-5)+
                   log(data_dhs_12[names_dhs,13]+1e-5)+log(data_dhs_12[names_dhs,14]+1e-5)+log(data_dhs_12[names_dhs,15]+1e-5)+
                   log(data_dhs_12[names_dhs,16]+1e-5)+log(data_dhs_12[names_dhs,17]+1e-5)+log(data_dhs_12[names_dhs,18]+1e-5)+
                   log(data_dhs_12[names_dhs,19]+1e-5)+log(data_dhs_12[names_dhs,20]+1e-5)+log(data_dhs_12[names_dhs,21]+1e-5)+
                   log(data_dhs_12[names_dhs,22]+1e-5)+log(data_dhs_12[names_dhs,23]+1e-5)+log(data_dhs_12[names_dhs,24]+1e-5)+
                   log(data_dhs_12[names_dhs,25]+1e-5)+log(data_dhs_12[names_dhs,26]+1e-5)+log(data_dhs_12[names_dhs,27]+1e-5)+
                   log(data_dhs_12[names_dhs,28]+1e-5)+log(data_dhs_12[names_dhs,29]+1e-5)+log(data_dhs_12[names_dhs,30]+1e-5)+
                   log(data_dhs_12[names_dhs,31]+1e-5)+log(data_dhs_12[names_dhs,32]+1e-5)+log(data_dhs_12[names_dhs,33]+1e-5)+
                   log(data_dhs_12[names_dhs,34]+1e-5)+log(data_dhs_12[names_dhs,35]+1e-5)+log(data_dhs_12[names_dhs,36]+1e-5)+
                   log(data_dhs_12[names_dhs,37]+1e-5)+log(data_dhs_12[names_dhs,38]+1e-5)+log(data_dhs_12[names_dhs,39]+1e-5)+
                   log(data_dhs_12[names_dhs,40]+1e-5)+log(data_dhs_12[names_dhs,41]+1e-5)+log(data_dhs_12[names_dhs,42]+1e-5)+
                   log(data_dhs_12[names_dhs,43]+1e-5)+log(data_dhs_12[names_dhs,44]+1e-5))
  #
  Rmatrix_400[1,i] = summary(model_200bp)$adj.r.squared
  Rmatrix_400[2,i] = summary(model_dhs)$adj.r.squared
  colnames(Rmatrix)[i] = sample
}


dat = t(Rmatrix_400)
matrix = melt(dat)
p = formatC(t.test(Rmatrix_400[1,],Rmatrix_400[2,], paired = TRUE)$p.value, format = "e", digits = 2)
ggplot(matrix, aes(x=Var2,y=value,fill=Var2, group=Var2)) +
  geom_boxplot(width=0.5) + theme_minimal() +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  ylab("adjusted R-squared") + xlab("TSS±400bp") +
  labs(title="given seperate features",
       subtitle = paste("paired t-test p-value=",p,sep="")) +
  theme(title = element_text(face="bold", size=12),legend.title=element_blank())
ggsave("TSS400bp_seperate_feature_adjRsqaure_boxplot.pdf")



###################
# rsquare average for TSS±400bp, 500bp, 2kb
###################

# TSS±400bp
a1 = R1_matrix[9:12,]
b1 = apply(a1,2,mean)
a2 = R2_matrix[9:12,]
b2 = apply(a2,2,mean)
p = formatC(t.test(b1, b2, paired = TRUE)$p.value, format = "e", digits = 2)

matrix = data.frame(rbind(transform(b1, type="normal"),
               transform(b2, type="dhs")))
colnames(matrix) = c("value","type")
ggplot(matrix,aes(x=type, y=value, fill = type)) +
  geom_boxplot(width=0.3) +
  theme_minimal() +
  scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  ylab("adjusted R-squared") + xlab("TSS±400bp") +
  labs(title ="average of bins' r-square",
       subtitle = paste("paired t-test p-value=",p,sep="")) +
  theme(title = element_text(face="bold", size=12),legend.title=element_blank())
ggsave("TSS400bp_adjRsqaure_average_boxplot.pdf")


# TSS±2kb
a1 = R1_matrix[1:20,]
b1 = apply(a1,2,mean)
a2 = R2_matrix[1:20,]
b2 = apply(a2,2,mean)
p = formatC(t.test(b1, b2, paired = TRUE)$p.value, format = "e", digits = 2)

matrix = data.frame(rbind(transform(b1, type="normal"),
                          transform(b2, type="dhs")))
colnames(matrix) = c("value","type")
ggplot(matrix,aes(x=type, y=value, fill = type)) +
  geom_boxplot(width=0.3) +
  theme_minimal() +
  scale_fill_manual(values = c("#00bfc4","#f8766d")) +
  ylab("adjusted R-squared") + xlab("TSS±2kb") +
  labs(title ="average of bins' r-square",
       subtitle = paste("paired t-test p-value=",p,sep="")) +
  theme(title = element_text(face="bold", size=12),legend.title=element_blank())
ggsave("TSS2kb_adjRsqaure_average_boxplot.pdf")
###################
# scatter plot, average
###################
r2_ave_200bp = apply(R1_matrix,2,mean)
r2_ave_dhs = apply(R2_matrix,2,mean)
#
r2_ave_200bp = transform(r2_ave_200bp, type="200bp",
                         t(data.frame(strsplit(names(r2_ave_200bp), split="_"))))
colnames(r2_ave_200bp) = c()
#
r2_ave = cbind(r2_ave_200bp, r2_ave_dhs[rownames(r2_ave_200bp)])
colnames(r2_ave) = c("rsquare_200bp", "type_200bp", "celltype_200bp", "timepoint_200bp", "rsquare_dhs")
r2_ave$celltype_200bp = factor(r2_ave$celltype_200bp)
r2_ave$timepoint_200bp = factor(r2_ave$timepoint_200bp)
#
ggplot(r2_ave, aes(x = rsquare_200bp, y = rsquare_dhs)) + 
  geom_point(size=3, aes(shape=timepoint_200bp, col = celltype_200bp)) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, col = "red", linetype = "dashed") + 
  xlab("adjusted R-square of 200bp bins") + ylab("adjusted R-square of dhs bins") +
  labs(title = "mean of adjusted R-square of 20 windows", 
       shape = "time point", colour = "tissue") +
  theme(title = element_text(size=12, face = "bold")) +
  scale_shape_manual(values = 12:18) +
  coord_cartesian(xlim=c(0.37,0.45), ylim=c(0.37,0.45))
ggsave("mean_rsquare_each_sample_scatter.pdf")


###################
# scatter plot, average 4 points (±400bp)
###################
r2_ave_200bp = apply(R1_matrix[9:12,],2,mean)
r2_ave_dhs = apply(R2_matrix[9:12,],2,mean)
#
r2_ave_200bp = transform(r2_ave_200bp, type="200bp",
                         t(data.frame(strsplit(names(r2_ave_200bp), split="_"))))
colnames(r2_ave_200bp) = c()
#
r2_ave = cbind(r2_ave_200bp, r2_ave_dhs[rownames(r2_ave_200bp)])
colnames(r2_ave) = c("rsquare_200bp", "type_200bp", "celltype_200bp", "timepoint_200bp", "rsquare_dhs")
r2_ave$celltype_200bp = factor(r2_ave$celltype_200bp)
r2_ave$timepoint_200bp = factor(r2_ave$timepoint_200bp)
#
ggplot(r2_ave, aes(x = rsquare_200bp, y = rsquare_dhs)) + 
  geom_point(size=3, aes(shape=timepoint_200bp, col = celltype_200bp)) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, col = "red", linetype = "dashed") + 
  xlab("adjusted R-square of 200bp bins") + ylab("adjusted R-square of dhs bins") +
  labs(title = "mean of adjusted R-square of 4 windows (±400bp)", 
       shape = "time point", colour = "tissue") +
  theme(title = element_text(size=12, face = "bold")) +
  scale_shape_manual(values = 12:18) +
  coord_cartesian(xlim=c(0.54,0.6), ylim=c(0.54,0.6))
ggsave("mean_rsquare_4points_each_sample_scatter.pdf")


###################
# for state proportion2
###################
# this one using exact 20 bins
# not fair

pdf("within_gene_expression_regression2.pdf")
for(i in 1:ncol(exp)){
  sample = colnames(exp)[i]
  sample = gsub("embryonic.facial.prominence_","embryonic-facial-prominence_",sample)
  sample = gsub("neural.tube_","neural-tube_",sample)
  
  R1=NULL
  R2=NULL
  for(j in 1:20){
    data1 = read.table(paste("./normal_window2/mm10_gene_",sample,"_normal_state_count_window_",
                             as.character(j),".txt",sep=""))
    data2 = read.table(paste("./dhs_window2/mm10_gene_",sample,"_dhs_state_count_window_",
                             as.character(j),".txt",sep=""), row.names=1)
    #
    r1 = NULL
    r2 = NULL
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
                  log(data1[,43]+1e-5))
    model2 = lm(log(exp[rownames(data2),11]+1e-5) ~ log(data2[,2]+1e-5)+log(data2[,3]+1e-5)+
                  log(data2[,4]+1e-5)+log(data2[,5]+1e-5)+log(data2[,6]+1e-5)+
                  log(data2[,7]+1e-5)+log(data2[,8]+1e-5)+log(data2[,9]+1e-5)+
                  log(data2[,10]+1e-5)+log(data2[,11]+1e-5)+log(data2[,12]+1e-5)+
                  log(data2[,13]+1e-5)+log(data2[,14]+1e-5)+log(data2[,15]+1e-5)+
                  log(data2[,16]+1e-5)+log(data2[,17]+1e-5)+log(data2[,18]+1e-5)+
                  log(data2[,19]+1e-5)+log(data2[,20]+1e-5)+log(data2[,21]+1e-5)+
                  log(data2[,22]+1e-5)+log(data2[,23]+1e-5)+log(data2[,24]+1e-5)+
                  log(data2[,25]+1e-5)+log(data2[,26]+1e-5)+log(data2[,27]+1e-5)+
                  log(data2[,28]+1e-5)+log(data2[,29]+1e-5)+log(data2[,30]+1e-5)+
                  log(data2[,31]+1e-5)+log(data2[,32]+1e-5)+log(data2[,33]+1e-5)+
                  log(data2[,34]+1e-5)+log(data2[,35]+1e-5)+log(data2[,36]+1e-5)+
                  log(data2[,37]+1e-5)+log(data2[,38]+1e-5)+log(data2[,39]+1e-5)+
                  log(data2[,40]+1e-5)+log(data2[,41]+1e-5)+log(data2[,42]+1e-5))
    r1 = summary(model1)$adj.r.squared
    r2 = summary(model2)$adj.r.squared
    R1[j] = round(r1,4)
    R2[j] = round(r2,4)
  }
  matrix = data.frame(rbind(cbind(1:20,R1,"normal_bins"),
                            cbind(1:20,R2,"dhs_bins")))
  colnames(matrix) = c("window","r", "type")
  matrix$window = factor(matrix$window, levels = 1:20)
  
  p = ggplot(matrix, aes(x=window, y=r, group=type)) + 
    geom_line(aes(col=type),size=1.5) +
    theme_minimal() + geom_point(size=1) +
    xlab("") + ylab("adjusted R-squared") +
    theme(title=element_text(face="bold", size=12),
          axis.text.x = element_blank()) + 
    labs(title=sample) +
    scale_x_discrete(name="TSS±2kb", breaks=c(1,5,10,15,20)) +
    theme(legend.position = c(0.15,0.85), 
          legend.text = element_text(face="bold", size=12))
  print(p)
}
dev.off()
#######################