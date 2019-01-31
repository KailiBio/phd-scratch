
# -- Kaili
# This script is for making aggregation within cell type.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/regression/")

exp = read.table("mm10_RNA_protein-coding_tpm_matrix_matched.txt", 
                 row.names = 1, header=TRUE)

pdf("within_gene_expression_regression.pdf")
for(i in 1:ncol(exp)){
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
    model1 = lm(exp[rownames(data1),11] ~ log(data1[,2]+1e-5)+log(data1[,3]+1e-5)+
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
                  log(data1[,37]+1e-5)+log(data1[,38]+1e-5))
    model2 = lm(exp[rownames(data2),11] ~ log(data2[,2]+1e-5)+log(data2[,3]+1e-5)+
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
                  log(data2[,37]+1e-5)+log(data2[,38]+1e-5))
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
    geom_line(aes(linetype=type)) +
    theme_minimal() + geom_point(size=0.5) +
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
