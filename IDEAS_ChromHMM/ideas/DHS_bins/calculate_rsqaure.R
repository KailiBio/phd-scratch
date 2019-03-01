
# -- Kaili
# This script is for calculating r-square.

args = commandArgs(trailingOnly=TRUE)
file = args[1]
type = args[2]
gene_name = args[3]

exp = read.table("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/mm10_RNA_protein-coding_tpm_matrix_matched.txt", row.names = 1,
                 header=TRUE)

data = read.table(file)
gene = as.character(data[1,1])
y = t(exp[gene_name,])
if(type=="dhs"){
  model = lm(log(y+1e-5) ~ log(data[,2]+1e-5)+log(data[,3]+1e-5)+log(data[,4]+1e-5)+
               log(data[,5]+1e-5)+log(data[,6]+1e-5)+
               log(data[,7]+1e-5)+log(data[,8]+1e-5)+log(data[,9]+1e-5)+
               log(data[,10]+1e-5)+log(data[,11]+1e-5)+log(data[,12]+1e-5)+
               log(data[,13]+1e-5)+log(data[,14]+1e-5)+log(data[,15]+1e-5)+
               log(data[,16]+1e-5)+log(data[,17]+1e-5)+log(data[,18]+1e-5)+
               log(data[,19]+1e-5)+log(data[,20]+1e-5)+log(data[,21]+1e-5)+
               log(data[,22]+1e-5)+log(data[,23]+1e-5)+log(data[,24]+1e-5)+
               log(data[,25]+1e-5)+log(data[,26]+1e-5)+log(data[,27]+1e-5)+
               log(data[,28]+1e-5)+log(data[,29]+1e-5)+log(data[,30]+1e-5)+
               log(data[,31]+1e-5)+log(data[,32]+1e-5)+log(data[,33]+1e-5)+
               log(data[,34]+1e-5)+log(data[,35]+1e-5)+log(data[,36]+1e-5)+
               log(data[,37]+1e-5)+log(data[,38]+1e-5)+log(data[,39]+1e-5)+
               log(data[,40]+1e-5)+log(data[,41]+1e-5)+log(data[,42]+1e-5)+
               log(data[,43]+1e-5)+log(data[,44]+1e-5))
}else{
  model = lm(log(y+1e-5) ~ log(data[,2]+1e-5)+log(data[,3]+1e-5)+log(data[,4]+1e-5)+
               log(data[,5]+1e-5)+log(data[,6]+1e-5)+
               log(data[,7]+1e-5)+log(data[,8]+1e-5)+log(data[,9]+1e-5)+
               log(data[,10]+1e-5)+log(data[,11]+1e-5)+log(data[,12]+1e-5)+
               log(data[,13]+1e-5)+log(data[,14]+1e-5)+log(data[,15]+1e-5)+
               log(data[,16]+1e-5)+log(data[,17]+1e-5)+log(data[,18]+1e-5)+
               log(data[,19]+1e-5)+log(data[,20]+1e-5)+log(data[,21]+1e-5)+
               log(data[,22]+1e-5)+log(data[,23]+1e-5)+log(data[,24]+1e-5)+
               log(data[,25]+1e-5)+log(data[,26]+1e-5)+log(data[,27]+1e-5)+
               log(data[,28]+1e-5)+log(data[,29]+1e-5)+log(data[,30]+1e-5)+
               log(data[,31]+1e-5)+log(data[,32]+1e-5)+log(data[,33]+1e-5)+
               log(data[,34]+1e-5)+log(data[,35]+1e-5)+log(data[,36]+1e-5)+
               log(data[,37]+1e-5)+log(data[,38]+1e-5)+log(data[,39]+1e-5)+
               log(data[,40]+1e-5)+log(data[,41]+1e-5)+log(data[,42]+1e-5)+
               log(data[,43]+1e-5)+log(data[,44]+1e-5)+log(data[,45]+1e-5))
}
r = round(summary(model)$adj.r.squared,4)
print(r)

