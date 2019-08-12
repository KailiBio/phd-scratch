
# -- Kaili
# This script is for calculating the state proportion in each samples.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/state_file/")

file = "ctcf_samples.chr1.state"
file = "DHS_v3_100-400bp.chr1.state"



colnames = c("lung_14.5", "liver_14.5", "stomach_0", "midbrain_0", "kidney_0", "liver_0",
             "intestine_0", "lung_0", "heart_0", "hindbrain_0", "forebrain_0")

data = read.table(file, row.names = 1, skip=1)

#############
n = ncol(data)-4
total = nrow(data)
matrix = matrix(nrow = n, ncol = length(table(data[,4])))
colnames(matrix) = names(table(data[,4]))
rownames(matrix) = colnames
for(i in 1:n){
  matrix[i,] = as.vector(table(data[,(i+3)]))/total
}


library(RColorBrewer)

pdf("ss.pdf")
split.screen(c(2,1))
screen(n=1)
barplot(matrix[,1:12], beside = TRUE, col = brewer.pal(11,"Set3"))
legend("topright", legend=colnames, col = brewer.pal(11,"Set3"), pch=15, bty="n", 
       yjust=0.6, lwd=1.5, y.intersp = 0.8)
screen(n=2)
barplot(matrix[,13:23], beside = TRUE, col = brewer.pal(11,"Set3"))
screen(n=3)
barplot(matrix[,24:34], beside = TRUE, col = brewer.pal(11,"Set3"))
screen(n=4)
barplot(matrix[,35:44], beside = TRUE, col = brewer.pal(11,"Set3"))
close.screen(all=TRUE)
dev.off()


barplot(matrix[,1:4], beside = TRUE)
barplot(matrix[,5:8], beside = TRUE)
barplot(matrix[,9:12], beside = TRUE)
barplot(matrix[,13:16], beside = TRUE)
barplot(matrix[,17:20], beside = TRUE)
barplot(matrix[,21:24], beside = TRUE)
barplot(matrix[,25:28], beside = TRUE)
barplot(matrix[,29:32], beside = TRUE)
barplot(matrix[,33:36], beside = TRUE)
barplot(matrix[,36:38], beside = TRUE)

################################







