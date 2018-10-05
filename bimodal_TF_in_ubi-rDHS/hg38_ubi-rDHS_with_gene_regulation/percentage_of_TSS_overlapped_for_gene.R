
# -- Kaili
# This script is for plotting tissue specificity index and finding the cut-off.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/")

######################

data1 = read.table("hg38_gene_all_TSS_count_merged.txt", header = FALSE, row.names = 2)
data2 = read.table("hg38_gene_overlapped_TSS_count_merged.txt", header = FALSE, row.names = 2)

matrix = data.frame(cbind(data1[,1], data2[rownames(data1),1],data2[rownames(data1),1]/data1[,1]))
rownames(matrix) = rownames(data1[,1])
colnames(matrix) = c("total", "overlapped", "percentage")

pdf("percentage_of_TSS_overlapped.pdf", width = 14)
boxplot(matrix[matrix$total==1,]$percentage, matrix[matrix$total==2,]$percentage, 
        matrix[matrix$total==3,]$percentage, matrix[matrix$total==4,]$percentage,
        matrix[matrix$total==5,]$percentage, matrix[matrix$total==6,]$percentage, 
        matrix[matrix$total==7,]$percentage, matrix[matrix$total==8,]$percentage, 
        matrix[matrix$total==9,]$percentage, matrix[matrix$total==10,]$percentage, 
        matrix[matrix$total>10 & matrix$total<21,]$percentage,
        matrix[matrix$total>20 & matrix$total<31,]$percentage,
        matrix[matrix$total>30 & matrix$total<41,]$percentage,
        matrix[matrix$total>40,]$percentage, 
        names = c("1","2","3","4","5","6","7","8","9","10","11-20","21-30","31-40","41-52"),
        xlab = "total number of TSSs of gene", ylab = "percentage of overlappted TSSs", col = "grey")
dev.off()

pdf("number_of_TSS_overlapped.pdf", width = 14)
boxplot(matrix[matrix$total==1,]$overlapped, matrix[matrix$total==2,]$overlapped, 
        matrix[matrix$total==3,]$overlapped, matrix[matrix$total==4,]$overlapped,
        matrix[matrix$total==5,]$overlapped, matrix[matrix$total==6,]$overlapped, 
        matrix[matrix$total==7,]$overlapped, matrix[matrix$total==8,]$overlapped, 
        matrix[matrix$total==9,]$overlapped, matrix[matrix$total==10,]$overlapped, 
        matrix[matrix$total>10 & matrix$total<21,]$overlapped,
        matrix[matrix$total>20 & matrix$total<31,]$overlapped,
        matrix[matrix$total>30 & matrix$total<41,]$overlapped,
        matrix[matrix$total>40,]$overlapped, 
        names = c("1","2","3","4","5","6","7","8","9","10","11-20","21-30","31-40","41-52"),
        xlab = "total number of TSSs of gene", ylab = "number of overlappted TSSs", col = "grey")
dev.off()

