
# -- Kaili
# This script is for making scatter plot for qvalue.

setwd("/data/zusers/fankaili/ccre/mm10_rnaseq/")
#setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/mm10_RNAseq/")


library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
sample = args[1]

# sample="midbrain_0_VS_facial_13.5"

########
input_file = paste("./compare_DEGs_M4_M18/",sample,"_matrix.txt",sep="")
dat = read.table(input_file, header = TRUE, row.names = 1)
num = table(dat$sig)

fc_type = cor(dat$log2FoldChange[1:10],dat$my_log2FoldChange[1:10])
fc_type = 1
if(fc_type <0){
  dat$log2FoldChange = - dat$log2FoldChange
}
p1 = ggplot(dat, aes(x=log2FoldChange, y=my_log2FoldChange, col=sig)) + geom_point(alpha=0.5) +
  theme_minimal() + theme(legend.position = "none") +
  coord_cartesian(xlim=c(-20,20), ylim=c(-20,20)) +
  geom_abline(slope = 1, col="blue", linetype="dashed") +
  labs(title = "log2FC")
  
p2 = ggplot(dat, aes(x=-log10(padj), y=-log10(my_padj), col=sig)) + geom_point(alpha=0.5) +
  theme_minimal() + theme(legend.position = c(0.8,0.2)) +
  coord_cartesian(xlim=c(0,300), ylim=c(0,300)) +
  geom_abline(slope = 1, col="blue", linetype="dashed") +
  labs(title = "qvalue")

p3 = ggplot(dat, aes(x=log2FoldChange, y=-log10(padj), col = sig)) + geom_point(alpha=0.5) +
  theme_minimal() + theme(legend.position = "none") +
  labs(title="M4's DE call") + coord_cartesian(xlim=c(-20,20), ylim=c(0,300))


p4 = ggplot(dat, aes(x=my_log2FoldChange, y=-log10(my_padj), col=sig)) + geom_point(alpha=0.5) +
  theme_minimal() + theme(legend.position = "none") +
  labs(title="M18's DE call") + coord_cartesian(xlim=c(-20,20), ylim=c(0,300))


p = grid.arrange(p1, p2, p3, p4, ncol=2, top = sample)
p

#out_file1 = paste("./compare_DEGs_M4_M18/", sample, "_readCounts.pdf", sep="")
out_file2 = paste("./compare_DEGs_M4_M18/", sample, "_readCounts.png", sep="")
#ggsave(out_file1, p, width = 11, height = 11)
ggsave(out_file2, p, width = 11, height = 11)
