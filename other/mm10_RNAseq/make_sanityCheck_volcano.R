
# -- Kaili
# This script is for making volcano plot for sanity check DE results.

setwd("/data/zusers/fankaili/ccre/mm10_rnaseq/")
#setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/mm10_RNAseq/")


library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
file = args[1]
sample = args[2]

# file="tmp.matrix5.txt"

########
dat = read.table(file, header = TRUE, row.names = 1)
num = table(dat$sig)

p1 = ggplot(dat, aes(x=log2FoldChange, y=-log10(padj), col = sig)) + geom_point(alpha=0.3) +
  theme_minimal() + theme(title = element_text(size=12, face="bold")) +
  labs(title=paste(sample,"\n","Junko's DE call",sep="")) + coord_cartesian(xlim=c(-20,20), ylim=c(0,180))


p2 = ggplot(dat, aes(x=my_log2FoldChange, y=-log10(my_padj), col=sig)) + geom_point(alpha=0.3) +
  theme_minimal() + theme(title = element_text(size=12, face="bold")) +
  labs(paste(sample,"\n",title="My DE call",sep="")) + coord_cartesian(xlim=c(-20,20), ylim=c(0,180))

p = grid.arrange(p1, p2, ncol=2, name="ss")
p
ggsave(paste("./volcano_check/",sample,".pdf", sep=""),p, width=12, height = 7)
ggsave(paste("./volcano_check/",sample,".png", sep=""),p, width=12, height = 7)


