
# -- Kaili
# This script is for making figs by PANTHER result.
# EXP: Rscript make_figs_from_PANTHER.R "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/"
# "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/GO_BP_fisher_clean.txt"
# "hg19_ubi-rDHS_BP_15" "10,921 ubi-rDHS (top 15)" 15
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE);
workDir = args[1]
inputFile = args[2]
outFile = args[3]
main = args[4]
n = args[5]

# workDir = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/"
# inputFile = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/GO/hg19_ubi-rDHS_BP_clean.txt"
# outFile = "hg19_ubi-rDHS_BP_15.pdf"
# main = "10,921 ubi-rDHS (top 15)"
# n=15

setwd(workDir)

data = read.csv(inputFile, header = TRUE, sep = "\t")
a=colnames(data)[1]
data_order = data[order(data[,8], decreasing = FALSE),]
dat=data_order[1:n,]
type=paste(unlist(strsplit(a, "[.]"))[1],unlist(strsplit(a, "[.]"))[2], 
           unlist(strsplit(a, "[.]"))[3], sep=" ")
colnames(dat) = c("id", "V2", "V3","V4","V5","fe","p","fdr")
#
go = c()  
for(i in 1:n){
  a = unlist(strsplit(as.character(dat[i,1]), " [(]"))[1] 
  go = c(go, a)
}
#
dat = cbind(dat, go)
#
ggplot(dat, aes(x=go, y=fe)) + geom_bar(stat = "identity") + coord_flip() + 
  labs(title=main, y="fold enrichment", x=type) + theme(axis.text.y = element_text(size=12)) + 
  scale_fill_manual(values="#41ae76") + 
  geom_text(aes(label = format(fdr, scientific = TRUE)), colour = "white", position = position_stack(vjust = 0.5), 
            size = 5, fontface = "bold") +
  scale_x_discrete(limits = go[n:1])
ggsave(paste(outFile,".pdf",sep=""), width = 9, height = 9)
ggsave(paste(outFile,".png",sep=""), width = 9, height = 9)


