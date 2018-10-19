
# -- Kaili
# This scirpt is for making boxplot for RAMPAGE signal between overlapped merged-TSSs and rest.
# EXP: Rscript make_rampage_signal_boxplot_overlapped_merged-TSS.R
#       "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/" "A172"
#       "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_between_tissue/"

args<-commandArgs(TRUE)
file_dir = args[1]
tissue = args[2]
out_dir = args[3]

# file_dir = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/"
# tissue = "sigmoid_colon_51_year"
# out_dir = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/"

#####################
setwd(file_dir)
library(ggplot2)

data = read.table(paste(tissue, "_rampage.txt", sep=""), row.names = 1)
overlapped = as.vector(read.table("/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/GRCh38_ubi-rOCR_overlapped_merged-TSS.bed")[,4])
matrix = transform(data, type = "no_overlapped")
colnames(matrix) = c("signal", "type")
matrix$type = as.vector(matrix$type)
matrix[overlapped,]$type = "overlapped"

pdf(paste(out_dir,"RAMPAGE_merged-TSS_",tissue,"_boxplot.pdf",sep=""))
ggplot(matrix, aes(y=log10(signal+0.1), x=type, fill=type)) +geom_boxplot() + 
  labs(title=paste("RAMPAGE signal of merged-TSSs in ", tissue, sep="")) + ylab(tissue)
dev.off()

#####################