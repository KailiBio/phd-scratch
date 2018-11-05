
# -- Kaili
# This script is for making scatter plot for given tissue RAMPAGE signal.
# EXP: Rscript make_rampage_signal_scatter_between_tissue_separate.R
#      "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_mergedTSS/" "A172" "K562"
#      "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue_mergedTSS/"

args<-commandArgs(TRUE)
file_dir = args[1]
tissue1 = args[2]
tissue2 = args[3]
out_dir = args[4]

# file_dir = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/"
# tissue1 = "sigmoid_colon_51_year"
# tissue2 = "stomach_51_year"
# out_dir = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/merged-TSS/"

#####################
setwd(file_dir)
library(ggplot2)

data1 = read.table(paste(tissue1, "_rampage.txt", sep=""), row.names = 1)
data2 = read.table(paste(tissue2, "_rampage.txt", sep=""), row.names = 1)
overlapped = as.vector(read.table("/data/zusers/fankaili/ccre/hg38_ubi-rDHS/merged-TSS/GRCh38_ubi-rOCR_overlapped_merged-TSS_list.txt")[,1])

matrix = cbind(data1, data2)
matrix = transform(matrix, type="no-overlapped")
colnames(matrix) = c("tissue1", "tissue2","type")
matrix$type = as.vector(matrix$type)
matrix[overlapped,]$type = "overlapped"

pdf(paste(out_dir, "RAMPAGE_merged-TSS_", tissue1,"_", tissue2, ".pdf", sep=""))
ggplot(matrix, aes(x=log10(tissue1+0.1), y=log10(tissue2+0.1), col=type)) +geom_count() +
  scale_size_area() + labs(title="RAMPAGE signal of merged-TSSs") + xlab(tissue1) + ylab(tissue2)
dev.off()
