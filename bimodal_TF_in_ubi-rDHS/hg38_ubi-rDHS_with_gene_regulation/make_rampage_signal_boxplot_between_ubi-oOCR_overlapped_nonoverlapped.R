# -- Kaili
# This script is for making scatter plot of RAMPAGE signal between tissues.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/rampage_tissue/")
library(ggplot2)

##############
files = list.files()
d = read.table("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/hg38_ubi-rOCR_overlapped_gene_all_TSS_uniqID.txt")
tss = unique(as.vector(d[,2]))

pdf("TSS_RAMPAGE_between_tissue_boxplot.pdf", onefile = TRUE)
for (i in 1:length(files)){
  data = read.table(files[i], row.names = 1)
  tissue = unlist(strsplit(files[i], "_rampage"))[1]
  matrix = transform(data[tss,], type="no-overlapped")
  rownames(matrix) = tss
  colnames(matrix) = c("tissue", "type")
  matrix$type = as.vector(matrix$type)
  matrix[as.vector(d[d$V3=="overlapped",2]),]$type = "overlapped"
  
  ggplot(matrix, aes(y=log10(tissue+0.1), x=type, fill=type)) +geom_boxplot() + 
    labs(title=paste("RAMPAGE signal of TSSs in ", tissue, sep="")) + ylab(tissue)
}
dev.off()

