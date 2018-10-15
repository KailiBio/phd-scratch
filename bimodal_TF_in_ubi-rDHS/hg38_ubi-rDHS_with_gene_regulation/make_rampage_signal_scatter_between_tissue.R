
# -- Kaili
# This script is for making scatter plot of RAMPAGE signal between tissues.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/rampage_tissue/")
library(ggplot2)

##############
files = list.files()
d = read.table("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/hg38_ubi-rOCR_overlapped_gene_all_TSS_uniqID.txt")
tss = unique(as.vector(d[,2]))

pdf("TSS_RAMPAGE_between_tissue_scatter.pdf", onefile = TRUE)
data1 = read.table(files[1], row.names = 1)
tissue1 = unlist(strsplit(files[1], "_rampage"))[1]
for(i in 2:length(files)){
  data2 = read.table(files[i], row.names = 1)
  tissue2 = unlist(strsplit(files[i], "_rampage"))[1]
  matrix = cbind(data1[tss,], data2[tss,])
  rownames(matrix) = tss
  # dat = matrix[matrix[,1]>1 | matrix[,2]>1,]
  # name = intersect(d, row.names(dat))
  # dat = transform(dat, type="no-overlapped")
  # colnames(dat) = c("tissue1", "tissue2", "type")
  # dat$type = as.vector(dat$type)
  # dat[d,]$type = "overlapped"
  matrix = transform(matrix, type="no-overlapped")
  colnames(matrix) = c("tissue1", "tissue2","type")
  matrix$type = as.vector(matrix$type)
  matrix[as.vector(d[d$V3=="overlapped",2]),]$type = "overlapped"
  
  ggplot(matrix, aes(x=log10(tissue1+0.1), y=log10(tissue2+0.1), col=type)) +geom_count() +
    scale_size_area() + labs(title="RAMPAGE signal of TSSs") + xlab(tissue1) + ylab(tissue2)
}
dev.off()

                                          