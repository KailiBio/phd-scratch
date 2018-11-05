
# -- Kaili
# This script is for making two-way boxplot for RNA-seq data between ubi-rOCRs overlapped vs. rOCR overlapped genes.

# workDir = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/gene_exp_comparison_file/"
# filePath = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/gene_exp_comparison_file/"
# ubi_file = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/GRCh38_ubi-rOCR_overlapped_gene_id.txt"
# title = "RNA-seq between\nubi-rOCRs vs. active-rOCRs overlapped genes"
# xlab = "ubi-rOCRs overlapped genes\nlog10(TPM+0.1)"
# ylab = "rest cell-type active-rOCRs overlapped genes\nlog10(TPM+0.1)"
# outFile = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/gene_exp_comparison_pdf/hg38_RNA_ubi-rOCR_rOCR_two-way_boxplot.pdf"

# workDir = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/TSS_exp_comparison_file/"
# filePath = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/TSS_exp_comparison_file/"
# ubi_file = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/GRCh38_ubi-rOCR_overlapped_TSS_uniqID_list.txt"
# title = "RAMPAGE\nubi-rOCRs vs. active-rOCRs overlapped TSSs"
# xlab = "ubi-rOCRs overlapped TSSs\nlog10(TPM+0.1)"
# ylab = "rest cell-type active-rOCRs overlapped TSSs\nlog10(TPM+0.1)"
# outFile = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/TSS_exp_comparison_pdf/hg38_TSS_ubi-rOCR_rOCR_two-way_boxplot.pdf"

# workDir = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/mergedTSS_exp_comparison_file/"
# filePath = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/mergedTSS_exp_comparison_file/"
# ubi_file = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/GRCh38_ubi-rOCR_overlapped_merged-TSS_list.txt"
# title = "RAMPAGE\nubi-rOCRs vs. active-rOCRs overlapped merged-TSSs"
# xlab = "ubi-rOCRs overlapped merged-TSSs\nlog10(TPM+0.1)"
# ylab = "rest cell-type active-rOCRs overlapped merged-TSSs\nlog10(TPM+0.1)"
# outFile = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/mergedTSS_exp_comparison_pdf/hg38_mergedTSS_ubi-rOCR_rOCR_two-way_boxplot.pdf"

# workDir = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/dnase_comparison_file/"
# filePath = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/dnase_comparison_file/"
# ubi_file = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/GRCh38_ubi-rOCRs_list.txt"
# title = "DNase-seq\nubi-rOCRs vs. active-rOCRs"
# xlab = "ubi-rOCRs\nz-score"
# ylab = "rest cell-type active-rOCRs\nz-score"
# outFile = "/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/dnase_comparison_pdf/hg38_DNase_ubi-rOCR_rOCR_two-way_boxplot.pdf"


args <- commandArgs(TRUE)
workDir = args[1]
filePath = args[2]
ubi_file = args[3]
title = args[4]
xlab = args[5]
ylab = args[6]
outFile = args[7]


setwd(workDir)
ubi = as.vector(read.table(ubi_file)[,1])

# get all median and sd
l = list.files(path=filePath)
for(i in 1:length(l)){
  data = read.table(paste(filePath, l[i], sep=""), header = F, row.names = 1)
  ubi_overlapped = log10(data[ubi,1]+0.1)
  active_overlapped = log10(data[! rownames(data)%in% ubi,1]+0.1)
  # ubi_overlapped = data[ubi,1]
  # active_overlapped = data[! rownames(data)%in% ubi,1]
  if(i==1){
    dat = c(median(ubi_overlapped), sd(ubi_overlapped), median(active_overlapped), sd(active_overlapped))
  }else{
    dat = rbind(dat, c(median(ubi_overlapped), sd(ubi_overlapped), median(active_overlapped), sd(active_overlapped)))
  }
}
dat = data.frame(dat) 
colnames(dat) = c("ubi_m", "ubi_sd", "res_m", "res_sd")

# make two-way boxplot
library(ggplot2)

ggplot(dat, aes(ubi_m, res_m)) + geom_point(col = "#525252") + 
  geom_errorbarh(aes(xmax = ubi_m + ubi_sd, xmin = ubi_m - ubi_sd), height = 0.1, col = "#ef8a62") +
  geom_errorbar(aes(ymax = res_m + res_sd, ymin = res_m - res_sd), width=0.1, col = "#67a9cf") +
  ylim(-1.5,0.7) + xlim(-1.5,0.7) + theme_minimal() + labs(title = title) + xlab(xlab) + ylab(ylab)
ggsave(outFile)


