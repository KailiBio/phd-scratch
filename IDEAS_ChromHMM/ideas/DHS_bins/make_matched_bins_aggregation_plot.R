
# -- Kaili
# This script is for making aggregation plot.

args = commandArgs(trailingOnly=TRUE)
prefix = args[1]
workDir = args[2]

# prefix = "embryonic-facial-prominence_14.5_ATAC"
# workDir = "/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/"

############
setwd(workDir)
library(ggplot2)

dhs_ocr = paste(prefix, "_dhs_OCR.txt", sep="")
dhs_gap = paste(prefix, "_flanking_gap.txt", sep="")
normal = paste(prefix, "_matched_normal.txt", sep="")
outfile = paste(prefix, "_matched_aggregation.pdf", sep="")

dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
dhs_gap_data = read.table(dhs_gap)
dhs_gap_data = transform(dhs_gap_data, type="dhs_gap", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, dhs_gap_data, normal_data)
colnames(matrix) = c("signal", "type", "loci")

ggplot(matrix, aes(x=loci, y=signal, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = prefix) + 
  annotate("text", x=25, y=max(matrix$signal)+0.1, label='bold("upstream")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$signal)+0.1, label='bold("bin")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$signal)+0.1, label='bold("downstream")', color = "black", size=5, parse = TRUE) +
  theme(axis.text.x = element_blank())
ggsave(outfile)

