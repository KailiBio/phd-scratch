
# -- Kaili
# This script is for making aggregation plot for each sample.

args = commandArgs(trailingOnly=TRUE)
sample = args[1]

# sample = "embryonic-facial-prominence_11.5"
# workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/"

############
setwd("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/")
library(ggplot2)
library(grid)
library(gridExtra)

assay_list = c("ATAC", "DNAme", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3", "H3K9ac", "H3K27me3", "H3K27ac", "H3K36me3")

# setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/ss/")
# sample = "embryonic-facial-prominence_11.5"
##################
# ATAC
##################
assay=assay_list[1]
dhs_ocr = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/", 
                sample, "_", assay, "_OCR.txt", sep="")
normal = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal/", 
               sample, "_", assay, "_normal.txt", sep="")
#
dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, normal_data)
colnames(matrix) = c("signal", "type", "loci")
#
p1 = ggplot(matrix, aes(x=loci, y=signal, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = assay) + 
  annotate("text", x=25, y=max(matrix$signal)+0.1, label='bold("upstream")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$signal)+0.1, label='bold("bin")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$signal)+0.1, label='bold("downstream")', color = "black", size=5, parse = TRUE) +
  theme(axis.text.x = element_blank(), legend.position="none")
##################
# DNAme
##################
assay=assay_list[2]
dhs_ocr = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/", 
                sample, "_", assay, "_OCR.txt", sep="")
normal = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal/", 
               sample, "_", assay, "_normal.txt", sep="")
#
dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, normal_data)
colnames(matrix) = c("percentage", "type", "loci")
#
p2 = ggplot(matrix, aes(x=loci, y=percentage, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = assay) + 
  annotate("text", x=25, y=max(matrix$percentage)+0.1, label='bold("upstream")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$percentage)+0.1, label='bold("bin")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$percentage)+0.1, label='bold("downstream")', color = "black", size=5, parse = TRUE) +
  theme(axis.text.x = element_blank(), legend.position="none")
##################
# H3K4me1
##################
assay=assay_list[3]
dhs_ocr = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/", 
                sample, "_", assay, "_OCR.txt", sep="")
normal = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal/", 
               sample, "_", assay, "_normal.txt", sep="")
#
dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, normal_data)
colnames(matrix) = c("signal", "type", "loci")
#
p3 = ggplot(matrix, aes(x=loci, y=signal, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = assay) + 
  annotate("text", x=25, y=max(matrix$signal)+0.1, label='bold("upstream")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$signal)+0.1, label='bold("bin")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$signal)+0.1, label='bold("downstream")', color = "black", size=5, parse = TRUE) +
  theme(axis.text.x = element_blank(), legend.position="none")
##################
# H3K4me2
##################
assay=assay_list[4]
dhs_ocr = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/", 
                sample, "_", assay, "_OCR.txt", sep="")
normal = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal/", 
               sample, "_", assay, "_normal.txt", sep="")
#
dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, normal_data)
colnames(matrix) = c("signal", "type", "loci")
#
p4 = ggplot(matrix, aes(x=loci, y=signal, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = assay) + 
  annotate("text", x=25, y=max(matrix$signal)+0.1, label='bold("upstream")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$signal)+0.1, label='bold("bin")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$signal)+0.1, label='bold("downstream")', color = "black", size=5, parse = TRUE) +
  theme(axis.text.x = element_blank(), legend.position="none")
##################
# H3K4me3
##################
assay=assay_list[5]
dhs_ocr = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/", 
                sample, "_", assay, "_OCR.txt", sep="")
normal = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal/", 
               sample, "_", assay, "_normal.txt", sep="")
#
dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, normal_data)
colnames(matrix) = c("signal", "type", "loci")
#
p5 = ggplot(matrix, aes(x=loci, y=signal, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = assay) + 
  annotate("text", x=25, y=max(matrix$signal)+0.1, label='bold("upstream")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$signal)+0.1, label='bold("bin")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$signal)+0.1, label='bold("downstream")', color = "black", size=5, parse = TRUE) +
  theme(axis.text.x = element_blank(), legend.position="none")
##################
# H3K9me3
##################
assay=assay_list[6]
dhs_ocr = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/", 
                sample, "_", assay, "_OCR.txt", sep="")
normal = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal/", 
               sample, "_", assay, "_normal.txt", sep="")
#
dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, normal_data)
colnames(matrix) = c("signal", "type", "loci")
#
p6 = ggplot(matrix, aes(x=loci, y=signal, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = assay) + 
  annotate("text", x=25, y=max(matrix$signal)+0.1, label='bold("upstream")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$signal)+0.1, label='bold("bin")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$signal)+0.1, label='bold("downstream")', color = "black", size=5, parse = TRUE) +
  theme(axis.text.x = element_blank(), legend.position="none")
##################
# H3K9ac
##################
assay=assay_list[7]
dhs_ocr = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/", 
                sample, "_", assay, "_OCR.txt", sep="")
normal = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal/", 
               sample, "_", assay, "_normal.txt", sep="")
#
dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, normal_data)
colnames(matrix) = c("signal", "type", "loci")
#
p7 = ggplot(matrix, aes(x=loci, y=signal, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = assay) + 
  annotate("text", x=25, y=max(matrix$signal)+0.1, label='bold("upstream")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$signal)+0.1, label='bold("bin")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$signal)+0.1, label='bold("downstream")', color = "black", size=5, parse = TRUE) +
  theme(axis.text.x = element_blank(), legend.position="none")
##################
# H3K27me3
##################
assay=assay_list[8]
dhs_ocr = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/", 
                sample, "_", assay, "_OCR.txt", sep="")
normal = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal/", 
               sample, "_", assay, "_normal.txt", sep="")
#
dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, normal_data)
colnames(matrix) = c("signal", "type", "loci")
#
p8 = ggplot(matrix, aes(x=loci, y=signal, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = assay) + 
  annotate("text", x=25, y=max(matrix$signal)+0.1, label='bold("upstream")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$signal)+0.1, label='bold("bin")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$signal)+0.1, label='bold("downstream")', color = "black", size=5, parse = TRUE) +
  theme(axis.text.x = element_blank(), legend.position="none")
##################
# H3K27ac
##################
assay=assay_list[9]
dhs_ocr = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/", 
                sample, "_", assay, "_OCR.txt", sep="")
normal = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal/", 
               sample, "_", assay, "_normal.txt", sep="")
#
dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, normal_data)
colnames(matrix) = c("signal", "type", "loci")
#
p9 = ggplot(matrix, aes(x=loci, y=signal, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = assay) + 
  annotate("text", x=25, y=max(matrix$signal)+0.1, label='bold("upstream")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$signal)+0.1, label='bold("bin")', color = "black", size=5, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$signal)+0.1, label='bold("downstream")', color = "black", size=5, parse = TRUE) +
  theme(axis.text.x = element_blank(), legend.position="none")
##################
# H3K36me3
##################
assay=assay_list[10]
dhs_ocr = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR/", 
                sample, "_", assay, "_OCR.txt", sep="")
normal = paste("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal/", 
               sample, "_", assay, "_normal.txt", sep="")
#
dhs_ocr_data = read.table(dhs_ocr)
dhs_ocr_data = transform(dhs_ocr_data, type="dhs_OCR", loci = 1:150)
normal_data = read.table(normal)
normal_data = transform(normal_data, type="normal_bins", loci = 1:150)
matrix = rbind(dhs_ocr_data, normal_data)
colnames(matrix) = c("signal", "type", "loci")
#
p10 = ggplot(matrix, aes(x=loci, y=signal, group=type, col = type)) + geom_line() + geom_point() +
  labs(title = assay) + 
  annotate("text", x=25, y=max(matrix$signal)+0.1, label='bold("upstream")', color = "black", size=3, parse = TRUE) +
  annotate("text", x=75, y=max(matrix$signal)+0.1, label='bold("bin")', color = "black", size=3, parse = TRUE) +
  annotate("text", x=125, y=max(matrix$signal)+0.1, label='bold("downstream")', color = "black", size=3, parse = TRUE) +
  theme(axis.text.x = element_blank())

##################
outfile = paste(sample, "_matched_aggregation.pdf", sep="")
lay = rbind(c(1,1,2,2,3,3),
            c(4,4,5,5,6,6),
            c(7,7,8,8,9,9),
            c(10,10,10,11,11,11))
p = grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol=3,
                 textGrob(sample,gp=gpar(fontsize=15,font=2)),
                 layout_matrix = lay)
ggsave(outfile,p, width = 12, height = 15)
##################
