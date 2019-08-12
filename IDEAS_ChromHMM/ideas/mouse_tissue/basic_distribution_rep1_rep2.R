
# -- Kaili
# This script is for making basic violin and histogram for rep1 & rep2.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/rep2/")

library(ggplot2)
library(data.table)
library(gridExtra)

mark_list = c("ATAC", "DNAme", "CTCF", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3", "H3K9ac", "H3K27me3", 
              "H3K27ac", "H3K36me3")

###############
# violin plot
###############
sample="lung_14.5"
sample="lung_0"

pdf(paste("reps_violin_",sample,".pdf",sep=""), width=12)
for(i in 1:11){
  mark = mark_list[i]
  dat_rep1 = fread(paste("rep1_signal/",sample,"_",mark,"_dhs.txt", sep=""))
  dat_rep2 = fread(paste("rep2_signal/",sample,"_",mark,"_dhs.txt", sep=""))
  #
  dat = rbind(transform(dat_rep1, rep="rep1"), transform(dat_rep2, rep="rep2"))
  colnames(dat) = c("signal","rep")
  p1 = ggplot(dat, aes(x= rep, y=signal)) + geom_violin() + 
    theme_minimal() + theme(text = element_text(face="bold", size=12)) +
    labs(title=mark) + xlab("")
  p2 = ggplot(dat, aes(x= rep, y=signal)) + geom_violin() + 
    theme_minimal() + theme(text = element_text(face="bold", size=12)) +
    labs(title=mark) + xlab("") +
    coord_cartesian(ylim=c(0,200))
  p3 = ggplot(dat, aes(x= rep, y=log10(signal+0.01))) + geom_violin() + 
    theme_minimal() + theme(text = element_text(face="bold", size=12)) +
    labs(title=mark) + xlab("")
  p = grid.arrange(p1, p2, p3, ncol=3)
  print(p)
}
dev.off()

###############
# histogram
###############
sample="lung_14.5"
sample="lung_0"

pdf(paste("reps_histogram_",sample,".pdf",sep=""), width = 12)
for(i in 1:11){
  mark = mark_list[i]
  dat_rep1 = fread(paste("rep1_signal/",sample,"_",mark,"_dhs.txt", sep=""))
  dat_rep2 = fread(paste("rep2_signal/",sample,"_",mark,"_dhs.txt", sep=""))
  colnames(dat_rep1) = colnames(dat_rep2) = c("signal")
  #
  p1 = ggplot(dat_rep1, aes(x=log10(signal+0.01))) + geom_histogram(binwidth = 0.01, aes(y=..density..)) +
    theme_minimal() + theme(text = element_text(face="bold", size=12)) +
    labs(title=paste(mark," rep1", sep=""))
  p2 = ggplot(dat_rep2, aes(x=log10(signal+0.01))) + geom_histogram(binwidth = 0.01, aes(y=..density..)) +
           theme_minimal() + theme(text = element_text(face="bold", size=12)) +
           labs(title=paste(mark," rep2", sep=""))
  p = grid.arrange(p1, p2, ncol=2)
  print(p)
}
dev.off()


###############