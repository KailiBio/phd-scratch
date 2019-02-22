
# -- Kaili
# This script is for making boxplot figues for comparing DNase&RNA&RAMPAGE signal.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/fig1/comprehensive_v28/")
library(ggplot2)
library(grid)
library(gridExtra)

cell = c("A172","Daoy","GM23248","GM23338","hepatocyte","HT1080","LHCN-M2","myotube",
         "NCI-H460","neural_progenitor_cell","RPMI-7951","SJCRH30","SJSA1",
         "skeletal_muscle_myoblast","SK-MEL-5","SK-N-DZ")
###################
# DNase
###################
DNase = read.table("sample_DNase_signal_comparison.txt")
colnames(DNase) = c("OCR", "signal", "ubi", "cell_line", "sum")
DNase=transform(DNase, group="a")
DNase$group = as.vector(DNase$group)
DNase[DNase$ubi=="non-ubi_active-rOCR",]$group="b"
p1 = ggplot(data=DNase, aes(x=cell_line, y=signal, fill=group)) + 
  geom_boxplot(width=0.3, outlier.shape = NA) + 
  theme_classic() + ylab("signal of DNase I") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        legend.position = "none") + coord_cartesian(ylim=c(1,6))
p1

# statistic
for(i in 1:length(cell)){
  a = as.numeric(DNase[DNase$cell_line==cell[i] & DNase$ubi=="ubi-rOCR",]$signal)
  b = as.numeric(DNase[DNase$cell_line==cell[i] & DNase$ubi=="non-ubi_active-rOCR",]$signal)
  p = wilcox.test(a,b)$p.value
  print(p)
}
###################
# RNA
###################
RNA = read.table("sample_RNA_signal_comparison.txt")
colnames(RNA) = c("OCR", "signal", "ubi", "cell_line", "sum")
RNA=transform(RNA, group="a")
RNA$group = as.vector(RNA$group)
RNA[RNA$ubi=="genes_whose_TSSs_overlap_other_active_rOCRs",]$group="b"
p2 = ggplot(data=RNA, aes(x=cell_line, y=log10(signal+0.1), fill=group)) + 
  geom_boxplot(width=0.3, outlier.shape = NA) + 
  theme_classic() + ylab("Expression of transcripts\nlog10(TPM+0.1)") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        legend.position = "none") + coord_cartesian(ylim=c(-1,5))
p2


# statistic
for(i in 1:length(cell)){
  a = as.numeric(RNA[RNA$cell_line==cell[i] & RNA$group=="a",]$signal)
  b = as.numeric(RNA[RNA$cell_line==cell[i] & RNA$group=="b",]$signal)
  p = wilcox.test(a,b)$p.value
  print(p)
}
# 2.158401e-148
###################
# RAMPAGE
###################
RAMPAGE = read.table("sample_RAMPAGE_signal_comparison.txt")
colnames(RAMPAGE) = c("OCR", "signal", "ubi", "cell_line", "sum")
RAMPAGE=transform(RAMPAGE, group="a")
RAMPAGE$group = as.vector(RAMPAGE$group)
RAMPAGE[RAMPAGE$ubi=="TSS_overlap_non-ubi_active-rOCRs",]$group="b"
p3 = ggplot(data=RAMPAGE, aes(x=cell_line, y=log10(signal+0.1), fill=group)) + 
  geom_boxplot(width=0.3, outlier.shape = NA) +
  theme_classic() + ylab("Expression of TSSs\nlog10(TPM+0.1)") +
  theme(axis.text.x = element_text(size=12, face="bold", angle=45, vjust=0.6, hjust=0.6), 
        axis.text.y = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x=element_blank(), legend.position = "none") +
  coord_cartesian(ylim=c(-1,3))
p3

# statistic
for(i in 1:length(cell)){
  a = as.numeric(RAMPAGE[RAMPAGE$cell_line==cell[i] & RAMPAGE$group=="a",]$signal)
  b = as.numeric(RAMPAGE[RAMPAGE$cell_line==cell[i] & RAMPAGE$group=="b",]$signal)
  p = wilcox.test(a,b)$p.value
  print(p)
}


###################
p = grid.arrange(p1, p2, p3, ncol=1,
                 layout_matrix = rbind(c(1),c(2),c(3),c(3)))
ggsave("fig1.comprehensive_v28.pdf", p)

###################

