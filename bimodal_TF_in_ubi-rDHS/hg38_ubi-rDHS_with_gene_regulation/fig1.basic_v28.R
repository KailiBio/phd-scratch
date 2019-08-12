
# -- Kaili
# This script is for making boxplot figues for comparing DNase&RNA&RAMPAGE signal.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/fig1/basic_v28/")
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
  geom_boxplot(width=0.6, outlier.shape = NA) +
  theme_classic() + ylab("DNase signal") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        legend.position = "none") + coord_cartesian(ylim=c(1,6)) +
  scale_fill_manual(values=c("#E73A2F", "#397AF2"))
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
  geom_boxplot(width=0.6, outlier.shape = NA) + 
  theme_classic() + ylab("Gene expression\nlog10(TPM+0.1)") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        legend.position = "none") + coord_cartesian(ylim=c(-1,5)) +
  scale_fill_manual(values=c("#E73A2F", "#397AF2"))
p2


# statistic
for(i in 1:length(cell)){
  a = as.numeric(RNA[RNA$cell_line==cell[i] & RNA$group=="a",]$signal)
  b = as.numeric(RNA[RNA$cell_line==cell[i] & RNA$group=="b",]$signal)
  p = wilcox.test(a,b)$p.value
  print(p)
}
# 5.931476e-128
###################
# RAMPAGE
###################
RAMPAGE = read.table("sample_RAMPAGE_signal_comparison.txt")
colnames(RAMPAGE) = c("OCR", "signal", "ubi", "cell_line", "sum")
RAMPAGE=transform(RAMPAGE, group="a")
RAMPAGE$group = as.vector(RAMPAGE$group)
RAMPAGE[RAMPAGE$ubi=="TSS_overlap_non-ubi_active-rOCRs",]$group="b"
p3 = ggplot(data=RAMPAGE, aes(x=cell_line, y=log10(signal+0.1), fill=group)) + 
  geom_boxplot(width=0.6, outlier.shape = NA) + 
  theme_classic() + ylab("TSS activity\nlog10(TPM+0.1)") +
  theme(axis.text.x = element_text(size=12, face="bold", angle=45, vjust=0.6, hjust=0.6), 
        axis.text.y = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x=element_blank(), legend.position = "none") +
  coord_cartesian(ylim=c(-1,3)) + scale_fill_manual(values=c("#E73A2F", "#397AF2"))
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
ggsave("fig1.basic_v28_newcolor.pdf", p)

###################
# lincRNA comparison
###################
setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/fig1/basic_v28/lincRNA/")
library(ggplot2)
library(grid)
library(gridExtra)

cell_list = c("A172","Daoy","GM23248","GM23338","hepatocyte","HT1080","LHCN-M2","myotube",
         "NCI-H460","neural_progenitor_cell","RPMI-7951","SJCRH30","SJSA1",
         "skeletal_muscle_myoblast","SK-MEL-5","SK-N-DZ")

###########
# gene
###########
RNA = read.table("sample_RNA_signal_comparison_PC_linc_new.txt")
colnames(RNA) = c("OCR", "signal", "ubi", "cell_line", "sum","PC","type")

RNA$type = as.character(RNA$type)
RNA[RNA$type=="other-PC",]$type = "PC_overlapping_non-ubi-rOCR"
RNA[RNA$type=="ubi-PC",]$type = "PC_overlapping_ubi-rOCR"
RNA[RNA$type=="other-lincRNA",]$type = "lincRNA_overlapping_non-ubi-rOCR"
RNA[RNA$type=="ubi-lincRNA",]$type = "lincRNA_overlapping_ubi-rOCR"

ggplot(data=RNA, aes(x=cell_line, y=log10(signal+0.1), fill=type)) + 
  geom_boxplot(width=0.5, outlier.shape = NA) + 
  theme_classic() + ylab("Gene expression\nlog10(TPM+0.1)") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12)) +
  scale_fill_manual(values = c("#00B0F0","#397AF2","#ffaaaa","#E73A2F")) +
  coord_cartesian(ylim=c(-1,3.2))

ggsave("PC_lincRNA_RNA.pdf", width = 12, height=4)
ggsave("PC_lincRNA_RNA.png", width = 12, height=4)

RNA2 = RNA[RNA$PC=="lincRNA",]
ggplot(data=RNA2, aes(x=cell_line, y=log10(signal+0.1), fill=type)) + 
  geom_boxplot(width=0.5, outlier.shape = NA) + 
  theme_classic() + ylab("Gene expression\nlog10(TPM+0.1)") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12)) +
  scale_fill_manual(values = c("#00B0F0","#ffaaaa")) +
  coord_cartesian(ylim=c(-1,2))

ggsave("lincRNA_RNA.pdf", width = 12, height=4)
ggsave("lincRNA_RNA.png", width = 12, height=4)

for(i in 1:length(cell_list)){
  print(cell[i])
  x = RNA2[RNA2$type=="ubi-lincRNA" & RNA2$cell_line==cell_list[i],]$signal
  y = RNA2[RNA2$type=="other-lincRNA" & RNA2$cell_line==cell_list[i],]$signal
  print(wilcox.test(x,y)$p.value)
}

###########
# TSS
###########

RAMPAGE = read.table("sample_RAMPAGE_signal_comparison_PC_linc_new.txt")
colnames(RAMPAGE) = c("OCR", "signal", "ubi", "cell_line", "sum","PC","type")

RAMPAGE$type = as.character(RAMPAGE$type)
RAMPAGE[RAMPAGE$type=="other-PC",]$type = "PC_overlapping_non-ubi-rOCR"
RAMPAGE[RAMPAGE$type=="ubi-PC",]$type = "PC_overlapping_ubi-rOCR"
RAMPAGE[RAMPAGE$type=="other-lincRNA",]$type = "lincRNA_overlapping_non-ubi-rOCR"
RAMPAGE[RAMPAGE$type=="ubi-lincRNA",]$type = "lincRNA_overlapping_ubi-rOCR"

ggplot(data=RAMPAGE, aes(x=cell_line, y=log10(signal+0.1), fill=type)) + 
  geom_boxplot(width=0.5, outlier.shape = NA) + 
  theme_classic() + ylab("TSS activity\nlog10(TPM+0.1)") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12)) +
  scale_fill_manual(values = c("#00B0F0","#397AF2","#ffaaaa","#E73A2F")) +
  coord_cartesian(ylim=c(-1,3.2))

ggsave("PC_lincRNA_RAMPAGE.pdf", width = 12, height=4)
ggsave("PC_lincRNA_RAMPAGE.png", width = 12, height=4)

RAMPAGE2 = RAMPAGE[RAMPAGE$PC=="lincRNA",]
ggplot(data=RAMPAGE2, aes(x=cell_line, y=log10(signal+0.1), fill=type)) + 
  geom_boxplot(width=0.5, outlier.shape = NA) + 
  theme_classic() + ylab("TSS activity\nlog10(TPM+0.1)") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12)) +
  scale_fill_manual(values = c("#00B0F0","#ffaaaa")) +
  coord_cartesian(ylim=c(-1,1.3))

ggsave("lincRNA_RAMPAGE.pdf", width = 12, height=4)
ggsave("lincRNA_RAMPAGE.png", width = 12, height=4)

for(i in 1:length(cell_list)){
  print(cell[i])
  x = RAMPAGE2[RAMPAGE2$type=="ubi-lincRNA" & RAMPAGE2$cell_line==cell_list[i],]$signal
  y = RAMPAGE2[RAMPAGE2$type=="other-lincRNA" & RAMPAGE2$cell_line==cell_list[i],]$signal
  print(wilcox.test(x,y)$p.value)
}

RAMPAGE3 = RAMPAGE[RAMPAGE$type=="ubi-lincRNA" | RAMPAGE$type=="other-PC",]
ggplot(data=RAMPAGE3, aes(x=cell_line, y=log10(signal+0.1), fill=type)) + 
  geom_boxplot(width=0.5, outlier.shape = NA) + 
  theme_minimal() + ylab("Expression of TSSs\nlog10(TPM+0.1)") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12)) +
  scale_fill_manual(values = c("#00B0F0","#ffaaaa")) +
  coord_cartesian(ylim=c(-1,1.3))

ggsave("ubi-lincRNA_other-PC_RAMPAGE.pdf", width = 12)
ggsave("ubi-lincRNA_other-PC_RAMPAGE.png", width = 12)

for(i in 1:length(cell_list)){
  print(cell[i])
  x = RAMPAGE3[RAMPAGE3$type=="ubi-lincRNA" & RAMPAGE3$cell_line==cell_list[i],]$signal
  y = RAMPAGE3[RAMPAGE3$type=="other-PC" & RAMPAGE3$cell_line==cell_list[i],]$signal
  print(wilcox.test(x,y)$p.value)
}
###################

