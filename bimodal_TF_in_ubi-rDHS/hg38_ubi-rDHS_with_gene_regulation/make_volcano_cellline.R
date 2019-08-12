
# -- Kaili
# This script is for making volcano plot from contigency table.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/TFmotif/")

library(ggplot2)
library(gridExtra)

do_fisher <- function(line){
  name = as.character(line[1])
  a = as.integer(line[2])
  c = as.integer(line[3])
  ubi = as.integer(line[4])
  non_ubi = as.integer(line[5])
  #
  b = ubi-a
  d = non_ubi-c
  #
  p1 = fisher.test(matrix(c(a,b,c,d),2,2),alternative = "two.sided")$p.value
  p2 = fisher.test(matrix(c(a,b,c,d),2,2),alternative = "greater")$p.value
  p3 = fisher.test(matrix(c(a,b,c,d),2,2),alternative = "less")$p.value
  x = a/ubi
  y = c/non_ubi
  r = x/y
  return(c(p1,p2,p3,x,y,r))
}
####################
# all rOCRs
####################
mark = "K562"
dat = read.table(paste("TFpeak_contigency_table_",mark,".txt",sep=""), header=TRUE)

p_vector = data.frame(t(apply(dat,1,do_fisher)))
colnames(p_vector) = c("two_sided","pre_ubi","pre_non_ubi","ubi_ratio","non_ubi_ratio","FC")
p_vector$TF=dat[,1]
p_vector$log10_pvalue = -log10(p_vector$two_sided)
p_vector$log10_pvalue = ifelse(p_vector$log10_pvalue==Inf, 350, p_vector$log10_pvalue)

p_vector[p_vector$FC<0,]

p1 = p_vector[p_vector$two_sided!=0,]
p2 = p1[order(-log10(p1$two_sided),decreasing=TRUE),]


pdf("TF_enrichment_in_rOCRs.pdf")
plot(log2(p_vector$FC), p_vector$log10_pvalue, pch=20, 
     xlab = "log2FC", ylab = "-log10(p-value)", 
     main = "the enrichment of TF peaks in ubi-rOCRs")
text(log2(p2[1,]$FC), y=-log10(p2[1,]$two_sided)+6.5, rownames(p2)[1], col="blue", cex=0.8)
text(log2(p2[2,]$FC)+0.2, y=-log10(p2[2,]$two_sided)+4, rownames(p2)[2], col="blue", cex=0.8)
text(log2(p2[3,]$FC), y=-log10(p2[3,]$two_sided)+6, rownames(p2)[3], col="blue", cex=0.8)
text(log2(p2[4,]$FC), y=-log10(p2[4,]$two_sided)+6, rownames(p2)[4], col="blue", cex=0.8)
text(log2(p2[5,]$FC)-0.35, y=-log10(p2[5,]$two_sided), rownames(p2)[5], col="blue", cex=0.8)
text(log2(p2[6,]$FC), y=-log10(p2[6,]$two_sided)+6, rownames(p2)[6], col="blue", cex=0.8)
text(log2(p2[7,]$FC), y=-log10(p2[7,]$two_sided)+6, rownames(p2)[7], col="blue", cex=0.8)
text(log2(p2[8,]$FC)+0.25, y=-log10(p2[8,]$two_sided)+2, rownames(p2)[8], col="blue", cex=0.8)
text(log2(p2[9,]$FC)-0.25, y=-log10(p2[9,]$two_sided), rownames(p2)[9], col="blue", cex=0.8)
text(log2(p2[10,]$FC), y=-log10(p2[10,]$two_sided)+6, rownames(p2)[10], col="blue", cex=0.8)
text(log2(p2[11,]$FC)+0.25, y=-log10(p2[11,]$two_sided)-2, rownames(p2)[11], col="blue", cex=0.8)
dev.off()


####################
# active rOCRs
####################

######
# K562
mark="K562"
dat = read.table(paste("TFpeak_contigency_table_active_",mark,".txt",sep=""), header=TRUE)

p_vector = data.frame(t(apply(dat,1,do_fisher)))
colnames(p_vector) = c("two_sided","pre_ubi","pre_non_ubi","ubi_ratio","non_ubi_ratio","FC")
p_vector$TF=dat[,1]
p_vector$log10_pvalue = -log10(p_vector$two_sided)
p_vector$log10_pvalue = ifelse(p_vector$log10_pvalue==Inf, 350, p_vector$log10_pvalue)

p1 = p_vector[p_vector$two_sided!=0,]
p2 = p1[order(-log10(p1$two_sided),decreasing=TRUE),]

p3 = p_vector[order(p_vector$FC,decreasing=TRUE),]

pdf(paste("TFpeak_enrichment_in_active_rOCRs_",mark,".pdf",sep=""))
plot(log2(p_vector$FC), p_vector$log10_pvalue, pch=20, 
     xlab = "log2FC", ylab = "-log10(p-value)", 
     main = paste("the enrichment of TF peaks in active rOCRs\n",mark,sep=""),
     ylim=c(0,365))

text(log2(p2[1,]$FC), y=-log10(p2[1,]$two_sided)+5, p2[1,]$TF, col="blue", cex=0.8)
text(log2(p2[2,]$FC)+0.2, y=-log10(p2[2,]$two_sided)+4, p2[2,]$TF, col="blue", cex=0.8)
text(log2(p2[3,]$FC)+0.2, y=-log10(p2[3,]$two_sided), p2[3,]$TF, col="blue", cex=0.8)
text(log2(p2[4,]$FC), y=-log10(p2[4,]$two_sided)+6, p2[4,]$TF, col="blue", cex=0.8)
text(log2(p2[5,]$FC), y=-log10(p2[5,]$two_sided)+6, p2[5,]$TF, col="blue", cex=0.8)
text(log2(p2[6,]$FC)-0.7, y=-log10(p2[6,]$two_sided)-1, p2[6,]$TF, col="blue", cex=0.8)
text(log2(p2[7,]$FC)+0.4, y=-log10(p2[7,]$two_sided), p2[7,]$TF, col="blue", cex=0.8)
text(log2(p2[8,]$FC), y=-log10(p2[8,]$two_sided)+6, p2[8,]$TF, col="blue", cex=0.8)
text(log2(p2[9,]$FC), y=-log10(p2[9,]$two_sided)-6, p2[9,]$TF, col="blue", cex=0.8)
text(log2(p2[10,]$FC)+0.6, y=-log10(p2[10,]$two_sided), p2[10,]$TF, col="blue", cex=0.8)
text(log2(p2[11,]$FC), y=-log10(p2[11,]$two_sided)+6, p2[11,]$TF, col="blue", cex=0.8)

text(log2(p3[1,]$FC), y=p3[1,]$log10_pvalue+6, p3[1,]$TF, col="blue", cex=0.8)
text(log2(p3[2,]$FC), y=p3[2,]$log10_pvalue-6, p3[2,]$TF, col="blue", cex=0.8)
text(log2(p3[3,]$FC), y=p3[3,]$log10_pvalue+10, p3[3,]$TF, col="blue", cex=0.8)
text(log2(p3[4,]$FC), y=p3[4,]$log10_pvalue+6, p3[4,]$TF, col="blue", cex=0.8)
#text(log2(p3[5,]$FC), y=p3[5,]$log10_pvalue+6, p3[5,]$TF, col="blue", cex=0.8)
#text(log2(p3[6,]$FC), y=p3[6,]$log10_pvalue+6, p3[6,]$TF, col="blue", cex=0.8)
#text(log2(p3[7,]$FC), y=p3[7,]$log10_pvalue+6, p3[7,]$TF, col="blue", cex=0.8)


TF_list = c("GATA1", "GATA2", "GATA3", "GATA4", "eGFP-GATA2", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", 
            "3xFLAG-IRF2", "eGFP-IRF1", "eGFP-IRF9")
for(t in 1:length(TF_list)){
  tf = TF_list[t]
  p = p_vector[p_vector$TF==tf,]
  if(nrow(p)!=0){
    for(j in 1:nrow(p)){
      text(log2(p[j,]$FC), y=p[j,]$log10_pvalue+6, p[j,]$TF, col="purple", cex=0.8)
    }
  }
}

dev.off()

######
# HepG2
mark="HepG2"
dat = read.table(paste("TFpeak_contigency_table_active_",mark,".txt",sep=""), header=TRUE)

p_vector = data.frame(t(apply(dat,1,do_fisher)))
colnames(p_vector) = c("two_sided","pre_ubi","pre_non_ubi","ubi_ratio","non_ubi_ratio","FC")
p_vector$TF=dat[,1]
p_vector$log10_pvalue = -log10(p_vector$two_sided)
p_vector$log10_pvalue = ifelse(p_vector$log10_pvalue==Inf, 350, p_vector$log10_pvalue)

p1 = p_vector[p_vector$two_sided!=0,]
p2 = p1[order(-log10(p1$two_sided),decreasing=TRUE),]

p3 = p_vector[order(p_vector$FC,decreasing=TRUE),]

pdf(paste("TFpeak_enrichment_in_active_rOCRs_",mark,".pdf",sep=""))
plot(log2(p_vector$FC), p_vector$log10_pvalue, pch=20, 
     xlab = "log2FC", ylab = "-log10(p-value)", 
     main = paste("the enrichment of TF peaks in active rOCRs\n",mark,sep=""),
     ylim=c(0,365))

text(log2(p2[1,]$FC), y=-log10(p2[1,]$two_sided)+6, p2[1,]$TF, col="blue", cex=0.8)
text(log2(p2[2,]$FC), y=-log10(p2[2,]$two_sided)+6, p2[2,]$TF, col="blue", cex=0.8)
text(log2(p2[3,]$FC), y=-log10(p2[3,]$two_sided)+6, p2[3,]$TF, col="blue", cex=0.8)
text(log2(p2[4,]$FC), y=-log10(p2[4,]$two_sided)-6, p2[4,]$TF, col="blue", cex=0.8)
text(log2(p2[5,]$FC), y=-log10(p2[5,]$two_sided)+6, p2[5,]$TF, col="blue", cex=0.8)
text(log2(p2[6,]$FC)-0.6, y=-log10(p2[6,]$two_sided), p2[6,]$TF, col="blue", cex=0.8)
text(log2(p2[7,]$FC), y=-log10(p2[7,]$two_sided)-6, p2[7,]$TF, col="blue", cex=0.8)
text(log2(p2[8,]$FC)+0.6, y=-log10(p2[8,]$two_sided), p2[8,]$TF, col="blue", cex=0.8)
text(log2(p2[9,]$FC), y=-log10(p2[9,]$two_sided)+6, p2[9,]$TF, col="blue", cex=0.8)
text(log2(p2[10,]$FC), y=-log10(p2[10,]$two_sided)+6, p2[10,]$TF, col="blue", cex=0.8)
text(log2(p2[11,]$FC)+0.3, y=-log10(p2[11,]$two_sided), p2[11,]$TF, col="blue", cex=0.8)

text(log2(p3[1,]$FC), y=p3[1,]$log10_pvalue+6, p3[1,]$TF, col="blue", cex=0.8)
text(log2(p3[2,]$FC), y=p3[2,]$log10_pvalue-6, p3[2,]$TF, col="blue", cex=0.8)
text(log2(p3[3,]$FC), y=p3[3,]$log10_pvalue+6, p3[3,]$TF, col="blue", cex=0.8)
text(log2(p3[4,]$FC), y=p3[4,]$log10_pvalue+10, p3[4,]$TF, col="blue", cex=0.8)
text(log2(p3[5,]$FC), y=p3[5,]$log10_pvalue-6, p3[5,]$TF, col="blue", cex=0.8)
text(log2(p3[6,]$FC), y=p3[6,]$log10_pvalue-6, p3[6,]$TF, col="blue", cex=0.8)
text(log2(p3[7,]$FC), y=p3[7,]$log10_pvalue+6, p3[7,]$TF, col="blue", cex=0.8)


TF_list = c("HNF4A", "HNF1A", "FOXA2", "CEBPB", "eGFP-CEBPG", "CEBPZ", "3xFLAG-CEBPA", "3xFLAG-CEBPG",
            "eGFP-CEBPB", "FOS","eGFP-FOS")
for(t in 1:length(TF_list)){
  tf = TF_list[t]
  p = p_vector[p_vector$TF==tf,]
  if(nrow(p)!=0){
    for(j in 1:nrow(p)){
      text(log2(p[j,]$FC), y=p[j,]$log10_pvalue+6, p[j,]$TF, col="purple", cex=0.8)
    }
  }
}

dev.off()

######
# H1
mark="H1"
dat = read.table(paste("TFpeak_contigency_table_active_",mark,".txt",sep=""), header=TRUE)

p_vector = data.frame(t(apply(dat,1,do_fisher)))
colnames(p_vector) = c("two_sided","pre_ubi","pre_non_ubi","ubi_ratio","non_ubi_ratio","FC")
p_vector$TF=dat[,1]
p_vector$log10_pvalue = -log10(p_vector$two_sided)
p_vector$log10_pvalue = ifelse(p_vector$log10_pvalue==Inf, 350, p_vector$log10_pvalue)

p1 = p_vector[p_vector$two_sided!=0,]
p2 = p1[order(-log10(p1$two_sided),decreasing=TRUE),]

p3 = p_vector[order(p_vector$FC,decreasing=TRUE),]

pdf(paste("TFpeak_enrichment_in_active_rOCRs_",mark,".pdf",sep=""))
plot(log2(p_vector$FC), p_vector$log10_pvalue, pch=20, 
     xlab = "log2FC", ylab = "-log10(p-value)", 
     main = paste("the enrichment of TF peaks in active rOCRs\n",mark,sep=""),
     ylim=c(0,365))

text(log2(p2[1,]$FC), y=-log10(p2[1,]$two_sided)+6, p2[1,]$TF, col="blue", cex=0.8)
text(log2(p2[2,]$FC), y=-log10(p2[2,]$two_sided)+6, p2[2,]$TF, col="blue", cex=0.8)
text(log2(p2[3,]$FC), y=-log10(p2[3,]$two_sided)+6, p2[3,]$TF, col="blue", cex=0.8)
text(log2(p2[4,]$FC), y=-log10(p2[4,]$two_sided)-6, p2[4,]$TF, col="blue", cex=0.8)
text(log2(p2[5,]$FC), y=-log10(p2[5,]$two_sided)+6, p2[5,]$TF, col="blue", cex=0.8)
text(log2(p2[6,]$FC), y=-log10(p2[6,]$two_sided)-6, p2[6,]$TF, col="blue", cex=0.8)
text(log2(p2[7,]$FC), y=-log10(p2[7,]$two_sided)+6, p2[7,]$TF, col="blue", cex=0.8)
text(log2(p2[8,]$FC), y=-log10(p2[8,]$two_sided)+6, p2[8,]$TF, col="blue", cex=0.8)
text(log2(p2[9,]$FC), y=-log10(p2[9,]$two_sided)+6, p2[9,]$TF, col="blue", cex=0.8)
text(log2(p2[10,]$FC), y=-log10(p2[10,]$two_sided)+6, p2[10,]$TF, col="blue", cex=0.8)
text(log2(p2[11,]$FC)+0.3, y=-log10(p2[11,]$two_sided), p2[11,]$TF, col="blue", cex=0.8)

text(log2(p3[1,]$FC), y=p3[1,]$log10_pvalue+6, p3[1,]$TF, col="blue", cex=0.8)
text(log2(p3[2,]$FC), y=p3[2,]$log10_pvalue-6, p3[2,]$TF, col="blue", cex=0.8)
text(log2(p3[3,]$FC), y=p3[3,]$log10_pvalue+6, p3[3,]$TF, col="blue", cex=0.8)
text(log2(p3[4,]$FC), y=p3[4,]$log10_pvalue-6, p3[4,]$TF, col="blue", cex=0.8)
text(log2(p3[5,]$FC), y=p3[5,]$log10_pvalue-6, p3[5,]$TF, col="blue", cex=0.8)
text(log2(p3[6,]$FC), y=p3[6,]$log10_pvalue-6, p3[6,]$TF, col="blue", cex=0.8)
text(log2(p3[7,]$FC), y=p3[7,]$log10_pvalue+6, p3[7,]$TF, col="blue", cex=0.8)


TF_list = c("NANOG", "SOX6", "3xFLAG-SOX13", "3xFLAG-SOX5", "SOX13", "eGFP-ZIC2", "3xFLAG-TEAD1", 
            "3xFLAG-TEAD2","eGFP-TEAD2", "TEAD4")
for(t in 1:length(TF_list)){
  tf = TF_list[t]
  p = p_vector[p_vector$TF==tf,]
  if(nrow(p)!=0){
    for(j in 1:nrow(p)){
      text(log2(p[j,]$FC), y=p[j,]$log10_pvalue+6, p[j,]$TF, col="purple", cex=0.8)
    }
  }
}

dev.off()


####################
# active rOCRs overlap TSS
####################

######
# K562
mark="K562"
dat = dat = read.table(paste("TFpeak_contigency_table_overlapTSS_",mark,".txt",sep=""), header=TRUE)

p_vector = data.frame(t(apply(dat,1,do_fisher)))
colnames(p_vector) = c("two_sided","pre_ubi","pre_non_ubi","ubi_ratio","non_ubi_ratio","FC")
p_vector$TF=dat[,1]
p_vector$log10_pvalue = -log10(p_vector$two_sided)
p_vector$log10_pvalue = ifelse(p_vector$log10_pvalue==Inf, 350, p_vector$log10_pvalue)

p1 = p_vector[p_vector$two_sided!=0,]
p2 = p1[order(-log10(p1$two_sided),decreasing=TRUE),]

p3_0 = p_vector[p_vector$non_ubi_ratio!=0,]
p3 = p3_0[order(p3_0$FC,decreasing=TRUE),]

pdf(paste("TFpeak_enrichment_in_active_rOCRs_overlapping_TSS_",mark,".pdf",sep=""))
plot(log2(p_vector$FC), p_vector$log10_pvalue, pch=20, 
     xlab = "log2FC", ylab = "-log10(p-value)", 
     main = paste("the enrichment of TF peaks in active rOCRs\n",mark,sep=""))

text(log2(p2[1,]$FC)+0.6, y=-log10(p2[1,]$two_sided), p2[1,]$TF, col="blue", cex=0.8)
text(log2(p2[2,]$FC), y=-log10(p2[2,]$two_sided)+4, p2[2,]$TF, col="blue", cex=0.8)
text(log2(p2[3,]$FC), y=-log10(p2[3,]$two_sided)+4, p2[3,]$TF, col="blue", cex=0.8)
text(log2(p2[4,]$FC), y=-log10(p2[4,]$two_sided)+4, p2[4,]$TF, col="blue", cex=0.8)
text(log2(p2[5,]$FC), y=-log10(p2[5,]$two_sided)+4, p2[5,]$TF, col="blue", cex=0.8)
text(log2(p2[6,]$FC), y=-log10(p2[6,]$two_sided)+4, p2[6,]$TF, col="blue", cex=0.8)
text(log2(p2[7,]$FC), y=-log10(p2[7,]$two_sided)+4, p2[7,]$TF, col="blue", cex=0.8)
text(log2(p2[8,]$FC), y=-log10(p2[8,]$two_sided)-4, p2[8,]$TF, col="blue", cex=0.8)
text(log2(p2[9,]$FC), y=-log10(p2[9,]$two_sided)-4, p2[9,]$TF, col="blue", cex=0.8)
text(log2(p2[10,]$FC), y=-log10(p2[10,]$two_sided)+4, p2[10,]$TF, col="blue", cex=0.8)
text(log2(p2[11,]$FC), y=-log10(p2[11,]$two_sided)+4, p2[11,]$TF, col="blue", cex=0.8)

text(log2(p3[1,]$FC), y=p3[1,]$log10_pvalue+4, p3[1,]$TF, col="blue", cex=0.8)
text(log2(p3[2,]$FC), y=p3[2,]$log10_pvalue+4, p3[2,]$TF, col="blue", cex=0.8)
text(log2(p3[3,]$FC), y=p3[3,]$log10_pvalue+4, p3[3,]$TF, col="blue", cex=0.8)
text(log2(p3[4,]$FC), y=p3[4,]$log10_pvalue+4, p3[4,]$TF, col="blue", cex=0.8)
text(log2(p3[5,]$FC)+0.1, y=p3[5,]$log10_pvalue-1, p3[5,]$TF, col="blue", cex=0.8)
text(log2(p3[6,]$FC), y=p3[6,]$log10_pvalue-3, p3[6,]$TF, col="blue", cex=0.8)
text(log2(p3[7,]$FC), y=p3[7,]$log10_pvalue+4, p3[7,]$TF, col="blue", cex=0.8)


TF_list = c("GATA1", "GATA2", "GATA3", "GATA4", "eGFP-GATA2", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", 
            "3xFLAG-IRF2", "eGFP-IRF1", "eGFP-IRF9")
for(t in 1:length(TF_list)){
  tf = TF_list[t]
  p = p_vector[p_vector$TF==tf,]
  if(nrow(p)!=0){
    for(j in 1:nrow(p)){
      text(log2(p[j,]$FC), y=p[j,]$log10_pvalue+4, p[j,]$TF, col="purple", cex=0.8)
    }
  }
}

dev.off()


######
# HepG2
mark="HepG2"
dat = dat = read.table(paste("TFpeak_contigency_table_overlapTSS_",mark,".txt",sep=""), header=TRUE)

p_vector = data.frame(t(apply(dat,1,do_fisher)))
colnames(p_vector) = c("two_sided","pre_ubi","pre_non_ubi","ubi_ratio","non_ubi_ratio","FC")
p_vector$TF=dat[,1]
p_vector$log10_pvalue = -log10(p_vector$two_sided)
p_vector$log10_pvalue = ifelse(p_vector$log10_pvalue==Inf, 350, p_vector$log10_pvalue)

p1 = p_vector[p_vector$two_sided!=0,]
p2 = p1[order(-log10(p1$two_sided),decreasing=TRUE),]

p3_0 = p_vector[p_vector$non_ubi_ratio!=0,]
p3 = p3_0[order(p3_0$FC,decreasing=TRUE),]

pdf(paste("TFpeak_enrichment_in_active_rOCRs_overlapping_TSS_",mark,".pdf",sep=""))
plot(log2(p_vector$FC), p_vector$log10_pvalue, pch=20, 
     xlab = "log2FC", ylab = "-log10(p-value)", 
     main = paste("the enrichment of TF peaks in active rOCRs\n",mark,sep=""))

text(log2(p2[1,]$FC), y=-log10(p2[1,]$two_sided)+4, p2[1,]$TF, col="blue", cex=0.8)
text(log2(p2[2,]$FC), y=-log10(p2[2,]$two_sided)+4, p2[2,]$TF, col="blue", cex=0.8)
text(log2(p2[3,]$FC), y=-log10(p2[3,]$two_sided)+4, p2[3,]$TF, col="blue", cex=0.8)
text(log2(p2[4,]$FC)+0.7, y=-log10(p2[4,]$two_sided), p2[4,]$TF, col="blue", cex=0.8)
text(log2(p2[5,]$FC), y=-log10(p2[5,]$two_sided)+4, p2[5,]$TF, col="blue", cex=0.8)
text(log2(p2[6,]$FC)+0.5, y=-log10(p2[6,]$two_sided), p2[6,]$TF, col="blue", cex=0.8)
text(log2(p2[7,]$FC), y=-log10(p2[7,]$two_sided)+4, p2[7,]$TF, col="blue", cex=0.8)
text(log2(p2[8,]$FC), y=-log10(p2[8,]$two_sided)+4, p2[8,]$TF, col="blue", cex=0.8)
text(log2(p2[9,]$FC), y=-log10(p2[9,]$two_sided)+4, p2[9,]$TF, col="blue", cex=0.8)
text(log2(p2[10,]$FC), y=-log10(p2[10,]$two_sided)+4, p2[10,]$TF, col="blue", cex=0.8)
text(log2(p2[11,]$FC)+0.5, y=-log10(p2[11,]$two_sided), p2[11,]$TF, col="blue", cex=0.8)

text(log2(p3[1,]$FC), y=p3[1,]$log10_pvalue+4, p3[1,]$TF, col="blue", cex=0.8)
text(log2(p3[2,]$FC), y=p3[2,]$log10_pvalue+4, p3[2,]$TF, col="blue", cex=0.8)
text(log2(p3[3,]$FC), y=p3[3,]$log10_pvalue+4, p3[3,]$TF, col="blue", cex=0.8)
text(log2(p3[4,]$FC), y=p3[4,]$log10_pvalue+4, p3[4,]$TF, col="blue", cex=0.8)
text(log2(p3[5,]$FC), y=p3[5,]$log10_pvalue+4, p3[5,]$TF, col="blue", cex=0.8)
#text(log2(p3[6,]$FC), y=p3[6,]$log10_pvalue+4, p3[6,]$TF, col="blue", cex=0.8)
text(log2(p3[7,]$FC), y=p3[7,]$log10_pvalue+4, p3[7,]$TF, col="blue", cex=0.8)


TF_list = c("HNF4A", "HNF1A", "FOXA2", "CEBPB", "eGFP-CEBPG", "CEBPZ", "3xFLAG-CEBPA", "3xFLAG-CEBPG",
            "eGFP-CEBPB", "FOS","eGFP-FOS")
for(t in 1:length(TF_list)){
  tf = TF_list[t]
  p = p_vector[p_vector$TF==tf,]
  if(nrow(p)!=0){
    for(j in 1:nrow(p)){
      text(log2(p[j,]$FC), y=p[j,]$log10_pvalue+4, p[j,]$TF, col="purple", cex=0.8)
    }
  }
}

dev.off()

######
# H1
mark="H1"
dat = dat = read.table(paste("TFpeak_contigency_table_overlapTSS_",mark,".txt",sep=""), header=TRUE)

p_vector = data.frame(t(apply(dat,1,do_fisher)))
colnames(p_vector) = c("two_sided","pre_ubi","pre_non_ubi","ubi_ratio","non_ubi_ratio","FC")
p_vector$TF=dat[,1]
p_vector$log10_pvalue = -log10(p_vector$two_sided)
p_vector$log10_pvalue = ifelse(p_vector$log10_pvalue==Inf, 350, p_vector$log10_pvalue)

p1 = p_vector[p_vector$two_sided!=0,]
p2 = p1[order(-log10(p1$two_sided),decreasing=TRUE),]

p3_0 = p_vector[p_vector$non_ubi_ratio!=0,]
p3 = p3_0[order(p3_0$FC,decreasing=TRUE),]

pdf(paste("TFpeak_enrichment_in_active_rOCRs_overlapping_TSS_",mark,".pdf",sep=""))
plot(log2(p_vector$FC), p_vector$log10_pvalue, pch=20, 
     xlab = "log2FC", ylab = "-log10(p-value)", 
     main = paste("the enrichment of TF peaks in active rOCRs\n",mark,sep=""))

text(log2(p2[1,]$FC), y=-log10(p2[1,]$two_sided)+3, p2[1,]$TF, col="blue", cex=0.8)
text(log2(p2[2,]$FC), y=-log10(p2[2,]$two_sided)+3, p2[2,]$TF, col="blue", cex=0.8)
text(log2(p2[3,]$FC), y=-log10(p2[3,]$two_sided)+3, p2[3,]$TF, col="blue", cex=0.8)
text(log2(p2[4,]$FC), y=-log10(p2[4,]$two_sided)+3, p2[4,]$TF, col="blue", cex=0.8)
text(log2(p2[5,]$FC), y=-log10(p2[5,]$two_sided)+3, p2[5,]$TF, col="blue", cex=0.8)
text(log2(p2[6,]$FC), y=-log10(p2[6,]$two_sided)+3, p2[6,]$TF, col="blue", cex=0.8)
text(log2(p2[7,]$FC), y=-log10(p2[7,]$two_sided)+3, p2[7,]$TF, col="blue", cex=0.8)
text(log2(p2[8,]$FC), y=-log10(p2[8,]$two_sided)+3, p2[8,]$TF, col="blue", cex=0.8)
text(log2(p2[9,]$FC), y=-log10(p2[9,]$two_sided)-3, p2[9,]$TF, col="blue", cex=0.8)
text(log2(p2[10,]$FC)+0.2, y=-log10(p2[10,]$two_sided), p2[10,]$TF, col="blue", cex=0.8)
text(log2(p2[11,]$FC), y=-log10(p2[11,]$two_sided)+3, p2[11,]$TF, col="blue", cex=0.8)

text(log2(p3[1,]$FC), y=p3[1,]$log10_pvalue+3, p3[1,]$TF, col="blue", cex=0.8)
text(log2(p3[2,]$FC), y=p3[2,]$log10_pvalue+3, p3[2,]$TF, col="blue", cex=0.8)
text(log2(p3[3,]$FC), y=p3[3,]$log10_pvalue+3, p3[3,]$TF, col="blue", cex=0.8)
text(log2(p3[4,]$FC), y=p3[4,]$log10_pvalue+3, p3[4,]$TF, col="blue", cex=0.8)
#text(log2(p3[5,]$FC)+0.1, y=p3[5,]$log10_pvalue-1, p3[5,]$TF, col="blue", cex=0.8)
text(log2(p3[6,]$FC), y=p3[6,]$log10_pvalue+3, p3[6,]$TF, col="blue", cex=0.8)
text(log2(p3[7,]$FC), y=p3[7,]$log10_pvalue+3, p3[7,]$TF, col="blue", cex=0.8)


TF_list = c("NANOG", "SOX6", "3xFLAG-SOX13", "3xFLAG-SOX5", "SOX13", "eGFP-ZIC2", "3xFLAG-TEAD1", 
            "3xFLAG-TEAD2","eGFP-TEAD2", "TEAD4")
for(t in 1:length(TF_list)){
  tf = TF_list[t]
  p = p_vector[p_vector$TF==tf,]
  if(nrow(p)!=0){
    for(j in 1:nrow(p)){
      text(log2(p[j,]$FC), y=p[j,]$log10_pvalue+3, p[j,]$TF, col="purple", cex=0.8)
    }
  }
}

dev.off()

####################
