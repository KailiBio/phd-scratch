
# -- Kaili
# This script is for doing linear regression for validating CTCF prediction power.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/ctcfstate_with_ctcfdata/")


dat = read.table("lung0_13signal_random100k.txt", header = TRUE, row.names = 1)

##################
# CTCF ~ 10 marks
##################
model0 = lm(dat$CTCF~dat$ATAC + dat$DNAme + dat$H3K4me1 + dat$H3K4me2 + dat$H3K4me3 + dat$H3K9me3 +
              dat$H3K9ac + dat$H3K27me3 + dat$H3K27ac + dat$H3K36me3)
summary(model0)$adj.r.squared

##################
# CTCF ~ 10 marks + CTCFmotif + CTCFaverage
##################
model1 = lm(dat$CTCF~dat$ATAC + dat$DNAme + dat$H3K4me1 + dat$H3K4me2 + dat$H3K4me3 + dat$H3K9me3 +
              dat$H3K9ac + dat$H3K27me3 + dat$H3K27ac + dat$H3K36me3 + dat$CTCFmotif + dat$CTCFaverage)
summary(model1)$adj.r.squared

##################
# CTCF ~ 10 marks + CTCFmotif
##################
model2 = lm(dat$CTCF~dat$ATAC + dat$DNAme + dat$H3K4me1 + dat$H3K4me2 + dat$H3K4me3 + dat$H3K9me3 +
              dat$H3K9ac + dat$H3K27me3 + dat$H3K27ac + dat$H3K36me3 + dat$CTCFmotif)
summary(model2)$adj.r.squared

##################
# CTCF ~ 10 marks + CTCFaverage
##################
model3 = lm(dat$CTCF~dat$ATAC + dat$DNAme + dat$H3K4me1 + dat$H3K4me2 + dat$H3K4me3 + dat$H3K9me3 +
              dat$H3K9ac + dat$H3K27me3 + dat$H3K27ac + dat$H3K36me3 + dat$CTCFaverage)
summary(model3)$adj.r.squared

##################
# CTCF ~ CTCFaverage
##################
model4 = lm(dat$CTCF~dat$CTCFaverage)
summary(model4)$adj.r.squared

##################
# CTCF ~ CTCFmotif
##################
model5 = lm(dat$CTCF~dat$CTCFmotif)
summary(model5)$adj.r.squared

##################
# CTCF ~ DNAme + CTCFmotif
##################
model6 = lm(dat$CTCF~dat$DNAme + dat$CTCFmotif)
summary(model6)$adj.r.squared

##################
# CTCFbinary ~ DNAme + CTCFmotif
##################
CTCFbinary = as.integer(dat$CTCF>10)
#
model7 = lm(CTCFbinary~dat$DNAme + dat$CTCFmotif)
summary(model7)$adj.r.squared

##################
# CTCFpeakSummit ~ DNAme
##################
model8 = lm(dat$CTCFpeakSummit~dat$DNAme + dat$CTCFmotif)
summary(model8)$adj.r.squared
##################