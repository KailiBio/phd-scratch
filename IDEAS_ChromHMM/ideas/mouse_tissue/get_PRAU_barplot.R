
# -- Kaili
# get PRAU barplot for all.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/state_ranked_AUC/")

library(ggplot2) 
library(gridExtra)

##################
# all imputation ways
##################
dat = data.frame(c(0.6, 0.49, 0.56, 0.66, 0.4, 0.59, 0.65, 0.71),
                 c("withRealData", "majority_vote", "CTCFaverageForMissing", "CTCFmotif+CTCFaverageForMissing",
                   "imputation", "imputation+CTCFmotif", "imputation+CTCFaverage", 
                   "imputation+CTCFmotif+CTCFaverage"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)

ggplot(dat, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity",  width = 0.5, 
           fill=c("#06DA93","#06DA93","#FFCD00","#FFCD00", "#FFCD00", "#FFCD00","#FFAAAA", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="liver_14.5")

ggsave("PRAU_for_all_liver_14.5.pdf", width=8, height = 5)
ggsave("PRAU_for_all_liver_14.5.png", width=8, height = 5)

######
dat = data.frame(c(0.63, 0.56, 0.6, 0.71, 0.41, 0.63, 0.71, 0.77),
                 c("withRealData", "majority_vote", "CTCFaverageForMissing", "CTCFmotif+CTCFaverageForMissing",
                   "imputation", "imputation+CTCFmotif", "imputation+CTCFaverage", 
                   "imputation+CTCFmotif+CTCFaverage"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)

ggplot(dat, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity",  width = 0.5, 
           fill=c("#06DA93","#06DA93","#FFCD00","#FFCD00", "#FFCD00", "#FFCD00","#FFAAAA", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="lung_14.5")

ggsave("PRAU_for_all_lung_14.5.pdf", width=8, height = 5)
ggsave("PRAU_for_all_lung_14.5.png", width=8, height = 5)


######
#----- unfinished
dat = data.frame(c(0.63, 0.56, 0.6, 0.71, 0.41, 0.63, 0.71, 0.77),
                 c("withRealData", "majority_vote", "CTCFaverageForMissing", "CTCFmotif+CTCFaverageForMissing",
                   "imputation", "imputation+CTCFmotif", "imputation+CTCFaverage", 
                   "imputation+CTCFmotif+CTCFaverage"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)

ggplot(dat, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity",  width = 0.5, 
           fill=c("#06DA93","#06DA93","#FFCD00","#FFCD00", "#FFCD00", "#FFCD00","#FFAAAA", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="lung_0")

ggsave("PRAU_for_all_lung_0.pdf", width=8, height = 5)
ggsave("PRAU_for_all_lung_0.png", width=8, height = 5)

##################
# PRAU for different para-fold
# liver 14.5
##################
dat = data.frame(c(0.7128989, 0.7131602, 0.7128445, 0.7112927, 0.7021849),
                 c("1-fold","100-fold","10k","1m", "directly-imputation"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)

ggplot(dat, aes(x = group, y=PRAU)) + geom_bar(stat = "identity", fill= "#AEAFAE", width = 0.5) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("")

ggsave("PRAU_for_different_fold_liver_14.5.pdf", width=7, height = 5)
ggsave("PRAU_for_different_fold_liver_14.5.png", width=7, height = 5)


##################
# PRAU for different para-fold
# lung 14.5
##################
dat = data.frame(c(0.7718592, 0.7725415, 0.7721105, 0.7680579, 0.7573383),
                 c("1-fold","100-fold","10k","1m","directly-imputation"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)

ggplot(dat, aes(x = group, y=PRAU)) + geom_bar(stat = "identity", fill= "#AEAFAE", width = 0.5) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("")

ggsave("PRAU_for_different_fold_lung_14.5.pdf", width=7, height = 5)
ggsave("PRAU_for_different_fold_lung_14.5.png", width=7, height = 5)

##################
# 8impute3
##################
## forebrain
dat = data.frame(c(0.6053042, 0.3230964, 0.5710588, 0.6642508, 0.7194593),
                 c("withRealData","imputation","imputation+CTCFmotif", "imputation+CTCFaverage", 
                   "imputation+CTCFmotif+CTCFaverage"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)

ggplot(dat, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#FFCD00", "#FFCD00", "#FFCD00", "#FFCD00", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="forebrain")

ggsave("PRAU_for_8impute11_forebrain_0.pdf", width=6, height = 8)
ggsave("PRAU_for_8impute11_forebrain_0.png", width=6, height = 8)

## midbrain
dat = data.frame(c(0.5769864, 0.3573072, 0.5486568, 0.653502, 0.707268),
                 c("withRealData","imputation","imputation+CTCFmotif", "imputation+CTCFaverage", 
                   "imputation+CTCFmotif+CTCFaverage"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)


ggplot(dat, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#FFCD00", "#FFCD00", "#FFCD00", "#FFCD00", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="midbrain")

ggsave("PRAU_for_8impute11_midbrain_0.pdf", width=6, height = 8)
ggsave("PRAU_for_8impute11_midbrain_0.png", width=6, height = 8)

## hindbrain
dat = data.frame(c(0.5980756, 0.3695073, 0.5609002, 0.6670622, 0.7200354),
                 c("withRealData","imputation","imputation+CTCFmotif", "imputation+CTCFaverage", 
                   "imputation+CTCFmotif+CTCFaverage"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)


ggplot(dat, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#FFCD00", "#FFCD00", "#FFCD00", "#FFCD00", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="hindbrain")

ggsave("PRAU_for_8impute11_hindbrain_0.pdf", width=6, height = 8)
ggsave("PRAU_for_8impute11_hindbrain_0.png", width=6, height = 8)
##################
# compare with Avocado & ChromImpute
##################

dat1 = data.frame(c(0.6941439, 0.6843176, 0.5952582, 0.71),
                 c("ChromImpute","ChomImpute_motif_average", "Avocado","IDEAS_motif_average"))
colnames(dat1) = c("PRAU", "group")
dat1$group = factor(dat1$group)

dat2 = data.frame(c(0.7135386, 0.7109106, 0.667829, 0.77),
                  c("ChromImpute","ChomImpute_motif_average", "Avocado","IDEAS_motif_average"))
colnames(dat2) = c("PRAU", "group")
dat2$group = factor(dat2$group)


p1 = ggplot(dat1, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#D642CA", "#FF7903", "#FF7903", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="liver_14.5")

p2 = ggplot(dat2, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#D642CA", "#FF7903", "#FF7903",  "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="lung_14.5")

p = grid.arrange(p1,p2, ncol=2)
ggsave("PRAU_for_chromimpute_avocado.pdf", p, width=9, height = 5)
ggsave("PRAU_for_chromimpute_avocado.png", p, width=9, height = 5)


##################
# compare with Avocado & ChromImpute
# in DHS-center bins 
##################
dat1 = data.frame(c(0.8057016, 0.7966459, 0.7087602, 0.71),
                  c("ChromImpute","ChomImpute_motif_average", "Avocado","IDEAS_motif_average"))
colnames(dat1) = c("PRAU", "group")
dat1$group = factor(dat1$group)

dat2 = data.frame(c(0.8406349, 0.8388087, 0.8130574, 0.77),
                  c("ChromImpute","ChomImpute_motif_average", "Avocado","IDEAS_motif_average"))
colnames(dat2) = c("PRAU", "group")
dat2$group = factor(dat2$group)


p1 = ggplot(dat1, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#D642CA", "#FF7903", "#FF7903", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="liver_14.5 (DHS-center bins)")

p2 = ggplot(dat2, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#D642CA", "#FF7903", "#FF7903",  "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="lung_14.5 (DHS-center bins)")

p = grid.arrange(p1,p2, ncol=2)
ggsave("PRAU_for_chromimpute_avocado_DHScenter_bins.pdf", p, width=9, height = 5)
ggsave("PRAU_for_chromimpute_avocado_DHScenter_bins.png", p, width=9, height = 5)

##################
# ChromImpute 8-3
##################
## forebrain
dat = data.frame(c(0.7194593, 0.6961614, 0.7085652, 0.8173914, 0.8298819),
                 c("IDEAS", "ChromImpute_25bp", "ChromImpute_motif_average_25bp", "ChromImpute_dhs",
                   "ChromImpute_motif_average_dhs"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)

ggplot(dat, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#FFCD00", "#FFCD00", "#FFCD00", "#FFCD00", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="forebrain")

ggsave("PRAU_for_8impute11_forebrain_0_chromimpute.pdf", width=6, height = 8)
ggsave("PRAU_for_8impute11_forebrain_0_chromimpute.png", width=6, height = 8)

## midbrain
dat = data.frame(c(0.707268, 0.7003482, 0.7195917, 0.7842273, 0.8016976),
                 c("IDEAS", "ChromImpute_25bp", "ChromImpute_motif_average_25bp", "ChromImpute_dhs",
                   "ChromImpute_motif_average_dhs"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)

ggplot(dat, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#FFCD00", "#FFCD00", "#FFCD00", "#FFCD00", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="midbrain")

ggsave("PRAU_for_8impute11_midbrain_0_chromimpute.pdf", width=6, height = 8)
ggsave("PRAU_for_8impute11_midbrain_0_chromimpute.png", width=6, height = 8)

## hindbrain
dat = data.frame(c(0.7200354, 0.6811362, 0.6950856, 0.8065607, 0.8185678),
                 c("IDEAS", "ChromImpute_25bp", "ChromImpute_motif_average_25bp", "ChromImpute_dhs",
                   "ChromImpute_motif_average_dhs"))
colnames(dat) = c("PRAU", "group")
dat$group = factor(dat$group)

ggplot(dat, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#FFCD00", "#FFCD00", "#FFCD00", "#FFCD00", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="hindbrain")

ggsave("PRAU_for_8impute11_hindbrain_0_chromimpute.pdf", width=6, height = 8)
ggsave("PRAU_for_8impute11_hindbrain_0_chromimpute.png", width=6, height = 8)

##################
# imputation for 66 samples
##################

##################
# Aug 07
# compare with Avocado & ChromImpute
# in DHS-center bins & with same resolution signal
##################
dat1 = data.frame(c(0.763494, 0.7101903, 0.71),
                  c("ChromImpute", "Avocado","IDEAS_motif_average"))
colnames(dat1) = c("PRAU", "group")
dat1$group = factor(dat1$group)

dat2 = data.frame(c(0.8002291, 0.8000357, 0.77),
                  c("ChromImpute","Avocado","IDEAS_motif_average"))
colnames(dat2) = c("PRAU", "group")
dat2$group = factor(dat2$group)


p1 = ggplot(dat1, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#D642CA", "#FF7903",  "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="liver_14.5 (DHS-center bins & same resolution)")

p2 = ggplot(dat2, aes(x = group, y=PRAU)) + 
  geom_bar(stat = "identity", width = 0.5, fill=c("#D642CA", "#FF7903", "#FF0000")) +
  geom_text(aes(label = PRAU), vjust=0) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1),
                          title = element_text(face="bold", size=12)) +
  xlab("") + labs(title="lung_14.5 (DHS-center bins & same resolution)")

p = grid.arrange(p1,p2, ncol=2)
ggsave("PRAU_for_chromimpute_avocado_DHScenter_bins_sameResolution.pdf", p, width=9, height = 5)
ggsave("PRAU_for_chromimpute_avocado_DHScenter_bins_sameResolution.png", p, width=9, height = 5)

##################





##################
##################


##################