
# -- Kaili
# This script is for making figures and statistic test for CTCF states, CTCF peaks and CTCF motif.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/")

data = read.table("liver_14.5_day_CTCF_peak_signal.txt")

####################
# histogram of CTCF signals in CTCF peaks
####################
summary(data[,5])
hist(data[,5], breaks=100, xlab = "CTCF signal", main = "CTCF signal in CTCF peaks\nliver_e14.5")
hist(data[,5], breaks=100, xlab = "CTCF signal", main = "CTCF signal in CTCF peaks\nliver_e14.5",
     ylim = c(0,5000))



####################
# fisher exact test
####################
dat = rbind(c(13812,39299), c(1051, 245838))
chisq.test(dat, correct = F)
fisher.test(dat)

dat = rbind(c(911,4556), c(2, 15706))
chisq.test(dat, correct = F)
fisher.test(dat)

dat = rbind(c(1655,17515), c(0, 30262))
chisq.test(dat, correct = F)
fisher.test(dat)

####################