
# -- Kaili
# This script is for making figures and statistic test for CTCF states, CTCF peaks and CTCF motif.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/")

####################
# histogram of CTCF signals in CTCF peaks
####################
data = read.table("liver_14.5_day_CTCF_peak_signal.txt")
summary(data[,5])
pdf("CTCFsignal_CTCFpeaks_liver_14.5_histogram.pdf")
hist(data[,5], breaks=100, xlab = "CTCF signal", main = "CTCF signal in CTCF peaks\nliver_e14.5", 
     freq = F)
lines(c(32,32),c(0,0.1), col="red", lwd=2, lty=2)
dev.off()


####################
# fisher exact test
####################
# CTCF peaks VS. CTCF motif
dat = rbind(c(16158,15506), c(4566, 14679))
chisq.test(dat, correct = F)
fisher.test(dat)

# CTCF peaks in CTCF states
dat = rbind(c(20088,23685), c(636,6500))
chisq.test(dat, correct = F)
fisher.test(dat)

# CTCF states in CTCF peaks
dat = rbind(c(40148,36040), c(2152,40489))
chisq.test(dat, correct = F)
fisher.test(dat)

# CTCF states in CTCF motif
dat = rbind(c(40147,74039), c(98025,13025000))
chisq.test(dat, correct = F)
fisher.test(dat)

####################
# CTCF signal boxplot
####################
library(ggplot2)

dat = data.frame(read.table("liver_14.5_day_ctcf_signal.txt"))
colnames(dat) = c("name", "signal", "group")
ggplot(dat, aes(x=group, y=signal, fill=group)) + geom_boxplot() +
  labs(title="CTCF signal (liver 14.5day)")
ggsave("liver_14.5_day_ctcf_signal.pdf", width=8, height = 7)

dat[,2] = log10(dat$signal+0.1)
ggplot(dat, aes(x=group, y=signal, fill=group)) + geom_boxplot() +
  labs(title="CTCF signal (liver 14.5day)") + ylab("log10(signal+0.1)")
ggsave("liver_14.5_day_ctcf_signal_log.pdf", width=8, height = 7)

####################
# running cut-off for CTCF signal in peaks
####################

with_motif = read.table("mm10_ctcf_peak_cutoff_with_ctcf_motif_num.txt")
pdf("mm10_ctcf_peak_cutoff_with_ctcf_motif_num.pdf")
bar <- barplot(with_motif$V2, col = "grey", border = "grey", 
        ylab = "num of high signal CTCF peaks with CTCF motif",
        xlab = "cut-off for high signal CTCF peaks",
        main = "high signal CTCF peaks contains CTCF motif")
axis(side=1, at = bar[c(1,50,100,150,200,250)], labels =c("0","50","100","150","200","250"))
lines(c(32,32), c(0,45000), col="red", lty=2, lwd=2)
dev.off()

with_state = read.table("mm10_ctcf_peak_cutoff_with_ctcf_state_num.txt")
pdf("mm10_ctcf_peak_cutoff_with_ctcf_state_num.pdf")
bar <- barplot(with_state$V2, col = "grey", border = "grey", 
               ylab = "num of high signal CTCF peaks with CTCF state",
               xlab = "cut-off for high signal CTCF peaks",
               main = "high signal CTCF peaks VS CTCF state")
axis(side=1, at = bar[c(1,50,100,150,200,250)], labels =c("0","50","100","150","200","250"))
lines(c(32,32), c(0,50000), col="red", lty=2, lwd=2)
dev.off()
####################
####################
