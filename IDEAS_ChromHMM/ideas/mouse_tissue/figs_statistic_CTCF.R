
# -- Kaili
# This script is for making figures and statistic test for CTCF states, CTCF peaks and CTCF motif.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/")

####################
# histogram of CTCF signals in CTCF peaks
####################
data = read.table("liver_14.5_day_CTCF_peak_signal.txt")
summary(data[,5])
hist(data[,5], breaks=100, xlab = "CTCF signal", main = "CTCF signal in CTCF peaks\nliver_e14.5")
hist(data[,5], breaks=100, xlab = "CTCF signal", main = "CTCF signal in CTCF peaks\nliver_e14.5",
     ylim = c(0,5000))


####################
# fisher exact test
####################
dat = rbind(c(12273,36518), c(2590, 148619))
chisq.test(dat, correct = F)
fisher.test(dat)

dat = rbind(c(911,4556), c(2, 15706))
chisq.test(dat, correct = F)
fisher.test(dat)

dat = rbind(c(1655,17515), c(0, 30262))
chisq.test(dat, correct = F)
fisher.test(dat)

####################
# CTCF signal boxplot
####################
dat = data.frame(read.table("liver_14.5_day_chr1_ctcf_signal.txt"))
colnames(dat) = c("name", "signal", "group")
library(ggplot2)
ggplot(dat, aes(x=group, y=signal, fill=group)) + geom_boxplot() +
  labs(title="CTCF signal (liver 14.5day)")
ggsave("liver_14.5_day_chr1_ctcf_signal.pdf", width=8, height = 7)
####################
# running cut-off for CTCF signal in peaks
####################

with_motif = read.table("mm10_ctcf_peak_cutoff_with_ctcf_motif_num.txt")
pdf("mm10_ctcf_peak_cutoff_with_ctcf_motif_num.pdf")
bar <- barplot(with_motif$V2, col = "grey", border = "grey", 
        ylab = "num of high signal CTCF peaks with CTCF motif",
        xlab = "cut-off for high signal CTCF peaks",
        main = "high signal CTCF peaks contains CTCF motif")
axis(side=1, at = bar[c(1,50,100,150,200)], labels =c("0","50","100","150","200"))
lines(c(50,50), c(0,40000), col="red", lty=2, lwd=2)
dev.off()

with_state = read.table("mm10_ctcf_peak_cutoff_with_ctcf_state_num.txt")
pdf("mm10_ctcf_peak_cutoff_with_ctcf_state_num.pdf")
bar <- barplot(with_state$V2, col = "grey", border = "grey", 
               ylab = "num of high signal CTCF peaks with CTCF state",
               xlab = "cut-off for high signal CTCF peaks",
               main = "high signal CTCF peaks VS CTCF state")
axis(side=1, at = bar[c(1,50,100,150,200)], labels =c("0","50","100","150","200"))
lines(c(50,50), c(0,4000), col="red", lty=2, lwd=2)
dev.off()
####################
####################
