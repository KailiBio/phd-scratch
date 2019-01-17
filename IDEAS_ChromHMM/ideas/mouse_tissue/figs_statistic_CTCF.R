
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
dat = rbind(c(20386,24267), c(338,5918))
chisq.test(dat, correct = F)
fisher.test(dat)

# CTCF states in CTCF peaks
dat = rbind(c(40951,36839), c(1349,39690))
chisq.test(dat, correct = F)
fisher.test(dat)

# CTCF states in CTCF motif
dat = rbind(c(40557,76071), c(97615,13022968))
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
with_motif = transform(with_motif, percent=round(with_motif[,2]/with_motif[,3],3)*100)
pdf("mm10_ctcf_peak_cutoff_with_ctcf_motif_num.pdf")
bar <- barplot(with_motif$percent, col = "grey", border = "grey", width=1, space = 0,
        ylab = "percentage of CTCF peaks that with CTCF motif",
        xlab = "CTCF signal of CTCF peaks",
        main = "proportion of CTCF peaks that have CTCF motif\n(with different signal cut-off)")
axis(side=1, at = bar[c(1,50,100,150,200,247)], labels =c("0","50","100","150","200","247"))
#lines(c(32,32), c(0,45000), col="red", lty=2, lwd=2)
dev.off()

with_state = read.table("mm10_ctcf_peak_cutoff_with_ctcf_state_num.txt")
with_state = transform(with_state, percent=round(with_state[,2]/with_state[,3],3)*100)
pdf("mm10_ctcf_peak_cutoff_with_ctcf_state_num.pdf")
bar <- barplot(with_state$percent, col = "grey", border = "grey", width=1, space = 0,
               ylab = "percentage of CTCF peaks that are CTCF state",
               xlab = "CTCF signal of CTCF peaks",
               main = "proportion of CTCF peaks that are CTCF state\n(with different signal cut-off)")
axis(side=1, at = bar[c(1,50,100,150,200,247)], labels =c("0","50","100","150","200","247"))
lines(c(12,12), c(0,50000), col="red", lty=2, lwd=2)
dev.off()
####################
# CTCF signal with CTCF peaks boxplot
####################
library(ggplot2)

data1 = data.frame(read.table("liver_14.5_day_ctcf_peak_rep1.2_signal.txt"))
data2 = data.frame(read.table("liver_14.5_day_ctcf_peak_rep2_signal.txt"))

colnames(data1) = c("name", "signal", "group")
ggplot(data1, aes(x=group, y=signal, fill=group)) + geom_boxplot() +
  labs(title="CTCF signal (liver 14.5day)\n(rep1,2)")
ggsave("liver_14.5_day_ctcf_peak_signal_rep1.2_log.pdf", width=8, height = 7)

colnames(data2) = c("name", "signal", "group")
ggplot(data2, aes(x=group, y=signal, fill=group)) + geom_boxplot() +
  labs(title="CTCF signal (liver 14.5day)\n(rep2)")
ggsave("liver_14.5_day_ctcf_peak_signal_rep2_log.pdf", width=8, height = 7)
####################
# CTCF peak length
####################
len = read.table("liver_14.5_day_CTCF_peak_length.bed")
summary(len[,7])
pdf("liver14.5_ctcf_peak_length.pdf")
hist(len[,7], breaks=78.5:778.5, freq = F, xlab = "length of CTCF peaks", 
     main = "length distribution of CTCF peaks\n(n=50,909)")
dev.off()

pdf("liver14.5_ctcf_peak_length_2.pdf")
hist(len[,7], breaks=78.5:778.5, freq = F, xlab = "length of CTCF peaks", 
     main = "length distribution of CTCF peaks\n(n=50,909)", ylim = c(0,0.005))
dev.off()

####################
# CTCF signal in CTCF states (histogram)
####################
ctcf_state_signal=read.table("mm10_ctcf_state_liver14.5_signal.txt", row.names = 1)
summary(ctcf_state_signal[,1])

pdf("mm10_ctcf_state_liver14.5_signal.pdf")
hist(ctcf_state_signal[,1], breaks=-0.5:243.5, freq = F, col="grey", 
     xlab = "CTCF signal", main = "distribution of CTCF signal in CTCF states")
dev.off()

####################
# running cut-off for CTCF states having CTCF peaks
####################
state_with_peak = read.table("mm10_ctcf_state_cutoff_with_ctcf_peak_num.txt")
state_with_peak = transform(state_with_peak, percent=round(state_with_peak[,2]/state_with_peak[,3],3)*100)
pdf("mm10_ctcf_state_cutoff_with_ctcf_peak_num.pdf")
bar <- barplot(state_with_peak$percent, col = "grey", border = "grey", width=1, space=0,
               ylab = "percentage of CTCF states that are CTCF peaks",
               xlab = "CTCF signal of CTCF states",
               main = "proportion of CTCF states that are CTCF peaks\n(with different signal cut-off)")
axis(side=1, at = bar[c(1,50,100,150,200,243)], labels =c("0","50","100","150","200","243"))
lines(c(43,43),c(0,100), lwd=2, lty=2,col="red")
dev.off()

####################
# CTCF states overlapped matrix
####################
matrix = read.table("ctcf_states_overlapped.txt", header = TRUE)
rownames(matrix) = colnames(matrix)
matrix_union = read.table("ctcf_states_union.txt", header = TRUE)
rownames(matrix_union) = colnames(matrix_union)


library(pheatmap)
pdf("CTCF_states_between_all_samples.pdf")
pheatmap(matrix/matrix_union, display_numbers = TRUE, main = "CTCF states in each cell-types")
dev.off()
####################




