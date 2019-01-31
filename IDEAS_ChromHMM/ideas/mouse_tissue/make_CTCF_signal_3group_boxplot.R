
# -- Kaili
# This script is for making states CTCF signal (with peak/motif).

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/state_CTCF_signal/")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

####################
# peak
####################
peak1 = read.table("liver14.5_11states_peak_with_signal.txt")
peak2 = read.table("liver14.5_11states_peak_without_signal.txt")
peak3 = read.table("liver14.5_11states_peak_non_with_signal.txt")

matrix_peak = data.frame(rbind(cbind(peak1, type="11states_withPeak"),
                          cbind(peak2, type="11states_withoutPeak"),
                          cbind(peak3, type="no11states_withPeak")))
colnames(matrix_peak) = c("bin", "signal", "type")

ggplot(matrix_peak, aes(x=type, y=signal, fill=type, group=type)) +
  geom_boxplot(width=0.3) +
  theme_minimal() +
  theme(title = element_text(face="bold", size=12),
        axis.text = element_text(face="bold", size=10)) +
  ylab("CTCF signal") + xlab("")

ggsave("CTCF_signal_11states_3group_peak.pdf")
ggsave("CTCF_signal_11states_3group_peak.png")

####################
# motif
####################
motif1 = read.table("liver14.5_11states_motif_with_signal.txt")
motif2 = read.table("liver14.5_11states_motif_without_signal.txt")
motif3 = read.table("liver14.5_11states_motif_non_with_signal.txt")

matrix_motif = data.frame(rbind(cbind(motif1,type="11states_withMotif"),
                               cbind(motif2,type="11states_withoutMotif"),
                               cbind(motif3,type="no11states_withMotif")))
colnames(matrix_motif) = c("bin", "signal", "type")

ggplot(matrix_motif, aes(x=type, y=signal, fill=type, group=type)) +
  geom_boxplot(width=0.3) +
  theme_minimal() +
  theme(title = element_text(face="bold", size=12),
        axis.text = element_text(face="bold", size=10)) +
  ylab("CTCF signal") + xlab("")

ggsave("CTCF_signal_11states_3group_motif.pdf")
ggsave("CTCF_signal_11states_3group_motif.png")

####################
# peak, 6 states
####################
peak1 = read.table("liver14.5_6states_peak_with_signal.txt")
peak2 = read.table("liver14.5_6states_peak_without_signal.txt")
peak3 = read.table("liver14.5_6states_peak_non_with_signal.txt")

matrix_peak = data.frame(rbind(cbind(peak1, type="6states_withPeak"),
                               cbind(peak2, type="6states_withoutPeak"),
                               cbind(peak3, type="no6states_withPeak")))
colnames(matrix_peak) = c("bin", "signal", "type")

ggplot(matrix_peak, aes(x=type, y=signal, fill=type, group=type)) +
  geom_boxplot(width=0.3) +
  theme_minimal() +
  theme(title = element_text(face="bold", size=12),
        axis.text = element_text(face="bold", size=10)) +
  ylab("CTCF signal") + xlab("") +
  labs(title="6 states: 21,29,30,40,42,11")

ggsave("CTCF_signal_6states_3group_peak.pdf")
ggsave("CTCF_signal_6states_3group_peak.png")

####################
# motif, 6 states
####################
motif1 = read.table("liver14.5_6states_motif_with_signal.txt")
motif2 = read.table("liver14.5_6states_motif_without_signal.txt")
motif3 = read.table("liver14.5_6states_motif_non_with_signal.txt")

matrix_motif = data.frame(rbind(cbind(motif1,type="6states_withMotif"),
                                cbind(motif2,type="6states_withoutMotif"),
                                cbind(motif3,type="no6states_withMotif")))
colnames(matrix_motif) = c("bin", "signal", "type")

ggplot(matrix_motif, aes(x=type, y=signal, fill=type, group=type)) +
  geom_boxplot(width=0.3) +
  theme_minimal() +
  theme(title = element_text(face="bold", size=12),
        axis.text = element_text(face="bold", size=10)) +
  ylab("CTCF signal") + xlab("") +
  labs(title="6 states: 21,29,30,40,42,11")

ggsave("CTCF_signal_6states_3group_motif.pdf")
ggsave("CTCF_signal_6states_3group_motif.png")
####################
