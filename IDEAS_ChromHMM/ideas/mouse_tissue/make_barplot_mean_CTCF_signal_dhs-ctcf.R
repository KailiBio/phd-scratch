
# -- Kaili
# This script is for making barplot of mean signal of CTCF in each sample.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/state_CTCF_signal/")


data = read.table("state_CTCF_signal_mean.txt")
colnames(data) = c("state", "sample", "mean")

ggplot(data, aes(x=state, y=mean)) + geom_bar(stat = "identity") +
  facet_wrap(~sample, scales = "free_x") +
  ylab("mean of CTCF signal") +
  theme(title = element_text(face="bold", size=12))

ggsave("dhs_ctcf_mean_signal_in_11samples.pdf")
ggsave("dhs_ctcf_mean_signal_in_11samples.png")
