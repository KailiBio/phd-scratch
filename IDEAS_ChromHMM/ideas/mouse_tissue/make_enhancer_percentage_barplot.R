
# -- Kaili
# This script is for making barplot to show the percentage of bins in Arjan's enhancer list.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/pool-bins/")

library(ggplot2)

data = read.table("bins_percentage_in_enhancer_list.txt")
colnames(data) = c("state","num","total","percentage")

ggplot(data, aes(x=state, y=percentage)) + geom_bar(stat="identity") +
  theme_minimal() + theme(title=element_text(face="bold", size=12)) +
  labs(title="percentage of bins in Arjan's enhancer list") + 
  scale_x_continuous(breaks=seq(0,46,by=2), labels=seq(0,46,by=2))
ggsave("state_percentage_in_enhancer_list_barplot.pdf", width=12)
ggsave("state_percentage_in_enhancer_list_barplot.png",width=12)
