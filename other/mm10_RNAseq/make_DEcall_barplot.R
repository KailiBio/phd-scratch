
# -- Kaili
# This script is for making comparison barplot for DE-calls.

library(ggplot2)
library(reshape2)


args = commandArgs(trailingOnly=TRUE)
workDir = args[1]
filename = args[2]
outfile = args[3]

setwd(workDir)

##################
dat = read.table(filename, header = FALSE, row.names = 1)
colnames(dat) = c("Junko's", "overlap","mine")

pdf(outfile)
for(i in 1:nrow(dat)){
  d = dat[i,]
  d2 = melt(d)
  #
  p = ggplot(d2, aes(x = variable, y = value)) + geom_bar(stat = "identity", fill = "#AEAFAE", width=0.5) +
    theme_minimal() + ylab("number of DE genes") +xlab("") + labs(title=rownames(d))
  print(p)
}
dev.off()


