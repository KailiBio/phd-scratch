
# -- Kaili
# This scirpt is for making barplot for comparison genes overlappe with ubi-rDHS or active-DHS.

# setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/closest_gene/")
setwd("/data/zusers/fankaili/ccre/tf/closest_gene/gene_exp_comparison_file/")

library("ggplot2")

args <- commandArgs(trailingOnly=TRUE);
filename = args[1]

#filename = "A172_53_year"
data = read.table(paste(filename,".txt", sep=""))


ggplot(data, aes(x=V3, y=log10(V2+0.1), fill=V3)) + geom_bar(stat="identity") +
  labs(titl=filename, x="", y="expression (log10)")
ggsave(paste("/data/zusers/fankaili/ccre/tf/closest_gene/gene_exp_comparison_pdf/",filename, ".pdf", sep=""))
