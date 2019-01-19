
# -- Kaili
# This scirpt is for making barplot for comparison genes overlappe with ubi-rDHS or active-DHS.

# setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/closest_gene/")
setwd("/data/zusers/fankaili/ccre/tf/closest_gene/gene_exp_comparison_file/")

library("ggplot2")

args <- commandArgs(trailingOnly=TRUE);
filename = args[1]
#filename = "GM12878"
#filename = "A172_53_year"
data = read.table(paste(filename,".txt", sep=""))
a = data[data$V3=="ubi-rDHS_overlapped",]$V2
b = data[data$V3=="active-rDHS_overlapped",]$V2
p = format(t.test(log10(a+0.1), log10(b+0.1))$p.value,scientific = TRUE, digits = 3)
t = paste(filename," (n=",as.character(length(a)),":",as.character(length(b)),
          ")\np=",as.character(p)," (t-test)",sep="")

ggplot(data, aes(x=V3, y=log10(V2+0.1), fill=V3)) + geom_boxplot() +
  labs(title=t, x="", y="TPM (log10)")
ggsave(paste("/data/zusers/fankaili/ccre/tf/closest_gene/gene_exp_comparison_pdf/",filename, ".pdf", sep=""))

