
# -- Kaili
# This scirpt is for making barplot for comparison genes overlappe with ubi-rDHS or active-DHS.


library("ggplot2")

args <- commandArgs(trailingOnly=TRUE);
filename = args[1]
inputDir = args[2]
outputDir = args[3]

setwd(inputDir)

data = read.table(filename)
a = data[data$V3=="ubi-rOCR",]$V2
b = data[data$V3=="active-rOCR",]$V2
p = format(t.test(log10(a+0.1), log10(b+0.1))$p.value,scientific = TRUE, digits = 3)
t = paste(filename," (n=",as.character(length(a)),":",as.character(length(b)),
          ")\np=",as.character(p)," (t-test)",sep="")

ggplot(data, aes(x=V3, y=log10(V2+0.1), fill=V3)) + geom_boxplot() +
  labs(title=t, x="", y="TPM (log10)")
ggsave(paste(outputDir,filename, ".pdf", sep=""))
