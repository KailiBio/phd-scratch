
# -- Kaili
# This script is for making histogram of numbers of high DNase samples for each rOCRs, then get the cutoff for ubi-rOCRs.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/hg38_rOCR/")

data = read.table("EDGE-DNase-Biosample-Counts.txt", row.names = 1)

pdf("hg38_rOCRs_high_DNase_histogram.pdf", width=8, height = 6)
hist(data[,1], breaks = max(data[,1]), col="grey", border = "grey", freq = F, main = "hg38 rOCRs",
     xlab="number of high DNase samples")
lines(c(580, 580), c(0,1), col="red", lwd=2, lty=2)
text(550,0.2,"580",col="red")
dev.off()

pdf("hg38_rOCRs_high_DNase_histogram2.pdf", width=8, height = 6)
hist(data[,1], breaks = max(data[,1]), col="grey", border = "grey", freq = F, ylim=c(0,0.005), 
     main = "hg38 rOCRs", xlab="number of high DNase samples")
lines(c(585, 585), c(0,1), col="red", lwd=2, lty=2)
text(550,0.004,"580",col="red")
dev.off()

####
data2=read.table("ccRE-DNase-Biosample-Counts.txt", row.names = 1)

pdf("hg19_rDHSs_high_DNase_histogram.pdf", width=8, height = 6)
hist(data2[,1], breaks = max(data2[,1]), col="grey", border = "grey", freq = F, main = "hg19 rDHSs",
     xlab="number of high DNase samples")
lines(c(450, 450), c(0,1), col="red", lwd=2, lty=2)
text(430,0.3,"450",col="red")
dev.off()

pdf("hg19_rDHSs_high_DNase_histogram2.pdf", width=8, height = 6)
hist(data2[,1], breaks = max(data2[,1]), col="grey", border = "grey", freq = F, main = "hg19 rDHSs",
     xlab="number of high DNase samples", ylim=c(0,0.005))
lines(c(450, 450), c(0,1), col="red", lwd=2, lty=2)
text(430,0.004,"450",col="red")
dev.off

# Apr25,2019
ggplot(data2, aes(x=V2)) + geom_histogram(aes(y=..density..),binwidth = 1, col="#8c8c8c") +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) +
  xlab("number of active biosamples") + ylab("density of regulatory elements") +
  geom_vline(xintercept = 450, color ="#E73A2F",linetype="dashed")
ggsave("hg19_rDHSs_high_DNase_histogram_new.png", width=5, height = 3.8)

ggplot(data2, aes(x=V2)) + geom_histogram(aes(y=..density..),binwidth = 1, col="#8c8c8c") +
  theme_minimal() + theme(title=element_text(face="bold",size=12)) +
  xlab("number of active biosamples") + ylab("density of regulatory elements") +
  geom_vline(xintercept = 450, color ="#E73A2F",linetype="dashed") +
  coord_cartesian(ylim=c(0,0.003))
ggsave("hg19_rDHSs_high_DNase_histogram_new2.png", width=6, height = 2)

### pie chart
install.packages("coolbutuseless/threed")
install.packages("coolbutuseless/ggthreed")
ggplot(df) + 
  geom_threedpie(aes(x = as.factor(value))) + 
  theme_void() + 
  theme(legend.position = 'bottom')

df <- data.frame(
  group = c("PLS", "ELS", "CTCF-only"),
  value = c(9009, 1826, 86)
)

library(googleVis)
op <- options(gvis.plot.tag = "chart")
pie <- gvisPieChart(df, options = list(title = "Sales per region", 
                                                 width = 1000, height = 500))
plot(pie)


library(plotrix)
group = c("PLS", "ELS", "CTCF-only")
value = c(9009, 1826, 86)
pdf("ss.pdf")
pie3D(value,labels = group,main="3D Pie chart of Sales per region", explode = 0)
dev.off()


region<-c("US","Europe","Japan","China","Others")
sales<-c(25000,12000 ,10000,5000,2000)
region_sales<-data.frame(region,sales)
pie3D(sales,labels = region,main="3D Pie chart of Sales per region")
