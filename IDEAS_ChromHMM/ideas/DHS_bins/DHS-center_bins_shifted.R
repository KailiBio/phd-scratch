
# -- Kaili
# This script is for making figures to show how much do DHS-center bins shift.


setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/dhs_bins_shift/")
library(ggplot2)
library(grid)
library(gridExtra)

library(MASS)
library(viridis)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
####################
# barplot for all
####################
num=c(689275,112104,110723,34277)


df = data.frame(cbind(matrix(c(689275,112104,110723,34277),ncol=1),
                      type=c("1_no_sift","2_left_end","3_right_end","4_both_end")))
colnames(df) = c("num", "type")
df$num = as.numeric(as.character(df$num))

ggplot(df, aes(x= type,y=num, fill=type)) + 
  geom_bar(stat = "identity", width = 0.5) +
  theme_classic() +
  theme(axis.title.x = element_blank(), legend.position = "none",
        title = element_text(face="bold", size=12), 
        axis.text.x = element_text(face="bold",size=12)) +
  labs(title="number of DHS-center bins that shifted") +
  geom_text(aes(label = num), size=4, face="bold") +
  scale_fill_manual(values = c("#FFCD00", "#00B0F0", "#06DA93", "#FF0000"))

ggsave("basci_num_dhs_bins_shifted.pdf")
####################
# left shift
####################

left_dat = read.table("mm10_rOCR_bins_left_shift.txt")
colnames(left_dat) = c("id", "left", "right", "real_length", "total_length")

left_dat$density <- get_density(left_dat$real_length, left_dat$total_length, n=1000)

p1 = ggplot(left_dat, aes(x=left)) + geom_histogram(binwidth = 1, fill = "#00B0F0") +
  theme_classic() + xlab("nucluetide of left shifted (nt)") + 
  theme(title = element_text(face="bold", size=12))

p2 = ggplot(left_dat, aes(x = real_length, y = total_length, color = density)) + 
  geom_point() + scale_color_viridis() + theme_classic() +
  xlab("length of rOCRs") + ylab("length of bins") +
  theme(title = element_text(face="bold", size=12)) + 
  coord_cartesian(xlim=c(150,400), ylim=c(150,400))

p = grid.arrange(p1, p2, ncol=2)
ggsave("left_shifted.pdf",p, width = 10, height = 6)
####################
# right shift
####################
right_dat = read.table("mm10_rOCR_bins_right_shift.txt")
colnames(right_dat) = c("id", "left", "right", "real_length", "total_length")

right_dat$density <- get_density(right_dat$real_length, right_dat$total_length, n=1000)

p1 = ggplot(right_dat, aes(x=right)) + geom_histogram(binwidth = 1, fill = "#06DA93") +
  theme_classic() + xlab("nucluetide of right shifted (nt)") + 
  theme(title = element_text(face="bold", size=12))

p2 = ggplot(right_dat, aes(x = real_length, y = total_length, color = density)) + 
  geom_point() + scale_color_viridis() + theme_classic() +
  xlab("length of rOCRs") + ylab("length of bins") +
  theme(title = element_text(face="bold", size=12)) + 
  coord_cartesian(xlim=c(150,400), ylim=c(150,400))

p = grid.arrange(p1, p2, ncol=2)
ggsave("right_shifted.pdf",p, width = 10, height = 6)

####################
# both
####################
both_dat = read.table("mm10_rOCR_bins_both_shift.txt")
colnames(both_dat) = c("id", "left", "right", "real_length", "total_length")

both_dat$density <- get_density(both_dat$real_length, both_dat$total_length, n=1000)

p1 = ggplot(both_dat, aes(x=right-left)) + geom_histogram(binwidth = 1, fill = "#FF0000") +
  theme_classic() + xlab("nucluetide of right shifted (nt)") + 
  theme(title = element_text(face="bold", size=12))

p2 = ggplot(both_dat, aes(x = real_length, y = total_length, color = density)) + 
  geom_point() + scale_color_viridis() + theme_classic() +
  xlab("length of rOCRs") + ylab("length of bins") +
  theme(title = element_text(face="bold", size=12)) + 
  coord_cartesian(xlim=c(150,400), ylim=c(150,400))

p = grid.arrange(p1, p2, ncol=2)
ggsave("both_shifted.pdf",p, width = 10, height = 6)

####################
# box plot
####################
left = transform(left_dat$id, left_dat$left, type="left")
right = transform(right_dat$id, right_dat$right, type="right")
both = transform(both_dat$id, abs(both_dat$left-both_dat$right), type="both")
colnames(left) = colnames(right) = colnames(both) = c("id", "bp", "type")
matrix = rbind(left, right, both)
colnames(matrix) = c("id", "bp", "type")

ggplot(matrix, aes(x=type, y=bp, fill=type)) +
  geom_boxplot(width = 0.3) +
  theme_minimal() + xlab("") + ylab("shifted nucleotide") +
  theme(title = element_text(face="bold", size=12), legend.position = "none",
        axis.text.x = element_text(face="bold", size=12)) +
  scale_fill_manual(values = c("#00B0F0", "#06DA93", "#FF0000"))
  
ggsave("shifted_bp_boxplot.pdf")

ggplot(matrix, aes(x=type, y=bp, fill=type)) + geom_violin() +
  geom_boxplot(width = 0.3) +
  theme_minimal() + xlab("") + ylab("shifted nucleotide") +
  theme(title = element_text(face="bold", size=12), legend.position = "none",
        axis.text.x = element_text(face="bold", size=12)) +
  scale_fill_manual(values = c("#00B0F0", "#06DA93", "#FF0000"))

ggsave("shifted_bp_violin.pdf")

####################
# histogram
####################
no_dat = read.table("mm10_rOCR_bins_no_shift.txt")
colnames(no_dat) = c("id", "left", "right", "real_length", "total_length")
no = transform(no_dat$id, no_dat$left, type="no")
colnames(no) = c("id", "bp", "type")
matrix2 = rbind(matrix, no)

ggplot(matrix2, aes(x = bp)) + 
  geom_histogram(binwidth = 1, aes(y=..density..), fill="#AEAFAE") +
  theme_minimal() +
  xlab("shifted nucleotide") + theme(title = element_text(face="bold", size=12)) +
  labs(title="shifted nucleotides for all the DHS-center bins")

ggsave("histogram_density_shifted_nt.pdf")
####################