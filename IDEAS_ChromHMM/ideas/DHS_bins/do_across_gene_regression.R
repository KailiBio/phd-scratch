
# -- Kaili
# This script is for making regression r-square across cell type.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/regression/")
library(reshape2)
library(ggplot2)

exp = read.table("mm10_RNA_protein-coding_tpm_matrix_matched.txt", row.names = 1,
                 header=TRUE)
###################
# 1. calculate gene exp sd, group gene by exp sd
###################
exp_sd = data.frame(apply(log10(exp+1e-5),1,sd))
colnames(exp_sd) = "sd"
exp_name_1 = rownames(exp_sd)[exp_sd$sd<=0.5]
exp_name_2 = rownames(exp_sd)[exp_sd$sd>0.5 & exp_sd$sd<=1]
exp_name_3 = rownames(exp_sd)[exp_sd$sd>1 & exp_sd$sd<=2]
exp_name_4 = rownames(exp_sd)[exp_sd$sd>2]

write.table(exp_name_1, "gene_list_sd_0-0.5.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(exp_name_2, "gene_list_sd_0.5-1.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(exp_name_3, "gene_list_sd_1-2.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(exp_name_4, "gene_list_sd_2.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
###################
# make figures
###################
dhs = read.table("gene_rsquare_matrix_dhs.txt", row.names = 1)
normal = read.table("gene_rsquare_matrix_normal.txt", row.names = 1)
colnames(dhs) = colnames(normal) = 1:20

names = intersect(rownames(dhs), rownames(normal))

## for (0,0.5]
a = intersect(exp_name_1, names)
dat_dhs = apply(dhs[a,],2,mean, na.rm=TRUE)
dat_200bp = apply(normal[a,],2,mean, na.rm=TRUE)
dat_dhs1 = transform(cbind(1:20, dat_dhs), type = "dhs")
dat_200bp1 = transform(cbind(1:20, dat_200bp), type = "200bp")
colnames(dat_dhs1) = colnames(dat_200bp1) = c("x", "r2", "type")
dat1 = data.frame(rbind(dat_dhs1, dat_200bp1))
p1 = ggplot(dat1, aes(x=x, y=r2, group=type, col = type)) + geom_line(size=1.2) +
  theme_minimal() + ylab("r-square") + xlab("0<Ø≤0.5 (n=13,329)") +
  theme(title = element_text(face="bold", size=12), legend.position = NaN) +
  coord_cartesian(ylim=c(0,0.4))
p1

## for (0.5,1]
a = intersect(exp_name_2, names)
dat_dhs = apply(dhs[a,],2,mean, na.rm=TRUE)
dat_200bp = apply(normal[a,],2,mean, na.rm=TRUE)
dat_dhs1 = transform(cbind(1:20, dat_dhs), type = "dhs")
dat_200bp1 = transform(cbind(1:20, dat_200bp), type = "200bp")
colnames(dat_dhs1) = colnames(dat_200bp1) = c("x", "r2", "type")
dat2 = data.frame(rbind(dat_dhs1, dat_200bp1))
p2 = ggplot(dat2, aes(x=x, y=r2, group=type, col = type)) + geom_line(size=1.2) +
  theme_minimal() + ylab("r-square") + xlab("0.5<Ø≤1 (n=3,335)") +
  theme(title = element_text(face="bold", size=12), legend.position = NaN)+
  coord_cartesian(ylim=c(0,0.4))
p2

## for (1,2]
a = intersect(exp_name_3, names)
dat_dhs = apply(dhs[a,],2,mean, na.rm=TRUE)
dat_200bp = apply(normal[a,],2,mean, na.rm=TRUE)
dat_dhs1 = transform(cbind(1:20, dat_dhs), type = "dhs")
dat_200bp1 = transform(cbind(1:20, dat_200bp), type = "200bp")
colnames(dat_dhs1) = colnames(dat_200bp1) = c("x", "r2", "type")
dat3 = data.frame(rbind(dat_dhs1, dat_200bp1))
p3 = ggplot(dat3, aes(x=x, y=r2, group=type, col = type)) + geom_line(size=1.2) +
  theme_minimal() + ylab("r-square") + xlab("1<Ø≤2 (n=4,054)") +
  theme(title = element_text(face="bold", size=12), legend.position = NaN)+
  coord_cartesian(ylim=c(0,0.4))
p3


## for (2,3]
a = intersect(exp_name_4, names)
dat_dhs = apply(dhs[a,],2,mean, na.rm=TRUE)
dat_200bp = apply(normal[a,],2,mean, na.rm=TRUE)
dat_dhs1 = transform(cbind(1:20, dat_dhs), type = "dhs")
dat_200bp1 = transform(cbind(1:20, dat_200bp), type = "200bp")
colnames(dat_dhs1) = colnames(dat_200bp1) = c("x", "r2", "type")
dat4 = data.frame(rbind(dat_dhs1, dat_200bp1))
p4 = ggplot(dat4, aes(x=x, y=r2, group=type, col = type)) + geom_line(size=1.2) +
  theme_minimal() + ylab("r-square") + xlab("Ø>2 (n=1,297)") +
  theme(title = element_text(face="bold", size=12), legend.position = c(0.7, 0.8)) +
  coord_cartesian(ylim=c(0,0.4))
p4


p = grid.arrange(p1, p2, p3, p4, ncol = 4)
ggsave("expression_regression_across_celltype.pdf", p, width = 9, height = 5)
ggsave("expression_regression_across_celltype.png", p, width = 9, height = 5)


###################
a = intersect(rownames(exp_sd), names)
dat_dhs = apply(dhs[a,],2,mean, na.rm=TRUE)
dat_200bp = apply(normal[a,],2,mean, na.rm=TRUE)
dat_dhs1 = transform(cbind(1:20, dat_dhs), type = "dhs")
dat_200bp1 = transform(cbind(1:20, dat_200bp), type = "200bp")
colnames(dat_dhs1) = colnames(dat_200bp1) = c("x", "r2", "type")
dat1 = data.frame(rbind(dat_dhs1, dat_200bp1))
ggplot(dat1, aes(x=x, y=r2, group=type, col = type)) + geom_line(size=1.2) +
  theme_minimal() + ylab("r-square") + xlab("") +
  theme(title = element_text(face="bold", size=12), legend.position = NaN) +
  coord_cartesian(ylim=c(0,0.4))
ggsave("expression_regression_across_celltype_all.pdf")
ggsave("expression_regression_across_celltype_all.png")
###################