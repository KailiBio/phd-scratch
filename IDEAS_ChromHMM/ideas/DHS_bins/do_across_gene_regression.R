
# -- Kaili
# This script is for making regression r-square across cell type.

setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_bins/regression/")
library(reshape2)

exp = read.table("mm10_RNA_protein-coding_tpm_matrix_matched.txt", row.names = 1,
                 header=TRUE)

###################
# calculate gene exp sd, group gene by exp sd
###################
exp_sd = data.frame(apply(log10(exp+1e-5),1,sd))
colnames(exp_sd) = "sd"
exp0_name = rownames(exp_sd)[exp_sd$sd<=0.5]
exp0_5_name = rownames(exp_sd)[exp_sd$sd>0.5 & exp_sd$sd<=1]
exp1_name = rownames(exp_sd)[exp_sd$sd>1 & exp_sd$sd<=2]
exp2_name = rownames(exp_sd)[exp_sd$sd>2 & exp_sd$sd<=3]
exp3_name = rownames(exp_sd)[exp_sd$sd>3]

###################
# make figures
###################
## for (0,0.5]

## for (0.5,1]

## for (1,2]




## for (2,3]







###################