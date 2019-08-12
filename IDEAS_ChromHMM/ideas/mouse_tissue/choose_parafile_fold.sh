#!/bin/bash

# -- Kaili
# This script is for choosing parafile fold for imputation.
# 1. 1 fold
# 2. 10k
# 3. 1m
# 4. calculate and compare PRAU. (1-fold, 100-fold, 10k, and 1m)

workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
signalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"

cd ${workDir}

# 1. 1 fold
cp 9impute11_average_motif_model.sh 9impute11_average_motif_model_1fold.sh
vim 9impute11_average_motif_model_1fold.sh
cp 9impute11_average_motif_model.parafile 9impute11_average_motif_model_1fold.parafile
vim 9impute11_average_motif_model_1fold.parafile
cp 9impute11_average_motif_model.input 9impute11_average_motif_model_1fold.input

nohup bash 9impute11_average_motif_model_1fold.sh > ./nohup/nohup.9impute11_average_motif_model_1fold.out 2>&1&
# z010 4243

# 2. 10k
awk '{FS=OFS}{if(NR==1){print $0}else{printf 10000*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 10000*$i};printf "\n"}}' \
/data/zusers/fankaili/ideas/dhs_ctcf/9sample_average_motif_result/9sample_average_motif.para0 > \
/data/zusers/fankaili/ideas/dhs_ctcf/9sample_average_motif_result/9sample_average_motif_10k.para0
#
cp 9impute11_average_motif_model.sh 9impute11_average_motif_model_10k.sh
vim 9impute11_average_motif_model_10k.sh
cp 9impute11_average_motif_model.parafile 9impute11_average_motif_model_10k.parafile
vim 9impute11_average_motif_model_10k.parafile
cp 9impute11_average_motif_model.input 9impute11_average_motif_model_10k.input

nohup bash 9impute11_average_motif_model_10k.sh > ./nohup/nohup.9impute11_average_motif_model_10k.out 2>&1&
# z010 4402


# 3. 1m
awk '{FS=OFS}{if(NR==1){print $0}else{printf 1000000*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 1000000*$i};printf "\n"}}' \
/data/zusers/fankaili/ideas/dhs_ctcf/9sample_average_motif_result/9sample_average_motif.para0 > \
/data/zusers/fankaili/ideas/dhs_ctcf/9sample_average_motif_result/9sample_average_motif_1m.para0
#
cp 9impute11_average_motif_model.sh 9impute11_average_motif_model_1m.sh
vim 9impute11_average_motif_model_1m.sh
cp 9impute11_average_motif_model.parafile 9impute11_average_motif_model_1m.parafile
vim 9impute11_average_motif_model_1m.parafile
cp 9impute11_average_motif_model.input 9impute11_average_motif_model_1m.input

nohup bash 9impute11_average_motif_model_1m.sh > ./nohup/nohup.9impute11_average_motif_model_1m.out 2>&1&
# z008 10419



###########
# 4. calculate and compare PRAU.
prefix="9impute11_average_motif_model"
mkdir state_bed_${prefix}
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${workDir}${prefix}_result/ ${prefix}. ${workDir}state_bed_${prefix}/ 11
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${workDir}${prefix}_result/${prefix}.para0 ${workDir}
#
prefix="9impute11_average_motif_model_1fold"
mkdir state_bed_${prefix}
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${workDir}${prefix}_result/ ${prefix}. ${workDir}state_bed_${prefix}/ 11
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${workDir}${prefix}_result/${prefix}.para0 ${workDir}
#
prefix="9impute11_average_motif_model_10k"
mkdir state_bed_${prefix}
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${workDir}${prefix}_result/ ${prefix}. ${workDir}state_bed_${prefix}/ 11
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${workDir}${prefix}_result/${prefix}.para0 ${workDir}
#
prefix="9impute11_average_motif_model_1m"
mkdir state_bed_${prefix}
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${workDir}${prefix}_result/ ${prefix}. ${workDir}state_bed_${prefix}/ 11
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${workDir}${prefix}_result/${prefix}.para0 ${workDir}


# 5. 0.1 fold
awk '{FS=OFS}{if(NR==1){print $0}else{printf 0.1*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 0.1*$i};printf "\n"}}' \
/data/zusers/fankaili/ideas/dhs_ctcf/9sample_average_motif_result/9sample_average_motif.para0 > \
/data/zusers/fankaili/ideas/dhs_ctcf/9sample_average_motif_result/9sample_average_motif_0.1fold.para0
#
cp 9impute11_average_motif_model.sh 9impute11_average_motif_model_0.1fold.sh
vim 9impute11_average_motif_model_0.1fold.sh
cp 9impute11_average_motif_model.parafile 9impute11_average_motif_model_0.1fold.parafile
vim 9impute11_average_motif_model_0.1fold.parafile
cp 9impute11_average_motif_model.input 9impute11_average_motif_model_0.1fold.input

nohup bash 9impute11_average_motif_model_0.1fold.sh > ./nohup/nohup.9impute11_average_motif_model_0.1fold.out 2>&1&
# z010 57152

## 6. no imputation
cp 9impute11_average_motif_model.sh 11sample_average_motif_model.sh
vim 11sample_average_motif_model.sh
cp 9impute11_average_motif_model.parafile 11sample_average_motif_model.parafile
vim 11sample_average_motif_model.parafile
cp 9impute11_average_motif_model.input 11sample_average_motif_model.input
#
nohup bash 11sample_average_motif_model.sh > ./nohup/nohup.11sample_average_motif_model.out 2>&1&
# z003 7811
####
prefix="11sample_average_motif_model"
mkdir state_bed_${prefix}
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${workDir}${prefix}_result/ ${prefix}. ${workDir}state_bed_${prefix}/ 11
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${workDir}${prefix}_result/${prefix}.para0 ${workDir}
