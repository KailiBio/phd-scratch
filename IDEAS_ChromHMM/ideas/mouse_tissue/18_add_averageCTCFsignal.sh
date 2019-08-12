#!/bin/bash

# -- Kaili
# This script is for adding average CTCF signal.
# 0. get average CTCF signal
# 1. add average-CTCF-signal as CTCF singnal for missing samples
# 2. add average-CTCF-signal as CTCF singnal for missing samples + CTCF motif signal
# 3. add average-CTCF-signal as a new mark
# 4. add average-CTCF-signal as a new mark + CTCF motif signal
# 5. 9to11+CTCFmotif+tCTCFaverage

workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
signalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"

cd ${workDir}

# 0. get average CTCF signal
awk '{FS=OFS="\t"}{print $4,0}' mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed > 9sample_sumCTCFsignal.txt
#
while read sample
do
    echo $sample
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $1,$2+a[$1]}}' ${signalDir}${sample}_CTCF_dhs.tab 9sample_sumCTCFsignal.txt > tmp.txt
    mv tmp.txt 9sample_sumCTCFsignal.txt
done < 9_ctcf_sample.txt
#
awk '{FS=OFS="\t"}{print $2/9}' 9sample_sumCTCFsignal.txt > 9sample_averageCTCFsignal.txt

# /data/zusers/fankaili/ideas/dhs_ctcf/9sample_averageCTCFsignal.txt

# 1. add average-CTCF-signal as CTCF singnal for missing samples
cp ctcf_9sample_to_11sample.sh 9to11_average.sh
vim 9to11_average.sh
cp ctcf_9sample_to_11sample.parafile 9to11_average.parafile
vim 9to11_average.parafile
cp ctcf_9sample_to_11sample.input 9to11_average.input
vim 9to11_average.input

nohup bash 9to11_average.sh > ./nohup/nohup.9to11_average.out 2>&1&
# z001 10155

# 2. add average-CTCF-signal as CTCF singnal for missing samples + CTCF motif signal
cp ctcf_9impute11_withMotif.sh 9to11_average_motif.sh
vim 9to11_average_motif.sh
cp ctcf_9impute11_withMotif.parafile 9to11_average_motif.parafile
vim 9to11_average_motif.parafile
cp ctcf_9impute11_withMotif.input 9to11_average_motif.input
echo "liver_14.5 CTCF /data/zusers/fankaili/ideas/dhs_ctcf/9sample_averageCTCFsignal.txt" >> 9to11_average_motif.input
echo "lung_14.5 CTCF /data/zusers/fankaili/ideas/dhs_ctcf/9sample_averageCTCFsignal.txt" >> 9to11_average_motif.input

nohup bash 9to11_average_motif.sh > ./nohup/nohup.9to11_average_motif.out 2>&1&
# z010 25381

# 3. add average-CTCF-signal as a new mark
## 1) make model
cp dhs_ctcf.sh 9sample_average.sh
vim 9sample_average.sh
cp dhs_ctcf.parafile 9sample_average.parafile
vim 9sample_average.parafile
cp dhs_ctcf.input 9sample_average.input
while read sample
do
    echo ${sample}
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/dhs_ctcf/9sample_averageCTCFsignal.txt" >> 9sample_average.input
done < 9_ctcf_sample.txt

nohup bash 9sample_average.sh > ./nohup/nohup.9sample_average.out 2>&1&
# z010 25629

## 2) imputation
awk '{FS=OFS}{if(NR==1){print $0}else{printf 100*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 100*$i};printf "\n"}}' \
/data/zusers/fankaili/ideas/dhs_ctcf/9sample_average_result/9sample_average.para0 > \
/data/zusers/fankaili/ideas/dhs_ctcf/9sample_average_result/9sample_average_100fold.para0
#
cp ctcf_9sample_impute_11sample.sh 9impute11_average_model.sh
vim 9impute11_average_model.sh
cp ctcf_9sample_impute_11sample.parafile 9impute11_average_model.parafile
vim 9impute11_average_model.parafile
cp ctcf_9sample_impute_11sample.input 9impute11_average_model.input
while read sample
do
    echo ${sample}
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/dhs_ctcf/9sample_averageCTCFsignal.txt" >> 9impute11_average_model.input
done < all_ctcf_sample.txt

nohup bash 9impute11_average_model.sh > ./nohup/nohup.9impute11_average_model.out 2>&1&
# z001 41157

# 4. add average-CTCF-signal as a new mark + CTCF motif signal
## 1) make model
cp ctcf_9sample_withMotif_model.sh 9sample_average_motif.sh
vim 9sample_average_motif.sh
cp ctcf_9sample_withMotif_model.parafile 9sample_average_motif.parafile
vim 9sample_average_motif.parafile
cp ctcf_9sample_withMotif_model.input 9sample_average_motif.input
while read sample
do
    echo ${sample}
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/dhs_ctcf/9sample_averageCTCFsignal.txt" >> 9sample_average_motif.input
done < 9_ctcf_sample.txt

nohup bash 9sample_average_motif.sh > ./nohup/nohup.9sample_average_motif.out 2>&1&
# z008 26388

## 2) imputation
awk '{FS=OFS}{if(NR==1){print $0}else{printf 100*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 100*$i};printf "\n"}}' \
/data/zusers/fankaili/ideas/dhs_ctcf/9sample_average_motif_result/9sample_average_motif.para0 > \
/data/zusers/fankaili/ideas/dhs_ctcf/9sample_average_motif_result/9sample_average_motif_100fold.para0
#
cp ctcf_9impute11_withMotif.sh 9impute11_average_motif_model.sh
vim 9impute11_average_motif_model.sh
cp ctcf_9impute11_withMotif.parafile 9impute11_average_motif_model.parafile
vim 9impute11_average_motif_model.parafile
cp ctcf_9impute11_withMotif.input 9impute11_average_motif_model.input
while read sample
do
    echo ${sample}
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/dhs_ctcf/9sample_averageCTCFsignal.txt" >> 9impute11_average_motif_model.input
done < all_ctcf_sample.txt

nohup bash 9impute11_average_motif_model.sh > ./nohup/nohup.9impute11_average_motif_model.out 2>&1&
# z001 41508


# 5. 9to11+CTCFmotif+tCTCFaverage
cp 9impute11_average_motif_model.sh 9to11_average_motif_model.sh
vim 9to11_average_motif_model.sh
cp 9impute11_average_motif_model.parafile 9to11_average_motif_model.parafile
vim 9to11_average_motif_model.parafile
cp 9impute11_average_motif_model.input 9to11_average_motif_model.input
vim 9to11_average_motif_model.input
# liver_14.5 CTCF /data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/liver_14.5_CTCF_dhs.txt
# lung_14.5 CTCF /data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/lung_14.5_CTCF_dhs.txt

nohup bash 9to11_average_motif_model.sh > ./nohup/nohup.9to11_average_motif_model.out 2>&1&
# z008 13515

prefix="9to11_average_motif_model"
bash ${scriptDir}make_each_sample_state_bed_tmp.sh /data/zusers/fankaili/ideas/dhs_ctcf/9to11_average_motif_model_result/9to11_average_motif_model.chr13.state /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_${prefix}/ 11
#
bash ${scriptDir}calculate_PRAU.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/9to11_average_motif_model_result/9to11_average_motif_model.para0 /data/zusers/fankaili/ideas/dhs_ctcf/
