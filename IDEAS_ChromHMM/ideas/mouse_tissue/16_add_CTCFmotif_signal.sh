#!/bin/bash

# -- Kaili
# This script is for adding CTCF motif as signal for CTCF imputation.

workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"

cd ${workDir}

# 1. make CTCF motif signal files
# column 6 in motif bed output is log-odds score.
intersectBed -a mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed -b /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/mm10_CTCF_motif_region.bed -wa -wb | cut -f 4,9 | sort -k1,1 -k2,2nr | awk 'BEGIN{FS=OFS="\t";bin=""}{if(NR==1){bin=$1;print $0}else{if(bin!=$1){bin=$1;print $0}}}' > CTCFmotif_bin_withScore.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$4]){print a[$4]}else{print 0}}}' CTCFmotif_bin_withScore.txt mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed > CTCFmotif_signal_OCR.txt

# 2. run 9sample model with CTCF motif
cp dhs_ctcf.sh ctcf_9sample_withMotif_model.sh
vim ctcf_9sample_withMotif_model.sh
cp dhs_ctcf.parafile ctcf_9sample_withMotif_model.parafile
vim ctcf_9sample_withMotif_model.parafile
cp dhs_ctcf.input ctcf_9sample_withMotif_model.input
while read sample
do
    echo ${sample}" motif /data/zusers/fankaili/ideas/dhs_ctcf/CTCFmotif_signal_OCR.txt" >> ctcf_9sample_withMotif_model.input
done < all_ctcf_sample.txt
vim ctcf_9sample_withMotif_model.input

nohup bash ctcf_9sample_withMotif_model.sh > ./nohup/nohup.ctcf_9sample_withMotif_model.out 2>&1&


# 3. imputation with CTCF motif
awk '{FS=OFS}{if(NR==1){print $0}else{printf 100*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 100*$i};printf "\n"}}' \
/data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_withMotif_model_result/ctcf_9sample_withMotif_model.para0 > \
/data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_withMotif_model_result/ctcf_9sample_withMotif_model_100fold.para0
#
cp ctcf_9sample_impute_11sample.sh ctcf_9impute11_withMotif.sh
vim ctcf_9impute11_withMotif.sh
cp ctcf_9sample_impute_11sample.parafile ctcf_9impute11_withMotif.parafile
vim ctcf_9impute11_withMotif.parafile
cp ctcf_9sample_impute_11sample.input ctcf_9impute11_withMotif.input
while read sample
do
    echo ${sample}" motif /data/zusers/fankaili/ideas/dhs_ctcf/CTCFmotif_signal_OCR.txt" >> ctcf_9impute11_withMotif.input
done < all_ctcf_sample.txt

nohup bash ctcf_9impute11_withMotif.sh > ./nohup/nohup.ctcf_9impute11_withMotif.out 2>&1&
# z001 65038
