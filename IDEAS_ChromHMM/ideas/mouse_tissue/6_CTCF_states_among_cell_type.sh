#!/bin/bash

# -- Kaili
# This script is for analyzing the cell-type specific of CTCF states.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

# 1. get CTCF states for all celltypes
mkdir ctcf_states
#
for j in {5..15}
do
    echo ${j}
    celltype=`head -1 /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.chr1.state | awk -v i="${j}" '{FS=" "}{print $i}'`
    #
    if [ -f ./ctcf_states/mm10_ctcf_state_${celltype}.bed ];then
        rm ./ctcf_states/mm10_ctcf_state_${celltype}.bed
    fi
    #
    for i in {1..19} X Y M
    do
        awk -v j="${j}" '{FS=" ";OFS="\t"}{s=0;if($j==28 ||$j==40 || $j==37 || $j==24 || $j==41 || $j==21 || $j==10 || $j==22 || $j==36){s=1};if(s){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' \
        /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.chr${i}.state >> ./ctcf_states/mm10_ctcf_state_${celltype}.bed
    done
    cut -f 1-4 ./ctcf_states/mm10_ctcf_state_${celltype}.bed | sort -u | sort -k1,1 -k2,2n > ./ctcf_states/mm10_ctcf_state_${celltype}_sorted.bed
done

# 2. get CTCF states overlapped matrix.
echo -e "lung_14.5\tliver_14.5\tstomach_0\tmidbrain_0\tkidney_0\tliver_0\tintestine_0\tlung_0\theart_0\thindbrain_0\tforebrain_0" \
> ctcf_states_overlapped.txt
echo -e "lung_14.5\tliver_14.5\tstomach_0\tmidbrain_0\tkidney_0\tliver_0\tintestine_0\tlung_0\theart_0\thindbrain_0\tforebrain_0" \
> ctcf_states_union.txt
#
for type1 in lung_14.5 liver_14.5 stomach_0 midbrain_0 kidney_0 liver_0 intestine_0 lung_0 heart_0 hindbrain_0 forebrain_0
do
    for type2 in lung_14.5 liver_14.5 stomach_0 midbrain_0 kidney_0 liver_0 intestine_0 lung_0 heart_0 hindbrain_0 forebrain_0
    do
        num=`cat ./ctcf_states/mm10_ctcf_state_${type1}_sorted.bed ./ctcf_states/mm10_ctcf_state_${type2}_sorted.bed | cut -f 4 | sort | uniq -d | wc -l`
        echo -n -e $num"\t" >> ctcf_states_overlapped.txt
        union=`cat ./ctcf_states/mm10_ctcf_state_${type1}_sorted.bed ./ctcf_states/mm10_ctcf_state_${type2}_sorted.bed | cut -f 4 | sort -u | wc -l`
        echo -n -e $union"\t" >> ctcf_states_union.txt
    done
    echo -n -e "\n" >> ctcf_states_overlapped.txt
    echo -n -e "\n" >> ctcf_states_union.txt
done

# Rscript figs_statistic_CTCF.R







####################
# redo this heatmap using new IDEAS run (all rep1)
# 1. get CTCF states for all celltypes
mkdir ctcf_states_rep1
#
for j in {5..15}
do
    echo ${j}
    celltype=`head -1 /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.chr1.state | awk -v i="${j}" '{FS=" "}{print $i}'`
    #
    if [ -f ./ctcf_states_rep1/mm10_ctcf_state_${celltype}.bed ];then
        rm ./ctcf_states_rep1/mm10_ctcf_state_${celltype}.bed
    fi
    #
    for i in {1..19} X Y
    do
        awk -v j="${j}" '{FS=" ";OFS="\t"}{s=0;if($j==21 ||$j==29 || $j==30 || $j==40 || $j==38 || $j==42 || $j==27 || $j==11 || $j==28 || $j==26 || $j==31){s=1};if(s){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' \
        /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.chr${i}.state >> ./ctcf_states_rep1/mm10_ctcf_state_${celltype}.bed
    done
    cut -f 1-4 ./ctcf_states_rep1/mm10_ctcf_state_${celltype}.bed | sort -u | sort -k1,1 -k2,2n > ./ctcf_states_rep1/mm10_ctcf_state_${celltype}_sorted.bed
done

# 2. get CTCF states overlapped matrix.
echo -e "lung_14.5\tliver_14.5\tstomach_0\tmidbrain_0\tkidney_0\tliver_0\tintestine_0\tlung_0\theart_0\thindbrain_0\tforebrain_0" \
> ./ctcf_states_rep1/ctcf_states_overlapped.txt
echo -e "lung_14.5\tliver_14.5\tstomach_0\tmidbrain_0\tkidney_0\tliver_0\tintestine_0\tlung_0\theart_0\thindbrain_0\tforebrain_0" \
> ./ctcf_states_rep1/ctcf_states_union.txt
#
for type1 in lung_14.5 liver_14.5 stomach_0 midbrain_0 kidney_0 liver_0 intestine_0 lung_0 heart_0 hindbrain_0 forebrain_0
do
    for type2 in lung_14.5 liver_14.5 stomach_0 midbrain_0 kidney_0 liver_0 intestine_0 lung_0 heart_0 hindbrain_0 forebrain_0
    do
        num=`cat ./ctcf_states_rep1/mm10_ctcf_state_${type1}_sorted.bed ./ctcf_states_rep1/mm10_ctcf_state_${type2}_sorted.bed | cut -f 4 | sort | uniq -d | wc -l`
        echo -n -e $num"\t" >> ./ctcf_states_rep1/ctcf_states_overlapped.txt
        union=`cat ./ctcf_states_rep1/mm10_ctcf_state_${type1}_sorted.bed ./ctcf_states_rep1/mm10_ctcf_state_${type2}_sorted.bed | cut -f 4 | sort -u | wc -l`
        echo -n -e $union"\t" >> ./ctcf_states_rep1/ctcf_states_union.txt
    done
    echo -n -e "\n" >> ./ctcf_states_rep1/ctcf_states_overlapped.txt
    echo -n -e "\n" >> ./ctcf_states_rep1/ctcf_states_union.txt
done

# Rscript make_11states_jaccard_heatmap.R
