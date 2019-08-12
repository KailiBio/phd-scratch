#!/bin/bash

# -- Kaili
# This script is pipeline for calculating AUPR and making PR AUC curve.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}


# 1. from IDEAS output
## sort states by averageCTCF signal
prefix=""
stateFile=""
num=""
parafile=""
## 1) step 1, preperation
bash ${scriptDir}make_each_sample_state_bed_tmp.sh ${stateFile} /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_${prefix}/ ${num}
## 2) step 2, get AUPR
bash ${scriptDir}calculate_PRAU.sh ${prefix} ${parafile} /data/zusers/fankaili/ideas/dhs_ctcf/
## 3) step 3, get PR AUC
for file in `ls ./state_bed_${prefix}/*_state_sorted.bed`
do
    file0=${file%_state_sorted.bed}
    sample=${file0#./state_bed_${prefix}/}
    echo $sample
    #
    Rscript ${scriptDir}make_PRcurve_peak.R /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/${prefix}/ ${sample}_${prefix}_peak_PR_pos.txt ${sample}_${prefix}_peak_PR_neg.txt ${sample}_${prefix}_PRcurve_peak
done





# 2. ChromImpute
# calculate_PRAU_other_imputation.sh

# 3. Avocado
