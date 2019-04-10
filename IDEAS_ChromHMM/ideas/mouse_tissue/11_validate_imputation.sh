#!/bin/bash

# -- Kaili
# This script is for validating imputated states.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

# 1. state conservation
## 1) for all the runs, get state file together
get_state_together(){
    prefix=$1
    #
    head -1 ./${prefix}_result/${prefix}.chr1.state > ${prefix}.state
    for i in {1..19} X Y
    do
        awk '{if(NR>1){print $0}}' ./${prefix}_result/${prefix}.chr${i}.state >> ${prefix}.state
    done
}
#
get_state_together ctcf_8sample_impute_11sample
get_state_together ctcf_8sample_impute_11sample_2
get_state_together ctcf_8sample_to_11sample
get_state_together ctcf_8sample_to_11sample_2
#
get_state_together ctcf_9sample_impute_11sample
get_state_together ctcf_9sample_impute_11sample_2
get_state_together ctcf_9sample_to_11sample
get_state_together ctcf_9sample_to_11sample_2

## 2)
# 8 samples model
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/dhs_ctcf/ ctcf_8sample_impute_11sample.state ctcf_8sample_impute_11sample_2.state 11 47 "8 samples impute 11 samples (reproducibility)" "ARI_8sample_imputation_reproducibility"
#
Rscript  ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/dhs_ctcf/ ctcf_8sample_to_11sample.state ctcf_8sample_to_11sample_2.state 11 47 "8 samples to 11 samples (reproducibility)" "ARI_8sample_withdata_reproducibility"
#
Rscript  ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/dhs_ctcf/ ctcf_8sample_impute_11sample.state ctcf_8sample_to_11sample.state 11 47 "8 samples impute 11 samples vs. to 11 samples" "ARI_8sample_imputation_validation1"
#
Rscript  ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/dhs_ctcf/ ctcf_8sample_impute_11sample_2.state ctcf_8sample_to_11sample_2.state 11 47 "8 samples impute 11 samples vs. to 11 samples" "ARI_8sample_imputation_validation2"
####
# 9 samples model
Rscript  ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/dhs_ctcf/ ctcf_9sample_impute_11sample.state ctcf_9sample_impute_11sample_2.state 11 47 "9 samples impute 11 samples (reproducibility)" "ARI_9sample_imputation_reproducibility"
#
Rscript  ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/dhs_ctcf/ ctcf_9sample_to_11sample.state ctcf_9sample_to_11sample_2.state 11 47 "9 samples to 11 samples (reproducibility)" "ARI_9sample_withdata_reproducibility"
#
Rscript  ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/dhs_ctcf/ ctcf_9sample_impute_11sample.state ctcf_9sample_to_11sample.state 11 47 "9 samples impute 11 samples vs. to 11 samples" "ARI_9sample_imputation_validation1"
#
Rscript  ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/dhs_ctcf/ ctcf_9sample_impute_11sample_2.state ctcf_9sample_to_11sample_2.state 11 47 "9 samples impute 11 samples vs. to 11 samples" "ARI_9sample_imputation_validation2"


# Rscript make_boxplot_ARImatrix.R
