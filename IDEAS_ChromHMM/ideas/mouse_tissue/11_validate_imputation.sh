#!/bin/bash

# -- Kaili
# This script is for validating imputated states.
# 1. state conservation
# 2. ratio of CTCF peaks that can be capture by CTCF states
# 3. predicted peaks for impuated samples
# 4. difference between 9-impute-11 and 9-impute-66
# 5. gene expression prediction
# 6. how good the imputation it is?


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


# 2. ratio of CTCF peaks that can be capture by CTCF states
bash ${scriptDir}compare_ratio_of_ctcf_peak_be_captured.sh


# 3. predicted peaks for impuated samples


# 4. difference between 9-impute-11 and 9-impute-66
# get_state_together function above
get_state_together ctcf_9sample_impute
awk '{FS=OFS}{print $1,$2,$3,$4,$52,$48,$67,$55,$35,$44,$31,$51,$17,$24,$10,$71}' ctcf_9sample_impute.state > ctcf_9sample_impute_11sampleOnly.state
#
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/dhs_ctcf/ ctcf_9sample_impute_11sample.state ctcf_9sample_impute_11sampleOnly.state 11 47 "9-impute-11 vs. 9-impute-66" "compare_9impute11_9impute66_ARI"

# 5. gene expression prediction
bash ${scriptDir}gene_expression_regression_addCTCF.sh

# 6. how good the imputation it is?
bash ${scriptDir}validate_majority_vote.sh
