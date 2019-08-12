#!/bin/bash

# -- Kaili
# This script is for validating imputation resulting using or comparing rep2.
# 1. ARI matrix
# 2. ARI bar for one sample
# 3. recall ratio

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

# 1. ARI matrix
get_state_together(){
    prefix=$1
    #
    head -1 ./${prefix}_result/${prefix}.chr1.state > ${prefix}.state
    for i in {1..19} X Y
    do
        awk '{if(NR>1){print $0}}' ./${prefix}_result/${prefix}.chr${i}.state >> ${prefix}.state
    done
}
## 1) rep1 vs. rep1-impute-rep2
get_state_together ctcf_samples_rep1-impute-rep2
get_state_together ctcf_samples
#
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ ctcf_samples.state ctcf_samples_rep1-impute-rep2.state 11 44 "rep1 vs. rep1-impute-rep2" "compare_rep1_rep1-impute-rep2_ARI"
## 2) rep1 vs. rep1-impute-rep1
get_state_together ctcf_samples_rep1-impute-rep1
#
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ ctcf_samples.state ctcf_samples_rep1-impute-rep1.state 11 42 "rep1 vs. rep1-impute-rep1" "compare_rep1_rep1-impute-rep1_ARI"
# 3) rep1 vs. rep1-impute-rep1 (given 9 samples CTCF)
get_state_together ctcf_samples_rep1-impute-rep1_given9samples
#
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ ctcf_samples.state ctcf_samples_rep1-impute-rep1_given9samples.state 11 42 "rep1 vs. rep1-impute-rep1\n(given 9 sample CTCF)" "compare_rep1_rep1-impute-rep1_given9samples_ARI"
# 4) rep1 vs. rep1-impute-rep2 (given 9 samples CTCF)
get_state_together ctcf_samples_rep1-impute-rep2_given9samples
#
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ ctcf_samples.state ctcf_samples_rep1-impute-rep2_given9samples.state 11 42 "rep1 vs. rep1-impute-rep2\n(given 9 sample CTCF)" "compare_rep1_rep1-impute-rep2_given9samples_ARI"
#
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ ctcf_samples_rep1-impute-rep1_given9samples.state ctcf_samples_rep1-impute-rep2_given9samples.state 11 42 "rep1-impute-rep1 vs. rep1-impute-rep2\n(given 9 sample CTCF)" "compare_rep1-impute-rep1_rep1-impute-rep2_given9samples_ARI"
# 5) rep1 vs. rep1-impute-rep2 (given 9 samples CTCF, 10k)
get_state_together ctcf_samples_rep1-impute-rep2_given9samples_10k
#
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ ctcf_samples.state ctcf_samples_rep1-impute-rep2_given9samples_10k.state 11 42 "rep1 vs. rep1-impute-rep2\n(given 9 sample CTCF, 10k)" "compare_rep1_rep1-impute-rep2_given9samples_10k_ARI"
# 6) rep1-impute-rep2 (given 9 samples CTCF) vs. rep1-impute-rep2 (given 9 samples CTCF, 10k)
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns.R /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ ctcf_samples_rep1-impute-rep2_given9samples.state ctcf_samples_rep1-impute-rep2_given9samples_10k.state 11 42 "rep1-impute-rep2(given 9 sample CTCF): 100fold vs. 10k" "compare_rep1-impute-rep2_given9samples_100_10k_ARI"


# 2. ARI bar for one sample
## 1) rep1-impute-rep1 lung_14.5
get_state_together ctcf_samples_rep1-impute-rep1_oneSample
awk '{FS=OFS}{print $1,$2,$3,$4,$5,$16}' ctcf_samples.state > ctcf_samples_lung14.5.state
#
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns_oneSample.R /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ ctcf_samples_lung14.5.state ctcf_samples_rep1-impute-rep1_oneSample.state 1 42 "rep1 vs. rep1-impute-rep1\n(lung_14.5)" "compare_rep1_rep1-impute-rep1_ARI_lung_14.5"
## 2) rep1-impute-rep1 lung_14.5
get_state_together ctcf_samples_rep1-impute-rep2_oneSample
#
Rscript ${scriptDir}get_ARI_heatmap_between_IDEASruns_oneSample.R /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ ctcf_samples_lung14.5.state ctcf_samples_rep1-impute-rep2_oneSample.state 1 42 "rep1 vs. rep2-impute-rep1\n(lung_14.5)" "compare_rep1_rep2-impute-rep1_ARI_lung_14.5"

# 3. recall ratio
bash ${scriptDir}compare_ratio_of_ctcf_peak_be_captured_rep2.sh
