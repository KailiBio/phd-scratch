#!/bin/bash

# -- Kaili
# This script is for calculating state transition rate.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/transition_rate/"

cd ${workDir}

######
normal_state_path="/data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/"
dhs_state_path="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_v3_100-400bp_result/"
for i in {1..19} X Y
do
    echo "chr"${i}
    # normal bins states
    awk '{FS=" ";OFS="\t"}{print $2,$3,$4,$5}' ${normal_state_path}run_IDEAS_8hm_atac_dname_pvalue.chr${i}.state | sort -k1,1 -k2,2n \
    > normal_bins_chr${i}_intestine_16.5_state.txt
    python ${scriptDir}calculate_transition_probability.py ${workDir}normal_bins_chr${i}_intestine_16.5_state.txt \
    "normal" ${workDir}intestine_16.5_chr${i}_normal_transitionCount.txt \
    ${workDir}intestine_16.5_chr${i}_normal_transitionProb.txt
    # dhs bins states
    awk '{FS=" ";OFS="\t"}{print $2,$3,$4,$5}' ${dhs_state_path}DHS_v3_100-400bp.chr${i}.state | sort -k1,1 -k2,2n \
    > dhs_bins_chr${i}_intestine_16.5_state.txt
    python ${scriptDir}calculate_transition_probability.py ${workDir}dhs_bins_chr${i}_intestine_16.5_state.txt \
    "dhs" ${workDir}intestine_16.5_chr${i}_dhs_transitionCount.txt \
    ${workDir}intestine_16.5_chr${i}_dhs_transitionProb.txt
done

## make figures
# Rscript make_state_transition_probability_figures.R
