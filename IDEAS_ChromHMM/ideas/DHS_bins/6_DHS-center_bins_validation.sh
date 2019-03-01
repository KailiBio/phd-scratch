#!/bin/bash

# -- Kaili
# This script is for validating DHS-bins results.
# 1. aggregation plot
# 2. states transition rate
# 3. new aggregation plot (with corresponding region)
# 4. state combination for close OCR bins (<100bp)
# 5. gene expression prediction
# 6. conservation scores
# 7. state conversation between samples
# 8. validated enhancer prediction

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/"

cd ${workDir}

# 1. aggregation plot
bash ${scriptDir}make_DHSbins_marks_aggregation_plot.sh

# 2. states transition rate
bash ${scriptDir}calculate_state_transition_rate.sh

# 3. new aggregation plot (with corresponding region)
bash ${scriptDir}make_DHSbins_marks_aggregation_plot_corresponding.sh
# nohup bash /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/make_DHSbins_marks_aggregation_plot_corresponding.sh > ./nohup.match_aggregation_plot.out 2>&1&

# 4. state combination for close OCR bins (<100bp)
bash ${scriptDir}check_state_combination.sh

# 5. gene expression prediction
bash ${scriptDir}gene_expression_prediction.sh
bash ${scriptDir}gene_expression_prediction_across_gene.sh

# 6. conservation scores
bash ${scriptDir}make_conservation_score_aggregation.sh

# 7. state conversation between samples
### dhs bins
head -1 ./DHS_v3_100-400bp_result/DHS_v3_100-400bp.chr1.state > DHS_v3_100-400bp.state
for i in {1..19} X Y
do
    echo ${i}
    awk '{if(NR>1){print $0}}' ./DHS_v3_100-400bp_result/DHS_v3_100-400bp.chr${i}.state >> DHS_v3_100-400bp.state
done
#
# state_conservation_dhs_bins.sh
# z018 8288336

### normal bins
head -1 /data/zusers/fankaili/ideas/dhs_bins/normal_bins/66samples_10marks_normal_bins_result/66samples_10marks_normal_bins.chr1.state > normal_200bp.state
for i in {1..19} X Y
do
    echo ${i}
    awk '{if(NR>1){print $0}}' /data/zusers/fankaili/ideas/dhs_bins/normal_bins/66samples_10marks_normal_bins_result/66samples_10marks_normal_bins.chr${i}.state >> normal_200bp.state
done
#
# state_conservation_normal_bins.sh
# z018 8288839


# 8. validated enhancer prediction
