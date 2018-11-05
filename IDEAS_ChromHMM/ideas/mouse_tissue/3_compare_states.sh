#!/bin/bash

# -- Kaili
# This script is for comparing states result form different IDEAS run (with or without CTCF).
# 1. get mean mark signal for each state. (matrix)
# 2. correlation between states in two runs.

#--------------------------------------------------------------------
# 1. get mean mark signal for each state. (matrix)

## 1) 66samples_10marks_pvalue

# make slurm script
# python /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/get_state_mean_signal_chr.py \
# /data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/ \
# run_IDEAS_8hm_atac_dname_pvalue.chr /data/zusers/fankaili/ideas/run_ideas_p_value/signal/ \
# /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input \
# /data/zusers/fankaili/ideas/compare_results/66samples_10marks_pvalue_state_signal_sum_chr${SLURM_ARRAY_TASK_ID}.txt \
# /data/zusers/fankaili/ideas/compare_results/66samples_10marks_pvalue_state_signal_num_chr${SLURM_ARRAY_TASK_ID}.txt \
# ${SLURM_ARRAY_TASK_ID}

# python /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/get_state_signal_mean_chr_mark_sample.py \
# /data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/run_IDEAS_8hm_atac_dname_pvalue.chr \
# /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input \
# /data/zusers/fankaili/ideas/compare_results/66samples_10marks/66samples_10marks_state_signal_sum_H3K4me2_ \
# /data/zusers/fankaili/ideas/compare_results/66samples_10marks/66samples_10marks_state_signal_num_H3K4me2_ \
# H3K4me2 ${SLURM_ARRAY_TASK_ID}

# cd /home/fankaili/IDEAS_results/
# sbatch 66samples_10marks_H3K4me2.sh
###########

python /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/get_state_signal_mean_mark.py \
/data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/run_IDEAS_8hm_atac_dname_pvalue.chr \
/data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input \
/data/zusers/fankaili/ideas/compare_results/66samples_10marks/66samples_10marks_state_signal_mean_${mark}.txt ${mark} 66

cd /home/fankaili/IDEAS_results/
sbatch 66samples_10marks_signal_mean.sh

## 2) e14.5p0_21biosamples_11marks

# python /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/get_state_signal_mean_chr_mark_sample.py \
# /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1_result/e14.5p0_8hm_ATAC_DNAme_CTCF_1.chr \
# /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input \
# /data/zusers/fankaili/ideas/compare_results/e14.5p0_21biosamples_11marks/e14.5p0_21biosamples_11marks_state_signal_sum_H3K4me2_ \
# /data/zusers/fankaili/ideas/compare_results/e14.5p0_21biosamples_11marks/e14.5p0_21biosamples_11marks_signal_num_H3K4me2_ \
# H3K4me2 ${SLURM_ARRAY_TASK_ID}
#
#
# cd /home/fankaili/IDEAS_results/
# sbatch e14.5p0_21biosamples_11marks_H3K4me2.sh

python /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/get_state_signal_mean_mark.py \
/data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1_result/e14.5p0_8hm_ATAC_DNAme_CTCF_1.chr \
/data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input \
/data/zusers/fankaili/ideas/compare_results/e14.5p0_21biosamples_11marks/e14.5p0_21biosamples_11marks_state_signal_mean_${mark}.txt \
${mark} 21

cd /home/fankaili/IDEAS_results/
sbatch e14.5p0_21biosamples_11marks_signal_mean.sh

## 3) e14.5p0_21biosamples_10marks

# python /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/get_state_signal_mean_chr_mark_sample.py \
# /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/e14.5p0_8hm_ATAC_DNAme_result/e14.5p0_8hm_ATAC_DNAme.chr \
# /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/e14.5p0_8hm_ATAC_DNAme.input \
# /data/zusers/fankaili/ideas/compare_results/e14.5p0_21biosamples_10marks/e14.5p0_21biosamples_10marks_state_signal_sum_H3K4me2_ \
# /data/zusers/fankaili/ideas/compare_results/e14.5p0_21biosamples_10marks/e14.5p0_21biosamples_10marks_signal_num_H3K4me2_ \
# H3K4me2 ${SLURM_ARRAY_TASK_ID}
#
# cd /home/fankaili/IDEAS_results/
# sbatch e14.5p0_21biosamples_10marks_H3K4me2.sh

python /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/get_state_signal_mean_mark.py \
/data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/e14.5p0_8hm_ATAC_DNAme_result/e14.5p0_8hm_ATAC_DNAme.chr \
/data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/e14.5p0_8hm_ATAC_DNAme.input \
/data/zusers/fankaili/ideas/compare_results/e14.5p0_21biosamples_10marks/e14.5p0_21biosamples_10marks_state_signal_mean_${mark}.txt \
${mark} 21

cd /home/fankaili/IDEAS_results/
sbatch e14.5p0_21biosamples_10marks_signal_mean.sh


## 4) DHS_v1_100-300bp
python /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/get_state_signal_mean_mark.py \
/data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/v1_100_300bp/DHS_v1_100-300bp_result/DHS_v1_100-300bp.chr \
/data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/v1_100_300bp/DHS_v1_100-300bp.input \
/data/zusers/fankaili/ideas/compare_results/DHS_v1_100-300bp/DHS_v1_100-300bp_state_signal_mean_${mark}.txt \
${mark} 66

cd /home/fankaili/IDEAS_results/
sbatch DHS_v1_100-300bp_signal_mean.sh


## 5) DHS_v2_1-300bp
python /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/get_state_signal_mean_mark.py \
/data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/v2_1_300bp/DHS_v2_1-300bp_result/DHS_v2_1-300bp.chr \
/data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/v2_1_300bp/DHS_v2_1-300bp.input \
/data/zusers/fankaili/ideas/compare_results/DHS_v2_1-300bp/DHS_v2_1-300bp_state_signal_mean_${mark}.txt \
${mark} 66

cd /home/fankaili/IDEAS_results/
sbatch DHS_v2_1-300bp_signal_mean.sh



#--------------------------------------------------------------------
# 2. correlation between states in two runs.

# cd /data/zusers/fankaili/ideas/compare_results/e14.5p0_21biosamples_10marks/
# n=`wc -l e14.5p0_21biosamples_10marks_state_signal_mean_ATAC.txt| awk '{print $1-2}'`
# #
# echo "state" > e14.5p0_21biosamples_10marks_state_signal_mean_matrix.txt
# for state in $(seq 0 $n)
# do
#     echo $state >> e14.5p0_21biosamples_10marks_state_signal_mean_matrix.txt
# done
# #
# for mark in ATAC DNAme H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3 H3K27ac H3K27me3 H3K36me3
# do
#     echo ${mark}
#     awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$1]}}' e14.5p0_21biosamples_10marks_state_signal_mean_${mark}.txt \
#     e14.5p0_21biosamples_10marks_state_signal_mean_matrix.txt > temp.txt
#     mv temp.txt e14.5p0_21biosamples_10marks_state_signal_mean_matrix.txt
# done

# Rscript make_states_mean_signal_heatmap.R
# Rscript make_state_correlation.R

### run locally
## 1) between 21 samples, with/without CTCF
Rscript make_state_correlation.R "e14.5p0_8hm_ATAC_DNAme.para0" "e14.5p0_8hm_ATAC_DNAme_CTCF_1.para0" \
"21_biosamples_without_CTCF" "21_biosamples_with_CTCF" "21 biosamples in e14.5&p0, 8HM+ATAC+DNAme" "e14.5p0_correlation.pdf" "CTCF"

## 2) 66 samples without CTCF & 21 samples with CTCF
Rscript make_state_correlation.R "run_IDEAS_8hm_atac_dname_pvalue.para0" "e14.5p0_8hm_ATAC_DNAme_CTCF_1.para0" \
"66_biosamples_without_CTCF" "21_biosamples_with_CTCF" "states comparison(38 states vs. 37 states)" "66samples_21samplesCTCF_correlation.pdf" "CTCF"

## 3) spearman

## 4) hcluster
# Rscript compare_state_euclidean_Hierarchical.R
