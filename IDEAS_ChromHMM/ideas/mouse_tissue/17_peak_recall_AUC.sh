#!/bin/bash

# -- Kaili
# This script is for calculating AUC.
# 1. 9impute11
# 2. 9to11
# 3. 9sample majority vote

workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
signalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"

cd ${workDir}
mkdir CTCF_signal_based_curve
mkdir state_ranked_PR
mkdir state_ranked_PR2

# 0. calculate states with CTCF peak in each sample
while read sample
do
    echo ${sample}
    #
    intersectBed -a mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed -b ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -wa -u > ./state_ranked_PR/${sample}_CTCFpeak_pos.txt
    intersectBed -a mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed -b ./peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed -v > ./state_ranked_PR/${sample}_CTCFpeak_neg.txt
done < all_ctcf_sample.txt
#####
while read sample
do
    echo ${sample}
    #
    intersectBed -a mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed -b ./peaks_validation_9To11/${sample}_ctcf_peak_sorted.bed -wa -u > ./state_ranked_PR2/${sample}_CTCFpeak_pos.txt
    intersectBed -a mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed -b ./peaks_validation_9To11/${sample}_ctcf_peak_sorted.bed -v > ./state_ranked_PR2/${sample}_CTCFpeak_neg.txt
done < all_ctcf_sample.txt


# 1. 9impute11
prefix="9impute11"
bash ${scriptDir}make_CTCF_signal_based_curve.sh ${prefix}
#
bash ${scriptDir}make_PRcurve_peak_state-ranking.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_impute_11sample_result/ctcf_9sample_impute_11sample.para0
#
bash ${scriptDir}make_PRcurve_peak_state-ranking2.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_impute_11sample_result/ctcf_9sample_impute_11sample.para0


# 2. 9to11
prefix="9to11"
bash ${scriptDir}make_CTCF_signal_based_curve.sh ${prefix}
#
bash ${scriptDir}make_PRcurve_peak_state-ranking.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_to_11sample_result/ctcf_9sample_to_11sample.para0
#
bash ${scriptDir}make_PRcurve_peak_state-ranking2.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9sample_to_11sample_result/ctcf_9sample_to_11sample.para0


# 3. 9sample majority vote
prefix="majority-vote"
bash ${scriptDir}make_CTCFsignal_based_curve.sh ${prefix}
##
mkdir ./state_ranked_PR/${prefix}/
mkdir ./state_ranked_PR2/${prefix}/
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5}else{if(a[$4]){print $4,a[$4]}else{print $4,0}}}' ./state_bed_9sample/CTCFstates_count.txt mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed > ./state_ranked_PR/${prefix}/all_state_counted.txt
#
for file in `ls ./state_bed_9impute11/*_state_sorted.bed`
do
    file0=${file%_state_sorted.bed}
    sample=${file0#./state_bed_9impute11/}
    echo $sample
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $4,a[$4]}}' ./state_ranked_PR/${prefix}/all_state_counted.txt ./state_ranked_PR/${sample}_CTCFpeak_pos.txt > ./state_ranked_PR/${prefix}/${sample}_${prefix}_peak_PR_pos.txt
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $4,a[$4]}}' ./state_ranked_PR/${prefix}/all_state_counted.txt ./state_ranked_PR/${sample}_CTCFpeak_neg.txt > ./state_ranked_PR/${prefix}/${sample}_${prefix}_peak_PR_neg.txt
    #
    Rscript ${scriptDir}make_PRcurve_peak.R /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/${prefix}/ ${sample}_${prefix}_peak_PR_pos.txt ${sample}_${prefix}_peak_PR_neg.txt ${sample}_${prefix}_PRcurve_peak
done
#
for file in `ls ./state_bed_9impute11/*_state_sorted.bed`
do
    file0=${file%_state_sorted.bed}
    sample=${file0#./state_bed_9impute11/}
    echo $sample
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $4,a[$4]}}' ./state_ranked_PR/${prefix}/all_state_counted.txt ./state_ranked_PR2/${sample}_CTCFpeak_pos.txt > ./state_ranked_PR2/${prefix}/${sample}_${prefix}_peak_PR_pos.txt
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $4,a[$4]}}' ./state_ranked_PR/${prefix}/all_state_counted.txt ./state_ranked_PR2/${sample}_CTCFpeak_neg.txt > ./state_ranked_PR2/${prefix}/${sample}_${prefix}_peak_PR_neg.txt
    #
    Rscript ${scriptDir}make_PRcurve_peak.R /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR2/${prefix}/ ${sample}_${prefix}_peak_PR_pos.txt ${sample}_${prefix}_peak_PR_neg.txt ${sample}_${prefix}_PRcurve_peak_2
done


# 4. 9-impute-11 withMotif
mkdir state_bed_9impute11motif
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9impute11_withMotif_result/ ctcf_9impute11_withMotif. /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_9impute11motif/ 11
#
prefix="9impute11motif"
bash ${scriptDir}make_CTCF_signal_based_curve.sh ${prefix}
#
bash ${scriptDir}make_PRcurve_peak_state-ranking.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9impute11_withMotif_result/ctcf_9impute11_withMotif.para0
#
bash ${scriptDir}make_PRcurve_peak_state-ranking2.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/ctcf_9impute11_withMotif_result/ctcf_9impute11_withMotif.para0


# 5. 9to11_average
mkdir state_bed_9to11_average
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/dhs_ctcf/9to11_average_result/ 9to11_average. /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_9to11_average/ 11
#
prefix="9to11_average"
bash ${scriptDir}make_CTCF_signal_based_curve.sh ${prefix}
#
bash ${scriptDir}make_PRcurve_peak_state-ranking.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/9to11_average_result/9to11_average.para0
#
bash ${scriptDir}make_PRcurve_peak_state-ranking2.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/9to11_average_result/9to11_average.para0


# 6. 9to11_average_motif
mkdir state_bed_9to11_average_motif
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/dhs_ctcf/9to11_average_motif_result/ 9to11_average_motif. /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_9to11_average_motif/ 11
#
prefix="9to11_average_motif"
bash ${scriptDir}make_CTCF_signal_based_curve.sh ${prefix}
#
bash ${scriptDir}make_PRcurve_peak_state-ranking.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/9to11_average_motif_result/9to11_average_motif.para0
#
bash ${scriptDir}make_PRcurve_peak_state-ranking2.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/9to11_average_motif_result/9to11_average_motif.para0


# 7. 9impute11_average_model
mkdir state_bed_9impute11_average_model
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/dhs_ctcf/9impute11_average_model_result/ 9impute11_average_model. /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_9impute11_average_model/ 11
#
prefix="9impute11_average_model"
bash ${scriptDir}make_CTCF_signal_based_curve.sh ${prefix}
#
bash ${scriptDir}make_PRcurve_peak_state-ranking.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/9impute11_average_model_result/9impute11_average_model.para0
#
bash ${scriptDir}make_PRcurve_peak_state-ranking2.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/9impute11_average_model_result/9impute11_average_model.para0


# 8. 9impute11_average_motif_model
mkdir state_bed_9impute11_average_motif_model
bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh /data/zusers/fankaili/ideas/dhs_ctcf/9impute11_average_motif_model_result/ 9impute11_average_motif_model. /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_9impute11_average_motif_model/ 11
#
prefix="9impute11_average_motif_model"
bash ${scriptDir}make_CTCF_signal_based_curve.sh ${prefix}
#
bash ${scriptDir}make_PRcurve_peak_state-ranking.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/9impute11_average_motif_model_result/9impute11_average_motif_model.para0
#
bash ${scriptDir}make_PRcurve_peak_state-ranking2.sh ${prefix} /data/zusers/fankaili/ideas/dhs_ctcf/9impute11_average_motif_model_result/9impute11_average_motif_model.para0


# Rscript get_PRAU_barplot.R
