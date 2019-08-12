#!/bin/bash

# -- Kaili
# This script is for calculating AUPR using tissue-specific CTCF peak as postive.
# DHSbins resolution

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/imputation_comparison/"

cd ${workDir}


# 0. get positive bins in all regions
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' /data/zusers/fankaili/ideas/dhs_ctcf/tissue_specific_peak/${sample}_specific_11CTCFpeaks.txt /data/zusers/fankaili/ideas/dhs_ctcf/peaks_validation_9To11/${sample}_ctcf_peak_center_withSignal.bed > tmp.bed
    intersectBed -a /data/zusers/fankaili/ideas/dhs_ctcf/mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed -b tmp.bed -wa -u > state_ranked_${sample}_specific_CTCFpeak_pos.txt
done
rm tmp.bed

# 1. chromImpute , same resolution
outputDir="./chromimpute_dhsResolution/output/"
#
sample="liver_14.5"
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_specific_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_specific_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
# calculate AUPR
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_specific_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute_dhsResolution/chromimpute_DHSbins_specific_PRAU.txt
# make PR AUC
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_specific_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt ChromImpute_${sample}_specific_PRcurve_peak
#
sample="lung_14.5"
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_specific_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_specific_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
# calculate AUPR
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_specific_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute_dhsResolution/chromimpute_DHSbins_specific_PRAU.txt
# make PR AUC
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_specific_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt ChromImpute_${sample}_specific_PRcurve_peak

# 2. chromImpute
outputDir="./chromimpute/output/"
#
sample="liver_14.5"
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_specific_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_specific_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
# calculate AUPR
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_specific_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_specific_PRAU.txt
# make PR AUC
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_specific_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt ChromImpute_${sample}_specific_PRcurve_peak
#
sample="lung_14.5"
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_specific_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_specific_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
# calculate AUPR
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_specific_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_specific_PRAU.txt
# make PR AUC
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_specific_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt ChromImpute_${sample}_specific_PRcurve_peak

# 3. Avocado, same resolution
outputDir="./avocado_dhsResolution/"
#
sample="liver_14.5"
bigWigAverageOverBed ${outputDir}avocado_predict_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_specific_CTCFpeak_pos.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt
# calculate AUPR
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./avocado_dhsResolution/avocado_DHSbins_specific_PRAU.txt
# make PR AUC
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt Avocado_${sample}_specific_PRcurve_peak
#
sample="lung_14.5"
bigWigAverageOverBed ${outputDir}avocado_predict_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_specific_CTCFpeak_pos.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt
# calculate AUPR
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./avocado_dhsResolution/avocado_DHSbins_specific_PRAU.txt
# make PR AUC
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt Avocado_${sample}_specific_PRcurve_peak

# 3. Avocado
outputDir="./avocado/"
#
sample="liver_14.5"
bigWigAverageOverBed ${outputDir}avocado_predict_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_specific_CTCFpeak_pos.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt
# calculate AUPR
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ${outputDir}avocado_DHSbins_specific_PRAU.txt
# make PR AUC
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt Avocado_${sample}_specific_PRcurve_peak
#
sample="lung_14.5"
bigWigAverageOverBed ${outputDir}avocado_predict_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_specific_CTCFpeak_pos.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt
# calculate AUPR
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ${outputDir}avocado_DHSbins_specific_PRAU.txt
# make PR AUC
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_specific_pos.txt avocado_predict_${sample}_CTCF_DHSbins_specific_neg.txt Avocado_${sample}_specific_PRcurve_peak


# 5. IDEAS
## 1) specific pos
cp /data/zusers/fankaili/ideas/imputation_comparison/state_ranked_*_specific_CTCFpeak_pos.txt /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/
## 2) calculate AUPR and make PR AUC
prefix="9impute11_average_motif_model"
parafile="/data/zusers/fankaili/ideas/dhs_ctcf/9impute11_average_motif_model_result/9impute11_average_motif_model.para0"
bash ${scriptDir}calculate_AUPR_specific.sh ${prefix} ${parafile} /data/zusers/fankaili/ideas/dhs_ctcf/
## 3) step 3, get PR AUC
for file in `ls /data/zusers/fankaili/ideas/dhs_ctcf/state_bed_${prefix}/*_state_sorted.bed`
do
    file0=${file%_state_sorted.bed}
    sample=${file0#/data/zusers/fankaili/ideas/dhs_ctcf/state_bed_${prefix}/}
    echo $sample
    #
    Rscript ${scriptDir}make_PRcurve_peak.R /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/${prefix}/ ${sample}_${prefix}_peak_PR_specific_pos.txt ${sample}_${prefix}_peak_PR_neg.txt ${sample}_${prefix}_PRcurve_specific_peak
done
