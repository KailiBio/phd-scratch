#!/bin/bash

# -- Kaili
# This script is for calculating AUPR and making PR UAC for ChromImpute and Avocado in the DHS-bins resolution.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/imputation_comparison/"

cd ${workDir}

# 1. ChromImpute
outputDir="./chromimpute_dhsResolution/output/"
#
sample="liver_14.5"
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute_dhsResolution/chromimpute_DHSbins_PRAU.txt
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt ChromImpute_${sample}_PRcurve_peak
#
sample="lung_14.5"
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute_dhsResolution/chromimpute_DHSbins_PRAU.txt
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt ChromImpute_${sample}_PRcurve_peak


# 2. avocado
outputDir="./avocado_dhsResolution/"
#
sample="liver_14.5"
bigWigAverageOverBed ${outputDir}avocado_predict_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_pos.txt avocado_predict_${sample}_CTCF_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ${outputDir}avocado_DHSbins_PRAU.txt
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_pos.txt avocado_predict_${sample}_CTCF_DHSbins_neg.txt Avocado_${sample}_CTCF_PRcurve_peak
#
sample="lung_14.5"
bigWigAverageOverBed ${outputDir}avocado_predict_${sample}_CTCF.bw mm10_OCR_chr19_bins.bed ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}avocado_predict_${sample}_CTCF_DHSbins.tab > ${outputDir}avocado_predict_${sample}_CTCF_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_pos.txt avocado_predict_${sample}_CTCF_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ${outputDir}avocado_DHSbins_PRAU.txt
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} avocado_predict_${sample}_CTCF_DHSbins_pos.txt avocado_predict_${sample}_CTCF_DHSbins_neg.txt Avocado_${sample}_CTCF_PRcurve_peak
