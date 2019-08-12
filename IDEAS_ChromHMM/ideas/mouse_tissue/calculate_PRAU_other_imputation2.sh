#!/bin/bash

# -- Kaili
# This script is for calculating PRAU for other computation tools.
###
# use the predicted signal in 25 bp resolution to calculate signal for each OCR-center bins.
# OCR-center bins with a CTCF peak summit are regarded as postive.
# bins are ranked by average predicted signal.

# 1. chromImpute1
# 2. chromImpute2
# 3. chromImpute3
#4. chromImpute4
#5. Avocado


scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/imputation_comparison/"

cd ${workDir}

# 0. pos & neg
awk '$1=="chr19"' /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/liver_14.5_CTCFpeak_pos.txt > state_ranked_liver_14.5_CTCFpeak_pos.txt
awk '$1=="chr19"' /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/liver_14.5_CTCFpeak_neg.txt > state_ranked_liver_14.5_CTCFpeak_neg.txt
#
awk '$1=="chr19"' /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/lung_14.5_CTCFpeak_pos.txt > state_ranked_lung_14.5_CTCFpeak_pos.txt
awk '$1=="chr19"' /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/lung_14.5_CTCFpeak_neg.txt > state_ranked_lung_14.5_CTCFpeak_neg.txt
#
awk '$1=="chr19"' /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/forebrain_0_CTCFpeak_pos.txt > state_ranked_forebrain_0_CTCFpeak_pos.txt
awk '$1=="chr19"' /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/forebrain_0_CTCFpeak_neg.txt > state_ranked_forebrain_0_CTCFpeak_neg.txt
#
awk '$1=="chr19"' /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/midbrain_0_CTCFpeak_pos.txt > state_ranked_midbrain_0_CTCFpeak_pos.txt
awk '$1=="chr19"' /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/midbrain_0_CTCFpeak_neg.txt > state_ranked_midbrain_0_CTCFpeak_neg.txt
#
awk '$1=="chr19"' /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/hindbrain_0_CTCFpeak_pos.txt > state_ranked_hindbrain_0_CTCFpeak_pos.txt
awk '$1=="chr19"' /data/zusers/fankaili/ideas/dhs_ctcf/state_ranked_PR/hindbrain_0_CTCFpeak_neg.txt > state_ranked_hindbrain_0_CTCFpeak_neg.txt


# 1. chromImpute1
outputDir="./chromimpute/output/"
#
sample="liver_14.5"
wigToBigWig ${outputDir}chr19_impute_liver_14.5_CTCF.wig ./chromimpute/mm10_chr19_chrom.size2 ${outputDir}chr19_impute_liver_14.5_CTCF.bigWig
bigWigAverageOverBed ${outputDir}chr19_impute_liver_14.5_CTCF.bigWig mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_liver_14.5_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_liver_14.5_CTCFpeak_pos.txt ${outputDir}chr19_impute_liver_14.5_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_liver_14.5_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_liver_14.5_CTCFpeak_neg.txt ${outputDir}chr19_impute_liver_14.5_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_liver_14.5_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_liver_14.5_DHSbins_pos.txt chromimpute_predict_liver_14.5_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e liver_14.5"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_PRAU.txt
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt ChromImpute_${sample}_PRcurve_peak
#
sample="lung_14.5"
wigToBigWig ${outputDir}chr19_impute_lung_14.5_CTCF.wig ./chromimpute/mm10_chr19_chrom.size2 ${outputDir}chr19_impute_lung_14.5_CTCF.bigWig
bigWigAverageOverBed ${outputDir}chr19_impute_lung_14.5_CTCF.bigWig mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_lung_14.5_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_lung_14.5_CTCFpeak_pos.txt ${outputDir}chr19_impute_lung_14.5_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_lung_14.5_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_lung_14.5_CTCFpeak_neg.txt ${outputDir}chr19_impute_lung_14.5_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_lung_14.5_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_lung_14.5_DHSbins_pos.txt chromimpute_predict_lung_14.5_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e lung_14.5"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_PRAU.txt
Rscript ${scriptDir}make_PRcurve_peak.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt ChromImpute_${sample}_PRcurve_peak

# 2. chromImpute2
outputDir="./chromimpute/output2/"
#
wigToBigWig ${outputDir}chr19_impute_liver_14.5_CTCF.wig ./chromimpute/mm10_chr19_chrom.size2 ${outputDir}chr19_impute_liver_14.5_CTCF.bigWig
bigWigAverageOverBed ${outputDir}chr19_impute_liver_14.5_CTCF.bigWig mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_liver_14.5_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_liver_14.5_CTCFpeak_pos.txt ${outputDir}chr19_impute_liver_14.5_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_liver_14.5_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_liver_14.5_CTCFpeak_neg.txt ${outputDir}chr19_impute_liver_14.5_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_liver_14.5_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_liver_14.5_DHSbins_pos.txt chromimpute_predict_liver_14.5_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e liver_14.5"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_PRAU.txt
#
wigToBigWig ${outputDir}chr19_impute_lung_14.5_CTCF.wig ./chromimpute/mm10_chr19_chrom.size2 ${outputDir}chr19_impute_lung_14.5_CTCF.bigWig
bigWigAverageOverBed ${outputDir}chr19_impute_lung_14.5_CTCF.bigWig mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_lung_14.5_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_lung_14.5_CTCFpeak_pos.txt ${outputDir}chr19_impute_lung_14.5_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_lung_14.5_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_lung_14.5_CTCFpeak_neg.txt ${outputDir}chr19_impute_lung_14.5_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_lung_14.5_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_lung_14.5_DHSbins_pos.txt chromimpute_predict_lung_14.5_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e lung_14.5"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_PRAU.txt


# 3. chromImpute3
outputDir="./chromimpute/output3/"
#
sample="forebrain_0"
wigToBigWig ${outputDir}chr19_impute_${sample}_CTCF.wig ./chromimpute/mm10_chr19_chrom.size2 ${outputDir}chr19_impute_${sample}_CTCF.bigWig
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bigWig mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_PRAU.txt
#
sample="midbrain_0"
wigToBigWig ${outputDir}chr19_impute_${sample}_CTCF.wig ./chromimpute/mm10_chr19_chrom.size2 ${outputDir}chr19_impute_${sample}_CTCF.bigWig
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bigWig mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_PRAU.txt
#
sample="hindbrain_0"
wigToBigWig ${outputDir}chr19_impute_${sample}_CTCF.wig ./chromimpute/mm10_chr19_chrom.size2 ${outputDir}chr19_impute_${sample}_CTCF.bigWig
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bigWig mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_PRAU.txt


#4. chromImpute4
outputDir="./chromimpute/output4/"
#
sample="forebrain_0"
wigToBigWig ${outputDir}chr19_impute_${sample}_CTCF.wig ./chromimpute/mm10_chr19_chrom.size2 ${outputDir}chr19_impute_${sample}_CTCF.bigWig
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bigWig mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_PRAU.txt
#
sample="midbrain_0"
wigToBigWig ${outputDir}chr19_impute_${sample}_CTCF.wig ./chromimpute/mm10_chr19_chrom.size2 ${outputDir}chr19_impute_${sample}_CTCF.bigWig
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bigWig mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_PRAU.txt
#
sample="hindbrain_0"
wigToBigWig ${outputDir}chr19_impute_${sample}_CTCF.wig ./chromimpute/mm10_chr19_chrom.size2 ${outputDir}chr19_impute_${sample}_CTCF.bigWig
bigWigAverageOverBed ${outputDir}chr19_impute_${sample}_CTCF.bigWig mm10_OCR_chr19_bins.bed ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ${outputDir}chr19_impute_${sample}_CTCF_DHSbins.tab > ${outputDir}chromimpute_predict_${sample}_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_DHSbins_pos.txt chromimpute_predict_${sample}_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_DHSbins_PRAU.txt


#5. Avocado
#
sample="liver_14.5"
awk '{FS=OFS="\t"}{print $1,$2,$3,$5}' ./avocado/avocado_predict_${sample}_CTCF.bed > ./avocado/avocado_predict_${sample}_CTCF.bedGraph
bedGraphToBigWig ./avocado/avocado_predict_${sample}_CTCF.bedGraph ./chromimpute/mm10_chr19_chrom.size2 ./avocado/avocado_predict_${sample}_CTCF.bigWig
bigWigAverageOverBed ./avocado/avocado_predict_${sample}_CTCF.bigWig mm10_OCR_chr19_bins.bed ./avocado/avocado_predict_${sample}_CTCF_DHSbins.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ./avocado/avocado_predict_${sample}_CTCF_DHSbins.tab > ./avocado/avocado_predict_${sample}_CTCF_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ./avocado/avocado_predict_${sample}_CTCF_DHSbins.tab > ./avocado/avocado_predict_${sample}_CTCF_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ./avocado/ avocado_predict_${sample}_CTCF_DHSbins_pos.txt avocado_predict_${sample}_CTCF_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./avocado/avocado_DHSbins_PRAU.txt
Rscript ${scriptDir}make_PRcurve_peak.R /data/zusers/fankaili/ideas/imputation_comparison/avocado/ avocado_predict_${sample}_CTCF_DHSbins_pos.txt avocado_predict_${sample}_CTCF_DHSbins_neg.txt Avocado_${sample}_CTCF_PRcurve_peak
#
sample="lung_14.5"
awk '{FS=OFS="\t"}{print $1,$2,$3,$5}' ./avocado/avocado_predict_${sample}_CTCF.bed > ./avocado/avocado_predict_${sample}_CTCF.bedGraph
bedGraphToBigWig ./avocado/avocado_predict_${sample}_CTCF.bedGraph ./chromimpute/mm10_chr19_chrom.size2 ./avocado/avocado_predict_${sample}_CTCF.bigWig
bigWigAverageOverBed ./avocado/avocado_predict_${sample}_CTCF.bigWig mm10_OCR_chr19_bins.bed ./avocado/avocado_predict_${sample}_CTCF_DHSbins.tab
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_pos.txt ./avocado/avocado_predict_${sample}_CTCF_DHSbins.tab > ./avocado/avocado_predict_${sample}_CTCF_DHSbins_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' state_ranked_${sample}_CTCFpeak_neg.txt ./avocado/avocado_predict_${sample}_CTCF_DHSbins.tab > ./avocado/avocado_predict_${sample}_CTCF_DHSbins_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ./avocado/ avocado_predict_${sample}_CTCF_DHSbins_pos.txt avocado_predict_${sample}_CTCF_DHSbins_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./avocado/avocado_DHSbins_PRAU.txt
Rscript ${scriptDir}make_PRcurve_peak.R /data/zusers/fankaili/ideas/imputation_comparison/avocado/ avocado_predict_${sample}_CTCF_DHSbins_pos.txt avocado_predict_${sample}_CTCF_DHSbins_neg.txt Avocado_${sample}_CTCF_PRcurve_peak
