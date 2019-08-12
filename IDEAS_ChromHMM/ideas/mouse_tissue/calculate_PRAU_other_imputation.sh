#!/bin/bash

# -- Kaili
# This script is for calculating PRAU for other computation tools.
###
# for each 25 bp bin, whole bin inside a CTCF peak are considered as positive.
# bins are ranked by predicted CTCF signal.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/imputation_comparison/"

cd ${workDir}

# 0. pos & neg
intersectBed -a mm10_chr19_25bp.bed -b /data/zusers/fankaili/ideas/dhs_ctcf/peaks_validation_9To11/liver_14.5_ctcf_peak_sorted.bed -f 1 -u > mm10_chr19_25bp_liver_14.5_pos.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' mm10_chr19_25bp_liver_14.5_pos.bed mm10_chr19_25bp.bed > mm10_chr19_25bp_liver_14.5_neg.bed
#
intersectBed -a mm10_chr19_25bp.bed -b /data/zusers/fankaili/ideas/dhs_ctcf/peaks_validation_9To11/lung_14.5_ctcf_peak_sorted.bed -f 1 -u > mm10_chr19_25bp_lung_14.5_pos.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' mm10_chr19_25bp_lung_14.5_pos.bed mm10_chr19_25bp.bed > mm10_chr19_25bp_lung_14.5_neg.bed
#
sample="forebrain_0"
intersectBed -a mm10_chr19_25bp.bed -b /data/zusers/fankaili/ideas/dhs_ctcf/peaks_validation_9To11/${sample}_ctcf_peak_sorted.bed -f 1 -u > mm10_chr19_25bp_${sample}_pos.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' mm10_chr19_25bp_${sample}_pos.bed mm10_chr19_25bp.bed > mm10_chr19_25bp_${sample}_neg.bed
#
sample="midbrain_0"
intersectBed -a mm10_chr19_25bp.bed -b /data/zusers/fankaili/ideas/dhs_ctcf/peaks_validation_9To11/${sample}_ctcf_peak_sorted.bed -f 1 -u > mm10_chr19_25bp_${sample}_pos.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' mm10_chr19_25bp_${sample}_pos.bed mm10_chr19_25bp.bed > mm10_chr19_25bp_${sample}_neg.bed
#
sample="hindbrain_0"
intersectBed -a mm10_chr19_25bp.bed -b /data/zusers/fankaili/ideas/dhs_ctcf/peaks_validation_9To11/${sample}_ctcf_peak_sorted.bed -f 1 -u > mm10_chr19_25bp_${sample}_pos.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' mm10_chr19_25bp_${sample}_pos.bed mm10_chr19_25bp.bed > mm10_chr19_25bp_${sample}_neg.bed


# 1. ChromImpute
## get bed signal
awk '{FS=OFS="\t"}{if(NR==FNR){a[NR-2]=$1}else{n=$2/25;print $0,a[n]}}' ./chromimpute/output/chr19_impute_liver_14.5_CTCF.wig mm10_chr19_25bp.bed > ./chromimpute/output/chromimpute_predict_liver_14.5.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[NR-2]=$1}else{n=$2/25;print $0,a[n]}}' ./chromimpute/output/chr19_impute_lung_14.5_CTCF.wig mm10_chr19_25bp.bed > ./chromimpute/output/chromimpute_predict_lung_14.5.bed
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_liver_14.5_pos.bed ./chromimpute/output/chromimpute_predict_liver_14.5.bed > ./chromimpute/output/chromimpute_predict_liver_14.5_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_liver_14.5_neg.bed ./chromimpute/output/chromimpute_predict_liver_14.5.bed > ./chromimpute/output/chromimpute_predict_liver_14.5_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}chromimpute/output/ chromimpute_predict_liver_14.5_pos.txt chromimpute_predict_liver_14.5_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e liver_14.5"\t"$tmp >> ./chromimpute/chromimpute_PRAU.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_lung_14.5_pos.bed ./chromimpute/output/chromimpute_predict_lung_14.5.bed > ./chromimpute/output/chromimpute_predict_lung_14.5_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_lung_14.5_neg.bed ./chromimpute/output/chromimpute_predict_lung_14.5.bed > ./chromimpute/output/chromimpute_predict_lung_14.5_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}chromimpute/output/ chromimpute_predict_lung_14.5_pos.txt chromimpute_predict_lung_14.5_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e lung_14.5"\t"$tmp >> ./chromimpute/chromimpute_PRAU.txt


# 2. Avocado
paste mm10_chr19_25bp.bed ./avocado/avocado_predict_liver_14.5_CTCF.txt > ./avocado/avocado_predict_liver_14.5_CTCF.bed
paste mm10_chr19_25bp.bed ./avocado/avocado_predict_lung_14.5_CTCF.txt > ./avocado/avocado_predict_lung_14.5_CTCF.bed
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_liver_14.5_pos.bed ./avocado/avocado_predict_liver_14.5_CTCF.bed > ./avocado/avocado_predict_liver_14.5_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_liver_14.5_neg.bed ./avocado/avocado_predict_liver_14.5_CTCF.bed > ./avocado/avocado_predict_liver_14.5_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}avocado/ avocado_predict_liver_14.5_pos.txt avocado_predict_liver_14.5_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e liver_14.5"\t"$tmp >> ./avocado/avocado_PRAU.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_lung_14.5_pos.bed ./avocado/avocado_predict_lung_14.5_CTCF.bed > ./avocado/avocado_predict_lung_14.5_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_lung_14.5_neg.bed ./avocado/avocado_predict_lung_14.5_CTCF.bed > ./avocado/avocado_predict_lung_14.5_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}avocado/ avocado_predict_lung_14.5_pos.txt avocado_predict_lung_14.5_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e lung_14.5"\t"$tmp >> ./avocado/avocado_PRAU.txt

# 3. ChromImpute - with CTCFmotif & CTCFaverage
## get bed signal
awk '{FS=OFS="\t"}{if(NR==FNR){a[NR-2]=$1}else{n=$2/25;print $0,a[n]}}' ./chromimpute/output2/chr19_impute_liver_14.5_CTCF.wig mm10_chr19_25bp.bed > ./chromimpute/output2/chromimpute_predict_liver_14.5.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[NR-2]=$1}else{n=$2/25;print $0,a[n]}}' ./chromimpute/output2/chr19_impute_lung_14.5_CTCF.wig mm10_chr19_25bp.bed > ./chromimpute/output2/chromimpute_predict_lung_14.5.bed
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_liver_14.5_pos.bed ./chromimpute/output2/chromimpute_predict_liver_14.5.bed > ./chromimpute/output2/chromimpute_predict_liver_14.5_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_liver_14.5_neg.bed ./chromimpute/output2/chromimpute_predict_liver_14.5.bed > ./chromimpute/output2/chromimpute_predict_liver_14.5_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}chromimpute/output2/ chromimpute_predict_liver_14.5_pos.txt chromimpute_predict_liver_14.5_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e liver_14.5"\t"$tmp >> ./chromimpute/chromimpute_PRAU.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_lung_14.5_pos.bed ./chromimpute/output2/chromimpute_predict_lung_14.5.bed > ./chromimpute/output2/chromimpute_predict_lung_14.5_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_lung_14.5_neg.bed ./chromimpute/output2/chromimpute_predict_lung_14.5.bed > ./chromimpute/output2/chromimpute_predict_lung_14.5_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}chromimpute/output2/ chromimpute_predict_lung_14.5_pos.txt chromimpute_predict_lung_14.5_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e lung_14.5"\t"$tmp >> ./chromimpute/chromimpute_PRAU.txt

# 4. chromImpute 8-11
outputDir="./chromimpute/output3/"
#
sample="forebrain_0"
awk '{FS=OFS="\t"}{if(NR==FNR){a[NR-2]=$1}else{n=$2/25;print $0,a[n]}}' ${outputDir}chr19_impute_${sample}_CTCF.wig mm10_chr19_25bp.bed > ${outputDir}chromimpute_predict_${sample}.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_pos.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_neg.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_pos.txt chromimpute_predict_${sample}_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_PRAU.txt
#
sample="midbrain_0"
awk '{FS=OFS="\t"}{if(NR==FNR){a[NR-2]=$1}else{n=$2/25;print $0,a[n]}}' ${outputDir}chr19_impute_${sample}_CTCF.wig mm10_chr19_25bp.bed > ${outputDir}chromimpute_predict_${sample}.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_pos.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_neg.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_pos.txt chromimpute_predict_${sample}_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_PRAU.txt
#
sample="hindbrain_0"
awk '{FS=OFS="\t"}{if(NR==FNR){a[NR-2]=$1}else{n=$2/25;print $0,a[n]}}' ${outputDir}chr19_impute_${sample}_CTCF.wig mm10_chr19_25bp.bed > ${outputDir}chromimpute_predict_${sample}.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_pos.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_neg.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_pos.txt chromimpute_predict_${sample}_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_PRAU.txt


# 5. chromImpute 8-11, with CTCFmotif & CTCFaverage
outputDir="./chromimpute/output4/"
#
sample="forebrain_0"
awk '{FS=OFS="\t"}{if(NR==FNR){a[NR-2]=$1}else{n=$2/25;print $0,a[n]}}' ${outputDir}chr19_impute_${sample}_CTCF.wig mm10_chr19_25bp.bed > ${outputDir}chromimpute_predict_${sample}.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_pos.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_neg.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_pos.txt chromimpute_predict_${sample}_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_PRAU.txt
#
sample="midbrain_0"
awk '{FS=OFS="\t"}{if(NR==FNR){a[NR-2]=$1}else{n=$2/25;print $0,a[n]}}' ${outputDir}chr19_impute_${sample}_CTCF.wig mm10_chr19_25bp.bed > ${outputDir}chromimpute_predict_${sample}.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_pos.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_neg.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_pos.txt chromimpute_predict_${sample}_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_PRAU.txt
#
sample="hindbrain_0"
awk '{FS=OFS="\t"}{if(NR==FNR){a[NR-2]=$1}else{n=$2/25;print $0,a[n]}}' ${outputDir}chr19_impute_${sample}_CTCF.wig mm10_chr19_25bp.bed > ${outputDir}chromimpute_predict_${sample}.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_pos.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_pos.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5}}}' mm10_chr19_25bp_${sample}_neg.bed ${outputDir}chromimpute_predict_${sample}.bed > ${outputDir}chromimpute_predict_${sample}_neg.txt
prau=`Rscript ${scriptDir}calculate_PRAU.R ${workDir}${outputDir} chromimpute_predict_${sample}_pos.txt chromimpute_predict_${sample}_neg.txt`
tmp=`awk '{print $2}' <<< $prau`
echo -e ${sample}"\t"$tmp >> ./chromimpute/chromimpute_PRAU.txt
