#!/bin/bash

# -- Kaili
# This script is for running new IDEAS adding rep2.
# 1. run DHS-bins without CTCF using rep1&rep2
# 2. imputation with rep2
# 3. validation with rep2
# 4. why so different between rep1 and rep1-impute-rep2?

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"

# 1. run DHS-bins without CTCF using rep1&rep2
cd /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/
## get .input file
cp DHS_v3_100-400bp.input DHS_reps.input
sed -i 's/rep1/rep2/g' DHS_reps.input
cat DHS_v3_100-400bp.input >> DHS_reps.input
## get .sh file
cp DHS_v3_100-400bp.sh DHS_reps.sh
vim DHS_reps.sh
## get .parafile
cp DHS_v3_100-400bp.parafile DHS_reps.parafile
vim DHS_reps.parafile
##
nohup bash DHS_reps.sh > ./nohup/nohup.DHS_reps.out 2>&1&

## comparison
bash ${scriptDir}gene_regression_rep1_rep1-2.sh

# 2. imputation with rep2
bash ${scriptDir}test_imputation_by_rep2.sh

# 3. validation with rep2
bash ${scriptDir}validate_imputation_by_rep2.sh

# 4. why so different between rep1 and rep1-impute-rep2?
## 1) correlation between rep1 and rep2 signals
cd /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/
#
if [ -f rep1_rep2_pearson.txt ]; then rm rep1_rep2_pearson.txt; fi
#
while read line
do
    for mark in ATAC DNAme H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 CTCF
    do
        r=$(Rscript ${scriptDir}calculate_correlation.R /data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/${line}_${mark}_normal.txt /data/zusers/fankaili/ideas/signal/rep2_signal_normal_bins/${line}_${mark}_normal.txt "pearson" | awk '{print $2}')
        echo -e $line"\t"${mark}"\t"$r >> rep1_rep2_pearson.txt
    done
done < CTCF_sample_list.txt
# Rscript make_pearson_barplot.R
#
if [ -f rep1_rep2_pearson_dhs.txt ]; then rm rep1_rep2_pearson_dhs.txt; fi
#
while read line
do
    for mark in ATAC DNAme H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 CTCF
    do
        r=$(Rscript ${scriptDir}calculate_correlation.R /data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/${line}_${mark}_dhs.txt /data/zusers/fankaili/ideas/signal/rep2_signal_dhs_bins/${line}_${mark}_dhs.txt "pearson" | awk '{print $2}')
        echo -e $line"\t"${mark}"\t"$r >> rep1_rep2_pearson_dhs.txt
    done
done < CTCF_sample_list.txt
#
if [ -f rep1_rep2_spearman.txt ]; then rm rep1_rep2_spearman.txt; fi
#
while read line
do
    for mark in ATAC DNAme H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 CTCF
    do
        r=$(Rscript ${scriptDir}calculate_correlation.R /data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/${line}_${mark}_normal.txt /data/zusers/fankaili/ideas/signal/rep2_signal_normal_bins/${line}_${mark}_normal.txt "spearman" | awk '{print $2}')
        echo -e $line"\t"${mark}"\t"$r >> rep1_rep2_spearman.txt
    done
done < CTCF_sample_list.txt

## 2) scatter plot: randomly pick 10,000 bins, signal between rep1 and rep2
cd /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/
#
mkdir random_10k_signal
# pick 10,000 regions
shuf -i 1-13237130 -n 10000 | awk '{print "R"$1}' > random_10k_bins.txt
# get signal
for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC CTCF
do
    echo ${mark}
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$5}}}' random_10k_bins.txt /data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/kidney_0_${mark}_normal.tab > tmp_rep1.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$5}}}' random_10k_bins.txt /data/zusers/fankaili/ideas/signal/rep2_signal_normal_bins/kidney_0_${mark}_normal.tab > tmp_rep2.txt
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,a[$1]}}' tmp_rep1.txt tmp_rep2.txt > ./random_10k_signal/kidney_0_${mark}_random_10k.txt
done
#
mark="DNAme"
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$6}}}' random_10k_bins.txt /data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/kidney_0_${mark}_normal.tab > tmp_rep1.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$6}}}' random_10k_bins.txt /data/zusers/fankaili/ideas/signal/rep2_signal_normal_bins/kidney_0_${mark}_normal.tab > tmp_rep2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,a[$1]}}' tmp_rep1.txt tmp_rep2.txt > ./random_10k_signal/kidney_0_${mark}_random_10k.txt
# make scatter
# Rscrip make_scatter_10k_random_reps.R
## 3) correlation between rep1 and rep2 in 2k windows
cd /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/
### get 2kb bins
cut -f 1-4 mm10_nearest_OCR-center_matched_normal.bed | sort -k1,1 -k2,2n > tmp.bed
bedtools merge -i tmp.bed | awk '{FS=OFS="\t"}{print $1,$2,$2}' > tmp_merge.bed
bedtools flank -i tmp_merge.bed -g /home/fankaili/genome/mm10.chrom.sizes.clean -b 1000 | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{print $0,"R"NR}' > random_2k_region.bed
## calculate signal
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo ${sample}","${assay}
    #
    bigWigAverageOverBed /data/projects/encode/data/${expID}/${fileID}.bigWig /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/random_2k_region.bed /data/zusers/fankaili/ideas/tmp_random_2k_signal/${sample}_${assay}_rep1_2k.tab
    awk '{print $5}' /data/zusers/fankaili/ideas/tmp_random_2k_signal/${sample}_${assay}_rep1_2k.tab > /data/zusers/fankaili/ideas/tmp_random_2k_signal/${sample}_${assay}_rep1_2k.txt
done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_8HM_filelist.txt
###
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo ${sample}","${assay}
    #
    bigWigAverageOverBed /data/projects/encode/data/${expID}/${fileID}.bigWig /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/random_2k_region.bed /data/zusers/fankaili/ideas/tmp_random_2k_signal/${sample}_${assay}_rep2_2k.tab
    awk '{print $5}' /data/zusers/fankaili/ideas/tmp_random_2k_signal/${sample}_${assay}_rep2_2k.tab > /data/zusers/fankaili/ideas/tmp_random_2k_signal/${sample}_${assay}_rep2_2k.txt
done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep2_8HM_filelist.txt
## calculate correlation
if [ -f random_10k_rep1_rep2_pearson.txt ]; then rm random_10k_rep1_rep2_pearson.txt; fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    #
    r=$(Rscript ${scriptDir}calculate_correlation.R /data/zusers/fankaili/ideas/tmp_random_2k_signal/${sample}_${assay}_rep1_2k.txt /data/zusers/fankaili/ideas/tmp_random_2k_signal/${sample}_${assay}_rep2_2k.txt "pearson" | awk '{print $2}')
    echo -e $sample"\t"${assay}"\t"$r >> random_10k_rep1_rep2_pearson.txt
done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_8HM_filelist.txt
# Rscript make_pearson_barplot.R

# 5. imputation with give partial CTCF data
bash ${scriptDir}test_imputation_by_rep2_givenCTCF.sh

# 6. imputation with give partial CTCF data
bash ${scriptDir}test_imputation_by_rep2_partial_model.sh
