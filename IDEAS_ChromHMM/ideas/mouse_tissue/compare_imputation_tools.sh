#!/bin/bash

# -- Kaili
# This script is for comparing IDEAS imputation with ChromImpute and Avocado.
# 1. run chromImpute
# 2. run Avocado
# 3. run IDEAS
# 4. comparison
# 5. add CTCFmotif & CTCFaverage
# 6. run chromImpute 8-11
# 7. make track hub

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/imputation_comparison/"

cd ${workDir}

# 0. prepare data
## 1) get chr19 regions, cut into 25bp bins
awk '{if($1=="chr19"){print $0}}' /data/zusers/fankaili/ideas/dhs_ctcf/mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed | sort -k1,1 -k2,2n > mm10_OCR_chr19_bins.bed
bedtools merge -i mm10_OCR_chr19_bins.bed > mm10_chr19.bed
bedtools makewindows -b mm10_chr19.bed -w 25 | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{print $0,"bin_"NR}' > mm10_chr19_25bp.bed
# 2,328,192 bins
## 2) calculate signal for 25bp bins
encode_data_path="/data/projects/encode/data/"
signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_chr19_25bp/"
mkdir ${signal_path}
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ideas/dhs_ctcf/all_ctcf_sample.txt tmp_ENCODE_mouse_rep1_filelist.txt | sort -k1,1 -k2,2 > ENCODE_mouse_rep1_CTCFsample_filelist1.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ideas/dhs_ctcf/all_ctcf_sample.txt ENCODE_mouse_rep1_DNAme_filelist.txt | sort -k1,1 -k2,2 > ENCODE_mouse_rep1_CTCFsample_filelist2.txt
#
if [ -f ENCODE_rep1_signal_25bp-bins_code.sh ]; then rm ENCODE_rep1_signal_25bp-bins_code.sh; fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "bigWigAverageOverBed "${encode_data_path}${expID}"/"${fileID}".bigWig /data/zusers/fankaili/ideas/imputation_comparison/mm10_chr19_25bp.bed "${signal_path}${sample}"_"${assay}"_chr19_25bp.tab" >> ENCODE_rep1_signal_25bp-bins_code.sh
    echo "awk '{print ""\$""5}' "${signal_path}${sample}"_"${assay}"_chr19_25bp.tab > "${signal_path}${sample}"_"${assay}"_chr19_25bp.txt" >> ENCODE_rep1_signal_25bp-bins_code.sh
done < ENCODE_mouse_rep1_CTCFsample_filelist1.txt
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "bigWigAverageOverBed /data/zusers/fankaili/ideas/signal/mouse_66samples_DNAme_bigWig/"${fileID}".bigWig /data/zusers/fankaili/ideas/imputation_comparison/mm10_chr19_25bp.bed "${signal_path}${sample}"_"${assay}"_chr19_25bp.tab" >> ENCODE_rep1_signal_25bp-bins_code.sh
    echo "awk '{print ""\$""6}' "${signal_path}${sample}"_"${assay}"_chr19_25bp.tab > "${signal_path}${sample}"_"${assay}"_chr19_25bp.txt" >> ENCODE_rep1_signal_25bp-bins_code.sh
done < ENCODE_mouse_rep1_CTCFsample_filelist2.txt
###
for i in {1..11}
do
    awk -v i="$i" '{if(NR>((i-1)*22) && NR<=(22*i)){print $0}}' ENCODE_rep1_signal_25bp-bins_code.sh > \
    /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_25bp-bins_code_${i}.sh
    nohup bash /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_25bp-bins_code_${i}.sh > \
    /data/zusers/fankaili/ideas/code/get_rep1_signal/nohup.ENCODE_rep1_signal_25bp-bins_code_${i}.out 2>&1&
done



# 1. run chromImpute
## 1) convert singal into wig
### has to be chr19_***.wig.gz format =_=#
for file in `ls /data/zusers/fankaili/ideas/signal/rep1_signal_chr19_25bp/*.txt`
do
    name0=${file%_chr19_25bp.txt}
    name=${name0#/data/zusers/fankaili/ideas/signal/rep1_signal_chr19_25bp/}
    echo $name
    #
    paste mm10_chr19_25bp.bed ${file} > tmp.txt
    #
    echo "track type=wiggle_0 name="${name} > /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${name}_25bp.wig
    echo "fixedStep  chrom=chr19 start=0 step=25 span=25" >> /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${name}_25bp.wig
    awk '{m=$2/25;a[m]=1;b[m]=$5}END{for(i=1;i<=2457263;i++){if(a[i]){print b[i]}else{print 0}}}' tmp.txt >> /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${name}_25bp.wig
    #
    gzip -c /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${name}_25bp.wig > /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data/chr19_${name}_25bp.wig.gz
done
## 2) run ChromImpute
## get master table
if [ -f ./chromimpute/master_table.txt ];then rm ./chromimpute/master_table.txt; fi
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC CTCF DNAme
    do
        echo -e ${sample}"\t"${mark}"\t"${sample}"_"${mark}"_25bp" >> ./chromimpute/master_table.txt
    done
done
vim ./chromimpute/master_table.txt
#
grep "chr19" /home/fankaili/genome/mm10.chrom.sizes > ./chromimpute/mm10_chr19_chrom.size
#
nohup bash ${scriptDir}run_chromimpute.sh > ./nohup/nohup.run_chromimpute.out 2>&1&
# z001 42264

# 2. run Avocado
## convert signal file into compressed numpy
python ${scriptDir}convert_signal_compressed_numpy.py
## run avocado using default
nohup python ${scriptDir}run_avocado.py > ./nohup/nohup.run_avocado.out 2>&1&
# z001 47900

# 3. run IDEAS
## 1) get signal
awk '{if($1=="chr19" || $1=="chr18"){print $0}}' /data/zusers/fankaili/ideas/dhs_ctcf/mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed | sort -k1,1 -k2,2n > ./IDEAS/mm10_OCR_chr18-19_bins.bed
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    echo ${sample}
    for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC CTCF
    do
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $5}}}' ./IDEAS/mm10_OCR_chr18-19_bins.bed /data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/${sample}_${mark}_dhs.tab > ./IDEAS/signal/${sample}_${mark}_dhs_signal_chr18-19.txt
    done
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $6}}}' ./IDEAS/mm10_OCR_chr18-19_bins.bed /data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/${sample}_DNAme_dhs.tab > ./IDEAS/signal/${sample}_DNAme_dhs_signal_chr18-19.txt
done
### get space delimited file for IDEAS
awk '{FS="\t";OFS=" "}{if(NR==FNR){a[$4]=$1;b[$4]=$2;c[$4]=$3}else{if(a[$1]){print a[$1],b[$1],c[$1],$1}}}' ./IDEAS/mm10_OCR_chr18-19_bins.bed /data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/forebrain_0_ATAC_dhs.tab > ./IDEAS/mm10_OCR_chr18-19_bins.txt
awk '{print $4}' ./IDEAS/mm10_OCR_chr18-19_bins.txt > ./IDEAS/mm10_OCR_chr18-19_bins_order.txt
### get CTCF average signal
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5}else{print a[$1]}}' /data/zusers/fankaili/ideas/dhs_ctcf/9sample_averageCTCFsignal.tab ./IDEAS/mm10_OCR_chr18-19_bins_order.txt > ./IDEAS/signal/averageCTCF_signal_chr18-19.txt
### get CTCF motif signal
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$1]){print a[$1]}else{print 0}}}' /data/zusers/fankaili/ideas/dhs_ctcf/CTCFmotif_bin_withScore.txt ./IDEAS/mm10_OCR_chr18-19_bins_order.txt > ./IDEAS/signal/CTCFmotif_signal_chr18-19.txt
###
cp -r /data/zusers/fankaili/ideas/dhs_ctcf/data ./IDEAS/
cp -r /data/zusers/fankaili/ideas/dhs_ctcf/bin ./IDEAS/
## 2) make files
# cp /data/zusers/fankaili/ideas/dhs_ctcf/11sample_average_motif_model.sh ./IDEAS/test_IDEAS_chr19.sh
# vim ./IDEAS/test_IDEAS_chr19.sh
# cp /data/zusers/fankaili/ideas/dhs_ctcf/11sample_average_motif_model.parafile ./IDEAS/test_IDEAS_chr19.parafile
# vim ./IDEAS/test_IDEAS_chr19.parafile
#
cp ./IDEAS/test_IDEAS_chr19.sh ./IDEAS/test_IDEAS_chr18-19.sh
vim ./IDEAS/test_IDEAS_chr18-19.sh
cp ./IDEAS/test_IDEAS_chr19.parafile ./IDEAS/test_IDEAS_chr18-19.parafile
vim ./IDEAS/test_IDEAS_chr18-19.parafile
###
if [ -f ./IDEAS/test_IDEAS_chr18-19.input ];then rm ./IDEAS/test_IDEAS_chr18-19.input;fi
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    echo ${sample}
    for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC DNAme
    do
        echo ${sample}" "${mark}" /data/zusers/fankaili/ideas/imputation_comparison/IDEAS/signal/"${sample}_${mark}_dhs_signal_chr18-19.txt >> ./IDEAS/test_IDEAS_chr18-19.input
    done
done
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 lung_0 stomach_0
do
    echo ${sample}" CTCF /data/zusers/fankaili/ideas/imputation_comparison/IDEAS/signal/"${sample}_CTCF_dhs_signal_chr18-19.txt >> ./IDEAS/test_IDEAS_chr18-19.input
done
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    echo ${sample}" CTCFmotif /data/zusers/fankaili/ideas/imputation_comparison/IDEAS/signal/CTCFmotif_signal_chr18-19.txt" >> ./IDEAS/test_IDEAS_chr18-19.input
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/imputation_comparison/IDEAS/signal/averageCTCF_signal_chr18-19.txt" >> ./IDEAS/test_IDEAS_chr18-19.input
done
####
nohup bash ./IDEAS/test_IDEAS_chr18-19.sh > ./nohup/nohup.test_IDEAS_chr18-19.out 2>&1&
# z008 39567
# chr19 z008 11766


# 4. comparison
bash ${scriptDir}calculate_PRAU_other_imputation.sh
bash ${scriptDir}calculate_PRAU_other_imputation2.sh

# 5. add CTCFmotif & CTCFaverage
## 1) CTCF motif
intersectBed -a mm10_chr19_25bp.bed -b /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/mm10_CTCF_motif_region.bed -wa -wb | cut -f 4,9 | sort -k1,1 -k2,2nr | awk 'BEGIN{FS=OFS="\t";bin=""}{if(NR==1){bin=$1;print $0}else{if(bin!=$1){bin=$1;print $0}}}' > mm10_chr19_25bp_bin_withCTCFmotifScore.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$4]){print a[$4]}else{print 0}}}' mm10_chr19_25bp_bin_withCTCFmotifScore.txt mm10_chr19_25bp.bed > mm10_chr19_25bp_bin_CTCFmotifScore.txt
#
python ${scriptDir}convert_file_npz.py /data/zusers/fankaili/ideas/imputation_comparison/mm10_chr19_25bp_bin_CTCFmotifScore.txt /data/zusers/fankaili/ideas/imputation_comparison/mm10_chr19_25bp_bin_CTCFmotifScore.npz

## 2) CTCF average
awk '{FS=OFS="\t"}{print $4,0}' mm10_chr19_25bp.bed > chr19_25bp_9sample_sumCTCFsignal.txt
#
while read sample
do
    echo $sample
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $1,$2+a[$1]}}' /data/zusers/fankaili/ideas/signal/rep1_signal_chr19_25bp/${sample}_CTCF_chr19_25bp.tab chr19_25bp_9sample_sumCTCFsignal.txt > tmp.txt
    mv tmp.txt chr19_25bp_9sample_sumCTCFsignal.txt
done < /data/zusers/fankaili/ideas/dhs_ctcf/9_ctcf_sample.txt
#
awk '{FS=OFS="\t"}{print $2/9}' chr19_25bp_9sample_sumCTCFsignal.txt > chr19_25bp_9sample_averageCTCFsignal.txt
#
python ${scriptDir}convert_file_npz.py /data/zusers/fankaili/ideas/imputation_comparison/chr19_25bp_9sample_averageCTCFsignal.txt /data/zusers/fankaili/ideas/imputation_comparison/chr19_25bp_9sample_averageCTCFsignal.npz

## 3) run imputation
nohup python ${scriptDir}run_avocado2.py > ./nohup/nohup.run_avocado2.out 2>&1&
# z001 47168
####
paste mm10_chr19_25bp.bed chr19_25bp_9sample_averageCTCFsignal.txt > tmp.txt
paste mm10_chr19_25bp.bed mm10_chr19_25bp_bin_CTCFmotifScore.txt > tmp2.txt
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    echo ${sample}
    #
    echo "track type=wiggle_0 name="${sample}"_CTCFaverage" > /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${sample}"_CTCFaverage"_25bp.wig
    echo "fixedStep  chrom=chr19 start=0 step=25 span=25" >> /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${sample}"_CTCFaverage"_25bp.wig
    awk '{m=$2/25;a[m]=1;b[m]=$5}END{for(i=1;i<=2457263;i++){if(a[i]){print b[i]}else{print 0}}}' tmp.txt >> /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${sample}"_CTCFaverage"_25bp.wig
    #
    gzip -c /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${sample}"_CTCFaverage"_25bp.wig > /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data/chr19_${sample}"_CTCFaverage"_25bp.wig.gz
    ###
    echo "track type=wiggle_0 name="${sample}"_CTCFmotif" > /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${sample}"_CTCFmotif"_25bp.wig
    echo "fixedStep  chrom=chr19 start=0 step=25 span=25" >> /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${sample}"_CTCFmotif"_25bp.wig
    awk '{m=$2/25;a[m]=1;b[m]=$5}END{for(i=1;i<=2457263;i++){if(a[i]){print b[i]}else{print 0}}}' tmp2.txt >> /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${sample}"_CTCFmotif"_25bp.wig
    #
    gzip -c /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data2/${sample}"_CTCFmotif"_25bp.wig > /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data/chr19_${sample}"_CTCFmotif"_25bp.wig.gz
done
###############
cp ./chromimpute/master_table.txt ./chromimpute/master_table2.txt
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    for mark in CTCFmotif CTCFaverage
    do
        echo -e ${sample}"\t"${mark}"\t"${sample}"_"${mark}"_25bp" >> ./chromimpute/master_table2.txt
    done
done
#
nohup bash ${scriptDir}run_chromimpute2.sh > ./nohup/nohup.run_chromimpute2.out 2>&1&
# z001 59247


# 6. run chromImpute 8-11
## 1) mark-only
if [ -f ./chromimpute/master_table3.txt ]; then rm ./chromimpute/master_table3.txt; fi
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC CTCF DNAme
    do
        echo -e ${sample}"\t"${mark}"\t"${sample}"_"${mark}"_25bp" >> ./chromimpute/master_table3.txt
    done
done
#
vim ./chromimpute/master_table3.txt
#
nohup bash ${scriptDir}run_chromimpute3.sh > ./nohup/nohup.run_chromimpute3.out 2>&1&
# z001 3912
## 2) with CTCFmotif & CTCFaverage
if [ -f ./chromimpute/master_table4.txt ]; then rm ./chromimpute/master_table4.txt; fi
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    for mark in CTCFmotif CTCFaverage
    do
        echo -e ${sample}"\t"${mark}"\t"${sample}"_"${mark}"_25bp" >> ./chromimpute/master_table4.txt
    done
    #
    for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC CTCF DNAme
    do
        echo -e ${sample}"\t"${mark}"\t"${sample}"_"${mark}"_25bp" >> ./chromimpute/master_table4.txt
    done
done
#
vim ./chromimpute/master_table4.txt
#
nohup bash ${scriptDir}run_chromimpute4.sh > ./nohup/nohup.run_chromimpute4.out 2>&1&
# z001 4062



# 7. make track hub
## 1) get bigWig file
cp /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/output/chr19_impute_liver_14.5_CTCF.bigWig /data/public_html_users/fankaili/IDEAS/other_imputation1/chromimpute_liver_14.5.bigWig
cp /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/output/chr19_impute_lung_14.5_CTCF.bigWig /data/public_html_users/fankaili/IDEAS/other_imputation1/chromimpute_lung_14.5.bigWig
#
cp /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/output2/chr19_impute_liver_14.5_CTCF.bigWig /data/public_html_users/fankaili/IDEAS/other_imputation1/chromimpute_CTCFmotif-average_liver_14.5.bigWig
cp /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/output2/chr19_impute_lung_14.5_CTCF.bigWig /data/public_html_users/fankaili/IDEAS/other_imputation1/chromimpute_CTCFmotif-average_lung_14.5.bigWig
####
cp /data/zusers/fankaili/ideas/imputation_comparison/avocado/avocado_predict_liver_14.5_CTCF.bigWig /data/public_html_users/fankaili/IDEAS/other_imputation1/avocado_liver_14.5.bigWig
cp /data/zusers/fankaili/ideas/imputation_comparison/avocado/avocado_predict_lung_14.5_CTCF.bigWig /data/public_html_users/fankaili/IDEAS/other_imputation1/avocado_lung_14.5.bigWig

## 2) get files
cp /data/public_html_users/fankaili/IDEAS/CTCF_signal_hub/genomes_e14.5p0_CTCF_signal.txt /data/public_html_users/fankaili/IDEAS/other_imputation1/genomes.txt
vim /data/public_html_users/fankaili/IDEAS/other_imputation1/genomes.txt
cp /data/public_html_users/fankaili/IDEAS/CTCF_signal_hub/hub_e14.5p0_CTCF_signal.txt /data/public_html_users/fankaili/IDEAS/other_imputation1/hub.txt
vim /data/public_html_users/fankaili/IDEAS/other_imputation1/hub.txt
vim /data/public_html_users/fankaili/IDEAS/other_imputation1/trackDb.txt
