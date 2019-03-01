#!/bin/bash

# -- Kaili
# This script is for making dhs-bins aggregation plot using rep1 data.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/"

mkdir ${workDir}
cd ${workDir}

# 0. preparation
## 1) merge bed file for normal bins and OCR bins
### normal bins
awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2]}' \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_up.bed > tmp_OCR-center_matched_normal_random_50windows.bed
awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2]}' \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_bins.bed >> tmp_OCR-center_matched_normal_random_50windows.bed
awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2]}' \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_down.bed >> tmp_OCR-center_matched_normal_random_50windows.bed
sort -k5,5 -k6,6n tmp_OCR-center_matched_normal_random_50windows.bed | cut -f 1-4 > \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_merged_sorted.bed
rm tmp_OCR-center_matched_normal_random_50windows.bed
### OCR bins
awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2]}' \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_OCR_random_50windows_up.bed > tmp_OCR-center_bins_v3_OCR_random_50windows.bed
awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2]}' \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_OCR_random_50windows_bins.bed >> tmp_OCR-center_bins_v3_OCR_random_50windows.bed
awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2]}' \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_OCR_random_50windows_down.bed >> tmp_OCR-center_bins_v3_OCR_random_50windows.bed
sort -k5,5 -k6,6n tmp_OCR-center_bins_v3_OCR_random_50windows.bed | cut -f 1-4 > \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_OCR_random_50windows_merged_sorted.bed
rm tmp_OCR-center_bins_v3_OCR_random_50windows.bed
## 2) get DNAme bigWig
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo ${sample}
    #
    gzip -dc /data/projects/encode/data/${expID}/${fileID}.bed.gz > tmp_${fileID}.bed
    bash /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/bedMethyl2bigWig.sh \
    tmp_${fileID}.bed /home/fankaili/genome/mm10.chrom.sizes /data/zusers/fankaili/ideas/signal/mouse_66samples_DNAme_bigWig/${fileID}.bigWig
    rm tmp_${fileID}.bed
done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_DNAme_filelist.txt

# while read line
# do
#     sample=`awk '{print $1}' <<< ${line}`
#     assay=`awk '{print $2}' <<< ${line}`
#     expID=`awk '{print $4}' <<< ${line}`
#     fileID=`awk '{print $5}' <<< ${line}`
#     #
#     echo ${sample}
#     #
#     echo "gzip -dc /data/projects/encode/data/"${expID}"/"${fileID}".bed.gz > tmp_"${fileID}".bed" >> get_DNAme_bw_code.sh
#     echo "bash /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/bedMethyl2bigWig.sh tmp_"${fileID}".bed /home/fankaili/genome/mm10.chrom.sizes /data/zusers/fankaili/ideas/signal/mouse_66samples_DNAme_bigWig/"${fileID}".bigWig"  >> get_DNAme_bw_code.sh
#     echo "rm tmp_"${fileID}".bed"  >> get_DNAme_bw_code.sh
# done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_DNAme_filelist.txt
#
# for i in {1..11}
# do
#     awk -v i="$i" '{if(NR>((i-1)*18) && NR<=(18*i)){print $0}}' get_DNAme_bw_code.sh > \
#     get_DNAme_bw_code_${i}.sh
#     nohup bash get_DNAme_bw_code_${i}.sh > nohup.get_DNAme_bw_code_${i}.out 2>&1&
# done

# 1. get signal
cat /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_8HM_filelist.txt \
/data/zusers/fankaili/ideas/ENCODE_mouse_rep1_ATAC_filelist.txt \
> /data/zusers/fankaili/ideas/tmp_ENCODE_mouse_rep1_8HM_ATAC_filelist.txt
#
get_signal_code(){
    bedFile=$1
    prefix=$2
    codefile=$3
    folder=$4
    encode_data_path="/data/projects/encode/data/"
    #
    if [ -f ${codefile} ]
    then
        rm ${codefile}
    fi
    #
    while read line
    do
        sample=`awk '{print $1}' <<< ${line}`
        assay=`awk '{print $2}' <<< ${line}`
        expID=`awk '{print $4}' <<< ${line}`
        fileID=`awk '{print $5}' <<< ${line}`
        #
        echo "bigWigAverageOverBed "${encode_data_path}${expID}"/"${fileID}".bigWig" ${bedFile} ${folder}$"/"${sample}"_"${assay}"_"${prefix}".tab" >> ${codefile}
        echo "awk '{FS=OFS=\"\\t\"}{if(NR==FNR){a[\$1]=\$5}else{print \$4,a[\$4]}}' "${folder}$"/"${sample}"_"${assay}"_"${prefix}".tab "${bedFile}" > "${folder}$"/"${sample}"_"${assay}"_"${prefix}".txt" >> ${codefile}
    done < /data/zusers/fankaili/ideas/tmp_ENCODE_mouse_rep1_8HM_ATAC_filelist.txt
    #
    while read line
    do
        sample=`awk '{print $1}' <<< ${line}`
        assay=`awk '{print $2}' <<< ${line}`
        expID=`awk '{print $4}' <<< ${line}`
        fileID=`awk '{print $5}' <<< ${line}`
        #
        echo "bigWigAverageOverBed /data/zusers/fankaili/ideas/signal/mouse_66samples_DNAme_bigWig/"${fileID}".bigWig" ${bedFile} ${folder}$"/"${sample}"_"${assay}"_"${prefix}".tab" >> ${codefile}
        echo "awk '{FS=OFS=\"\\t\"}{if(NR==FNR){if(\$3==0){a[\$1]=-10}else{a[\$1]=\$6}}else{print \$4,a[\$4]}}' "${folder}$"/"${sample}"_"${assay}"_"${prefix}".tab "${bedFile}" > "${folder}$"/"${sample}"_"${assay}"_"${prefix}".txt" >> ${codefile}
    done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_DNAme_filelist.txt
}

run_code(){
    codefile=$1
    codefile_name=$2
    #
    for i in {1..10}
    do
        awk -v i="$i" '{FS=OFS="\t"}{if(NR>((i-1)*132) && NR<=(132*i)){print $0}}' ${codefile} \
        > ./code/${codefile_name}_${i}.sh
        nohup bash ./code/${codefile_name}_${i}.sh > ./code/${codefile_name}_${i}.out 2>&1&
    done
}

mkdir code

### OCR bins
bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_OCR_random_50windows_merged_sorted.bed"
prefix="OCR"
codefile="get_signal_code_OCR_bins.sh"
codefile_name="get_signal_code_OCR_bins"
folder="signal_OCR_each_bin"
get_signal_code ${bedfile} ${prefix} ${codefile} ${folder}
mkdir ${folder}
run_code ${codefile} ${codefile_name}

###
# OCR matched normal bins
bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_merged_sorted.bed"
prefix="normal"
codefile="get_signal_code_normal_bins.sh"
codefile_name="get_signal_code_normal_bins"
folder="signal_normal_each_bin"
get_signal_code ${bedfile} ${prefix} ${codefile} ${folder}
mkdir ${folder}
run_code ${codefile} ${codefile_name}


# 2. make signal matrix
cd /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/
# OCR bins
mkdir signal_OCR
for file in `ls /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR_each_bin/*.txt`
do
    a=${file##/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR_each_bin/}
    filename=${a%_OCR.txt}
    echo ${filename}
    #
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%150!=1){printf "\t"$2}else{printf "\n"$2}}END{printf "\n"}' \
    ./signal_OCR_each_bin/${filename}_OCR.txt > temp_OCR.txt
    awk 'BEGIN{FS=OFS="\t";for(i=1;i<=150;i++){b[i]=0}}{for(i=1;i<=150;i++){if($i>=0){a[i]+=$i;b[i]+=1}}}END{for(i=1;i<=150;i++){print a[i]/b[i]}}' \
    temp_OCR.txt > ./signal_OCR/${filename}_OCR.txt
done
# normal bins
mkdir signal_normal
for file in `ls /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal_each_bin/*.txt`
do
    a=${file##/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_normal_each_bin/}
    filename=${a%_normal.txt}
    echo ${filename}
    #
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%150!=1){printf "\t"$2}else{printf "\n"$2}}END{printf "\n"}' \
    ./signal_normal_each_bin/${filename}_normal.txt > temp_normal.txt
    awk 'BEGIN{FS=OFS="\t";for(i=1;i<=150;i++){b[i]=0}}{for(i=1;i<=150;i++){if($i>=0){a[i]+=$i;b[i]+=1}}}END{for(i=1;i<=150;i++){print a[i]/b[i]}}' \
    temp_normal.txt > ./signal_normal/${filename}_normal.txt
done

# nohup bash ss1.sh > ./nohup.ss1.out 2>&1&
# nohup bash ss2.sh > ./nohup.ss2.out 2>&1&

# 3. make figures
cut -f 1 /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_DNAme_filelist.txt | sort -u | sort > \
/data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt
while read line
do
    echo $line
    Rscript ${scriptDir}make_aggregation_each_sample.R ${line}
done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt



################
# Feb27
# average aggregation
mkdir average_signal
#
for mark in ATAC DNAme H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3
do
    while read line
    do
        if [ $line == "embryonic-facial-prominence_11.5" ];then
            cat ./signal_OCR/${line}"_"${mark}_OCR.txt > ./average_signal/${mark}_average_signal_dhs.txt
            #
            cat ./signal_normal/${line}"_"${mark}_normal.txt > ./average_signal/${mark}_average_signal_normal.txt
        else
            paste ./signal_OCR/${line}"_"${mark}_OCR.txt ./average_signal/${mark}_average_signal_dhs.txt | awk '{print $1+$2}' > tmp_dhs.txt
            mv tmp_dhs.txt ./average_signal/${mark}_average_signal_dhs.txt
            #
            paste ./signal_normal/${line}"_"${mark}_normal.txt ./average_signal/${mark}_average_signal_normal.txt | awk '{print $1+$2}' > tmp_normal.txt
            mv tmp_normal.txt ./average_signal/${mark}_average_signal_normal.txt
        fi
    done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt
    awk '{print $1/66}' ./average_signal/${mark}_average_signal_dhs.txt > tmp_dhs.txt
    mv tmp_dhs.txt ./average_signal/${mark}_average_signal_dhs.txt
    #
    awk '{print $1/66}' ./average_signal/${mark}_average_signal_normal.txt > tmp_normal.txt
    mv tmp_normal.txt ./average_signal/${mark}_average_signal_normal.txt
done

# Rscript ${scriptDir}make_aggregation_average_all_sample.R
