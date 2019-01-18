#!/bin/bash

# -- Kaili
# This script is for making dhs-bins aggregation plot using rep1 data.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/"

mkdir ${workDir}
cd ${workDir}


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
    rm ${codefile}
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
        echo "bigWigAverageOverBed "${encode_data_path}${expID}"/"${fileID}".bigWig" ${bedFile} ${folder}$"/"${sample}"_"${assay}"_"${prefix}".tab" >> ${codefile}
        echo "awk '{FS=OFS=\"\\t\"}{if(NR==FNR){a[\$1]=\$6}else{print \$4,a[\$4]}}' "${folder}$"/"${sample}"_"${assay}"_"${prefix}".tab "${bedFile}" > "${folder}$"/"${sample}"_"${assay}"_"${prefix}".txt" >> ${codefile}
    done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_DNAme_filelist.txt
}

get_signal(){
    codefile=$1
    codefile_name=$2
    #
    awk '{FS=OFS="\t"}{if(NR<1*133){print $0}}' ${codefile} > ./code/${codefile_name}_1.sh
    awk '{FS=OFS="\t"}{if(NR>(1*132) && NR<(2*132+1)){print $0}}' ${codefile} > ./code/${codefile_name}_2.sh
    awk '{FS=OFS="\t"}{if(NR>(2*132) && NR<(3*132+1)){print $0}}' ${codefile} > ./code/${codefile_name}_3.sh
    awk '{FS=OFS="\t"}{if(NR>(3*132) && NR<(4*132+1)){print $0}}' ${codefile} > ./code/${codefile_name}_4.sh
    awk '{FS=OFS="\t"}{if(NR>(4*132) && NR<(5*132+1)){print $0}}' ${codefile} > ./code/${codefile_name}_5.sh
    awk '{FS=OFS="\t"}{if(NR>(5*132) && NR<(6*132+1)){print $0}}' ${codefile} > ./code/${codefile_name}_6.sh
    awk '{FS=OFS="\t"}{if(NR>(6*132) && NR<(7*132+1)){print $0}}' ${codefile} > ./code/${codefile_name}_7.sh
    awk '{FS=OFS="\t"}{if(NR>(7*132) && NR<(8*132+1)){print $0}}' ${codefile} > ./code/${codefile_name}_8.sh
    awk '{FS=OFS="\t"}{if(NR>(8*132) && NR<(9*132+1)){print $0}}' ${codefile} > ./code/${codefile_name}_9.sh
    awk '{FS=OFS="\t"}{if(NR>(9*132) && NR<(10*132+1)){print $0}}' ${codefile} > ./code/${codefile_name}_10.sh
    #
    for i in {1..10}
    do
        nohup bash ./code/${codefile_name}_${i}.sh > ./code/${codefile_name}_${i}.out 2>&1&
    done
}

mkdir code

### OCR bins
bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_OCR_random_50windows_bins.bed"
prefix="v3_OCR_bins"
codefile="get_signal_code_OCR_bins.sh"
codefile_name="get_signal_code_OCR_bins"
folder="signal_OCR_bins"
get_signal_code ${bedfile} ${prefix} ${codefile} ${folder}
mkdir ${folder}
get_signal ${codefile} ${codefile_name}

bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_OCR_random_50windows_up.bed"
prefix="v3_OCR_up"
codefile="get_signal_code_OCR_up.sh"
codefile_name="get_signal_code_OCR_up"
folder="signal_OCR_up"
get_signal_code ${bedfile} ${prefix} ${codefile} ${folder}
mkdir ${folder}
get_signal ${codefile} ${codefile_name}

bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_OCR_random_50windows_down.bed"
prefix="v3_OCR_down"
codefile="get_signal_code_OCR_down.sh"
codefile_name="get_signal_code_OCR_down"
folder="signal_OCR_down"
get_signal_code ${bedfile} ${prefix} ${codefile} ${folder}
mkdir ${folder}
get_signal ${codefile} ${codefile_name}


###
# OCR matched normal bins
bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_bins.bed"
prefix="matched_normal_bins"
codefile="get_signal_code_matched_normal_bins.sh"
codefile_name="get_signal_code_matched_normal_bins"
folder="signal_matched_normal_bins"
get_signal_code ${bedfile} ${prefix} ${codefile} ${folder}
mkdir ${folder}
get_signal ${codefile} ${codefile_name}

bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_up.bed"
prefix="matched_normal_up"
codefile="get_signal_code_matched_normal_up.sh"
codefile_name="get_signal_code_matched_normal_up"
folder="signal_matched_normal_up"
get_signal_code ${bedfile} ${prefix} ${codefile} ${folder}
mkdir ${folder}
get_signal ${codefile} ${codefile_name}

bedfile="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_matched_normal_random_50windows_down.bed"
prefix="matched_normal_down"
codefile="get_signal_code_matched_normal_down.sh"
codefile_name="get_signal_code_matched_normal_down"
folder="signal_matched_normal_down"
get_signal_code ${bedfile} ${prefix} ${codefile} ${folder}
mkdir ${folder}
get_signal ${codefile} ${codefile_name}


# 2. make signal matrix
cd /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/
# normal bins
mkdir signal_matched_normal
for file in `ls /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_matched_normal_up/*.txt`
do
    a=${file##/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_matched_normal_up/}
    filename=${a%_matched_normal_up.txt}
    echo ${filename}
    #
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    ./signal_matched_normal_up/${filename}_matched_normal_up.txt > temp1.txt
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    ./signal_matched_normal_bins/${filename}_matched_normal_bins.txt > temp2.txt
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    ./signal_matched_normal_down/${filename}_matched_normal_down.txt > temp3.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/100489}}' temp1.txt > \
    ./signal_matched_normal/${filename}_matched_normal.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/100489}}' temp2.txt >> \
    ./signal_matched_normal/${filename}_matched_normal.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/100489}}' temp3.txt >> \
    ./signal_matched_normal/${filename}_matched_normal.txt
done
# OCR bins
mkdir signal_OCR
for file in `ls /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR_up/*.txt`
do
    a=${file##/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/signal_OCR_up/}
    filename=${a%_v3_OCR_up.txt}
    echo ${filename}
    #
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    ./signal_OCR_up/${filename}_v3_OCR_up.txt > temp4.txt
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    ./signal_OCR_bins/${filename}_v3_OCR_bins.txt > temp5.txt
    awk '{FS=OFS="\t"}{if(NR==1){printf $2}else if(NR%50!=1){printf "\t"$2}else{printf "\n"$2}}' \
    ./signal_OCR_down/${filename}_v3_OCR_down.txt > temp6.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/100000}}' temp4.txt > \
    ./signal_OCR/${filename}_OCR.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/100000}}' temp5.txt >> \
    ./signal_OCR/${filename}_OCR.txt
    awk '{FS=OFS}{for(i=1;i<=50;i++){a[i]+=$i}}END{for(i=1;i<=50;i++){print a[i]/100000}}' temp6.txt >> \
    ./signal_OCR/${filename}_OCR.txt
done

# 3. make figures
cut -f 1 /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_DNAme_filelist.txt | sort -u | sort > \
/data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt
while read line
do
    Rscript ${scriptDir}make_aggregation_each_sample.R ${line}
done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt
