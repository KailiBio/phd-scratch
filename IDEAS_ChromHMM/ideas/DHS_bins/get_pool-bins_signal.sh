#!/bin/bash

# -- Kaili
# This script is for getting Histone marks signal for pool-bins.

state=$1

poolDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/pool-bins/"

########
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3
do
    echo ${mark}
    #
    mkdir ${poolDir}state${state}_${mark}
    cd ${poolDir}state${state}_${mark}
    #
    echo "ID" > state${state}_${mark}_signalMatrix.txt
    cut -f 4 ${poolDir}dhs_ctcf_9-66_state_${state}_pool-bins_withID.bed >> state${state}_${mark}_signalMatrix.txt
    #
    while read line
    do
        sample=`awk '{print $1}' <<< ${line}`
        hm=`awk '{print $2}' <<< ${line}`
        expID=`awk '{print $4}' <<< ${line}`
        fileID=`awk '{print $5}' <<< ${line}`
        #
        if [ ${hm} == ${mark} ]
        then
            echo ${sample}
            # get signal file
            bigWigAverageOverBed /data/projects/encode/data/${expID}/${fileID}.bigWig ${poolDir}dhs_ctcf_9-66_state_${state}_pool-bins_withID.bed state${state}_${mark}_${sample}.tab
            echo -e "ID\ta\tb\tc\t"${sample}"\td" > tmp_${state}_${mark}_${sample}.txt
            cat state${state}_${mark}_${sample}.tab >> tmp_${state}_${mark}_${sample}.txt
            # # merged
            awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $0,a[$1]}}' tmp_${state}_${mark}_${sample}.txt state${state}_${mark}_signalMatrix.txt > tmp_${state}_${mark}_${sample}_merged.txt
            mv tmp_${state}_${mark}_${sample}_merged.txt state${state}_${mark}_signalMatrix.txt
        fi
    done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_8HM_filelist.txt
    #
    rm tmp_${state}_${mark}_*.txt
done
