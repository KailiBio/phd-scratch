#!/bin/bash

# -- Kaili
# This script is for calculating z-score for given bin files.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/"

cd ${workDir}


#--------------------------------------------------------------------
## for v1_100-300bp_bins
cd ${workDir}v1_100_300bp/

# 1. do Normalization
if [ -f /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp_zscore.input ]
then
    rm /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp_zscore.input
else
    touch /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp_zscore.input
fi
#
while read line
do
    echo ${line}
    sample=`awk '{print $1}' <<< ${line}`
    mark=`awk '{print $2}' <<< ${line}`
    path=`awk '{print $3}' <<< ${line}`
    tab_path=${path//txt/tab}
    #
    if [ mark == "DNAme" ]
    then
        echo "python "${scriptDir}"zscore-normalization_norNotLog.py "${tab_path}" 6 > /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_zscore/"${sample}"_"${mark}".tab" >> v1_zscore_code.txt
    else
        echo "python "${scriptDir}"zscore-normalization_norNotLog.py "${tab_path}" 5 > /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_zscore/"${sample}"_"${mark}".tab" >> v1_zscore_code.txt
    fi
    echo "sed -i 's/ \t/\t/g' /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_zscore/"${sample}"_"${mark}".tab" >> v1_zscore_code.txt
    #
    echo "cut -f 2 /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_zscore/"${sample}"_"${mark}".tab > /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_zscore/"${sample}"_"${mark}".txt" >> v1_zscore_code.txt
    echo ${sample}" "${mark}" /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_zscore/"${sample}"_"${mark}".txt" >> \
    /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp_zscore.input
done < /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp.input

awk '{if(NR<=220){print $0}}' v1_zscore_code.txt > ./code/v1_zscore_code_1.txt
awk '{if(NR>220 && NR<=(220*2)){print $0}}' v1_zscore_code.txt > ./code/v1_zscore_code_2.txt
awk '{if(NR>220*2 && NR<=(220*3)){print $0}}' v1_zscore_code.txt > ./code/v1_zscore_code_3.txt
awk '{if(NR>220*3 && NR<=(220*4)){print $0}}' v1_zscore_code.txt > ./code/v1_zscore_code_4.txt
awk '{if(NR>220*4 && NR<=(220*5)){print $0}}' v1_zscore_code.txt > ./code/v1_zscore_code_5.txt
awk '{if(NR>220*5 && NR<=(220*6)){print $0}}' v1_zscore_code.txt > ./code/v1_zscore_code_6.txt
awk '{if(NR>220*6 && NR<=(220*7)){print $0}}' v1_zscore_code.txt > ./code/v1_zscore_code_7.txt
awk '{if(NR>220*7 && NR<=(220*8)){print $0}}' v1_zscore_code.txt > ./code/v1_zscore_code_8.txt
awk '{if(NR>220*8 && NR<=(220*9)){print $0}}' v1_zscore_code.txt > ./code/v1_zscore_code_9.txt

nohup bash ./code/v1_zscore_code_1.txt > ./code/nohup.v1_zscore_code_1.out 2>&1&
nohup bash ./code/v1_zscore_code_2.txt > ./code/nohup.v1_zscore_code_2.out 2>&1&
nohup bash ./code/v1_zscore_code_3.txt > ./code/nohup.v1_zscore_code_3.out 2>&1&
nohup bash ./code/v1_zscore_code_4.txt > ./code/nohup.v1_zscore_code_4.out 2>&1&
nohup bash ./code/v1_zscore_code_5.txt > ./code/nohup.v1_zscore_code_5.out 2>&1&
nohup bash ./code/v1_zscore_code_6.txt > ./code/nohup.v1_zscore_code_6.out 2>&1&
nohup bash ./code/v1_zscore_code_7.txt > ./code/nohup.v1_zscore_code_7.out 2>&1&
nohup bash ./code/v1_zscore_code_8.txt > ./code/nohup.v1_zscore_code_8.out 2>&1&
nohup bash ./code/v1_zscore_code_9.txt > ./code/nohup.v1_zscore_code_9.out 2>&1&

# 2. get all files
cp DHS_v1_100-300bp_nor.parafile DHS_v1_100-300bp_zscore.parafile
vim DHS_v1_100-300bp_zscore.parafile
cp DHS_v1_100-300bp_nor.sh DHS_v1_100-300bp_zscore.sh
vim DHS_v1_100-300bp_zscore.sh

# 3. run IDEAS
nohup bash DHS_v1_100-300bp_zscore.sh > ./nohup.DHS_v1_100-300bp_zscore.out 2>&1&



#--------------------------------------------------------------------
## for v2_1-300bp_bins
cd ${workDir}v2_1_300bp/

# 1. do Normalization
if [ -f /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/DHS_v2_1-300bp_zscore.input ]
then
    rm /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/DHS_v2_1-300bp_zscore.input
else
    touch /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/DHS_v2_1-300bp_zscore.input
fi
#
while read line
do
    echo ${line}
    sample=`awk '{print $1}' <<< ${line}`
    mark=`awk '{print $2}' <<< ${line}`
    path=`awk '{print $3}' <<< ${line}`
    tab_path=${path//txt/tab}
    #
    if [ mark == "DNAme" ]
    then
        echo "python "${scriptDir}"zscore-normalization_norNotLog.py "${tab_path}" 6 > /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_zscore/"${sample}"_"${mark}".tab" >> v2_zscore_code.txt
    else
        echo "python "${scriptDir}"zscore-normalization_norNotLog.py "${tab_path}" 5 > /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_zscore/"${sample}"_"${mark}".tab" >> v2_zscore_code.txt
    fi
    echo "sed -i 's/ \t/\t/g' /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_zscore/"${sample}"_"${mark}".tab" >> v2_zscore_code.txt
    #
    echo "cut -f 2 /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_zscore/"${sample}"_"${mark}".tab > /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_zscore/"${sample}"_"${mark}".txt" >> v2_zscore_code.txt
    echo ${sample}" "${mark}" /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_zscore/"${sample}"_"${mark}".txt" >> \
    /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/DHS_v2_1-300bp_zscore.input
done < /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/DHS_v2_1-300bp.input

awk '{if(NR<=220){print $0}}' v2_zscore_code.txt > ./code/v2_zscore_code_1.txt
awk '{if(NR>220 && NR<=(220*2)){print $0}}' v2_zscore_code.txt > ./code/v2_zscore_code_2.txt
awk '{if(NR>220*2 && NR<=(220*3)){print $0}}' v2_zscore_code.txt > ./code/v2_zscore_code_3.txt
awk '{if(NR>220*3 && NR<=(220*4)){print $0}}' v2_zscore_code.txt > ./code/v2_zscore_code_4.txt
awk '{if(NR>220*4 && NR<=(220*5)){print $0}}' v2_zscore_code.txt > ./code/v2_zscore_code_5.txt
awk '{if(NR>220*5 && NR<=(220*6)){print $0}}' v2_zscore_code.txt > ./code/v2_zscore_code_6.txt
awk '{if(NR>220*6 && NR<=(220*7)){print $0}}' v2_zscore_code.txt > ./code/v2_zscore_code_7.txt
awk '{if(NR>220*7 && NR<=(220*8)){print $0}}' v2_zscore_code.txt > ./code/v2_zscore_code_8.txt
awk '{if(NR>220*8 && NR<=(220*9)){print $0}}' v2_zscore_code.txt > ./code/v2_zscore_code_9.txt

nohup bash ./code/v2_zscore_code_1.txt > ./code/nohup.v2_zscore_code_1.out 2>&1&
nohup bash ./code/v2_zscore_code_2.txt > ./code/nohup.v2_zscore_code_2.out 2>&1&
nohup bash ./code/v2_zscore_code_3.txt > ./code/nohup.v2_zscore_code_3.out 2>&1&
nohup bash ./code/v2_zscore_code_4.txt > ./code/nohup.v2_zscore_code_4.out 2>&1&
nohup bash ./code/v2_zscore_code_5.txt > ./code/nohup.v2_zscore_code_5.out 2>&1&
nohup bash ./code/v2_zscore_code_6.txt > ./code/nohup.v2_zscore_code_6.out 2>&1&
nohup bash ./code/v2_zscore_code_7.txt > ./code/nohup.v2_zscore_code_7.out 2>&1&
nohup bash ./code/v2_zscore_code_8.txt > ./code/nohup.v2_zscore_code_8.out 2>&1&
nohup bash ./code/v2_zscore_code_9.txt > ./code/nohup.v2_zscore_code_9.out 2>&1&


# 2. get all files
cp DHS_v2_1-300bp_nor.parafile DHS_v2_1-300bp_zscore.parafile
vim DHS_v2_1-300bp_zscore.parafile
cp DHS_v2_1-300bp_nor.sh DHS_v2_1-300bp_zscore.sh
vim DHS_v2_1-300bp_zscore.sh

# 3. run IDEAS
nohup bash DHS_v2_1-300bp_zscore.sh > ./nohup.DHS_v2_1-300bp_zscore.out 2>&1&
