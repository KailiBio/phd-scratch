#!/bin/bash

# -- Kaili
# This script is for normalizing bins signal value, then rerun IDEAS.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/"

cd ${workDir}


#--------------------------------------------------------------------
## for v1_100-300bp_bins
cd ${workDir}v1_100_300bp/

# 1. do Normalization
if [ -f /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp_nor.input ]
then
    rm /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp_nor.input
else
    touch /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp_nor.input
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
        python ${scriptDir}zscore-normalization_norByLength.py ${tab_path} 6 > \
        /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_nor/${sample}"_"${mark}.tab
    else
        python ${scriptDir}zscore-normalization_norByLength.py ${tab_path} 5 > \
        /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_nor/${sample}"_"${mark}.tab
    fi
    sed -i 's/ \t/\t/g' /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_nor/${sample}"_"${mark}.tab
    #
    cut -f 2 /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_nor/${sample}"_"${mark}.tab > \
    /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_nor/${sample}"_"${mark}.txt
    echo ${sample}" "${mark}" /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/signal_nor/"${sample}"_"${mark}".txt" >> \
    /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp_nor.input
done < /data/zusers/fankaili/ideas/dhs_bins/v1_100_300bp/DHS_v1_100-300bp.input

# 2. get all files
cp DHS_v1_100-300bp.parafile DHS_v1_100-300bp_nor.parafile
vim DHS_v1_100-300bp_nor.parafile
cp DHS_v1_100-300bp.sh DHS_v1_100-300bp_nor.sh
vim DHS_v1_100-300bp_nor.sh

# 3. run IDEAS
nohup bash DHS_v1_100-300bp_nor.sh > ./nohup.DHS_v1_100-300bp_nor.out 2>&1&


#--------------------------------------------------------------------
## for v2_1-300bp_bins
cd ${workDir}v2_1_300bp/

# 1. do Normalization
if [ -f /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/DHS_v2_1-300bp_nor.input ]
then
    rm /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/DHS_v2_1-300bp_nor.input
else
    touch /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/DHS_v2_1-300bp_nor.input
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
        python ${scriptDir}zscore-normalization_norByLength.py ${tab_path} 6 > \
        /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_nor/${sample}"_"${mark}.tab
    else
        python ${scriptDir}zscore-normalization_norByLength.py ${tab_path} 5 > \
        /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_nor/${sample}"_"${mark}.tab
    fi
    sed -i 's/ \t/\t/g' /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_nor/${sample}"_"${mark}.tab
    #
    cut -f 2 /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_nor/${sample}"_"${mark}.tab > \
    /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_nor/${sample}"_"${mark}.txt
    echo ${sample}" "${mark}" /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/signal_nor/"${sample}"_"${mark}".txt" >> \
    /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/DHS_v2_1-300bp_nor.input
done < /data/zusers/fankaili/ideas/dhs_bins/v2_1_300bp/DHS_v2_1-300bp.input

# 2. get all files
cp DHS_v2_1-300bp.parafile DHS_v2_1-300bp_nor.parafile
vim DHS_v2_1-300bp_nor.parafile
cp DHS_v2_1-300bp.sh DHS_v2_1-300bp_nor.sh
vim DHS_v2_1-300bp_nor.sh

# 3. run IDEAS
nohup bash DHS_v2_1_300bp_nor.sh > ./nohup.DHS_v2_1_300bp_nor.out 2>&1&
