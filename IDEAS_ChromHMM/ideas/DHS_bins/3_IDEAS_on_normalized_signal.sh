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




#--------------------------------------------------------------------
## for v3_150-350bp_bins
cd ${workDir}v3_150_350bp/

# 1. do Normalization
if [ -f /data/zusers/fankaili/ideas/dhs_bins/v3_150_350bp/DHS_v3_150-350bp_nor.input ]
then
    rm /data/zusers/fankaili/ideas/dhs_bins/v3_150_350bp/DHS_v3_150-350bp_nor.input
else
    touch /data/zusers/fankaili/ideas/dhs_bins/v3_150_350bp/DHS_v3_150-350bp_nor.input
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
        echo "python "${scriptDir}"zscore-normalization_norByLength.py "${tab_path}" 6 > /data/zusers/fankaili/ideas/dhs_bins/v3_150_350bp/signal_nor/"${sample}"_"${mark}".tab" >> v3_zscore_length_code.txt
    else
        echo "python "${scriptDir}"zscore-normalization_norByLength.py "${tab_path}" 5 > /data/zusers/fankaili/ideas/dhs_bins/v3_150_350bp/signal_nor/"${sample}"_"${mark}".tab" >> v3_zscore_length_code.txt
    fi
    echo "sed -i 's/ \t/\t/g' /data/zusers/fankaili/ideas/dhs_bins/v3_150_350bp/signal_nor/"${sample}"_"${mark}".tab" >> v3_zscore_length_code.txt
    #
    echo "cut -f 2 /data/zusers/fankaili/ideas/dhs_bins/v3_150_350bp/signal_nor/"${sample}"_"${mark}".tab > /data/zusers/fankaili/ideas/dhs_bins/v3_150_350bp/signal_nor/"${sample}"_"${mark}".txt" >> v3_zscore_length_code.txt
    echo ${sample}" "${mark}" /data/zusers/fankaili/ideas/dhs_bins/v3_150_350bp/signal_nor/"${sample}"_"${mark}".txt" >> \
    /data/zusers/fankaili/ideas/dhs_bins/v3_150_350bp/DHS_v3_150-350bp_nor.input
done < /data/zusers/fankaili/ideas/dhs_bins/v3_150_3500bp/DHS_v3_150-350bp.input


awk '{if(NR<=220){print $0}}' v3_zscore_length_code.txt > ./code/v3_zscore_length_code_1.txt
awk '{if(NR>220 && NR<=(220*2)){print $0}}' v3_zscore_length_code.txt > ./code/v3_zscore_length_code_2.txt
awk '{if(NR>220*2 && NR<=(220*3)){print $0}}' v3_zscore_length_code.txt > ./code/v3_zscore_length_code_3.txt
awk '{if(NR>220*3 && NR<=(220*4)){print $0}}' v3_zscore_length_code.txt > ./code/v3_zscore_length_code_4.txt
awk '{if(NR>220*4 && NR<=(220*5)){print $0}}' v3_zscore_length_code.txt > ./code/v3_zscore_length_code_5.txt
awk '{if(NR>220*5 && NR<=(220*6)){print $0}}' v3_zscore_length_code.txt > ./code/v3_zscore_length_code_6.txt
awk '{if(NR>220*6 && NR<=(220*7)){print $0}}' v3_zscore_length_code.txt > ./code/v3_zscore_length_code_7.txt
awk '{if(NR>220*7 && NR<=(220*8)){print $0}}' v3_zscore_length_code.txt > ./code/v3_zscore_length_code_8.txt
awk '{if(NR>220*8 && NR<=(220*9)){print $0}}' v3_zscore_length_code.txt > ./code/v3_zscore_length_code_9.txt

nohup bash ./code/v3_zscore_length_code_1.txt > ./code/nohup.v3_zscore_length_code_1.out 2>&1&
nohup bash ./code/v3_zscore_length_code_2.txt > ./code/nohup.v3_zscore_length_code_2.out 2>&1&
nohup bash ./code/v3_zscore_length_code_3.txt > ./code/nohup.v3_zscore_length_code_3.out 2>&1&
nohup bash ./code/v3_zscore_length_code_4.txt > ./code/nohup.v3_zscore_length_code_4.out 2>&1&
nohup bash ./code/v3_zscore_length_code_5.txt > ./code/nohup.v3_zscore_length_code_5.out 2>&1&
nohup bash ./code/v3_zscore_length_code_6.txt > ./code/nohup.v3_zscore_length_code_6.out 2>&1&
nohup bash ./code/v3_zscore_length_code_7.txt > ./code/nohup.v3_zscore_length_code_7.out 2>&1&
nohup bash ./code/v3_zscore_length_code_8.txt > ./code/nohup.v3_zscore_length_code_8.out 2>&1&
nohup bash ./code/v3_zscore_length_code_9.txt > ./code/nohup.v3_zscore_length_code_9.out 2>&1&


# 2. get all files
cp DHS_v3_150-350bp.parafile DHS_v3_150-350bp_nor.parafile
vim DHS_v3_150-350bp_nor.parafile
cp DHS_v3_150-350bp.sh DHS_v3_150-350bp_nor.sh
vim DHS_v3_150-350bp_nor.sh

# 3. run IDEAS
nohup bash DHS_v3_150-350bp_nor.sh > ./nohup.DHS_v3_150-350bp_nor.out 2>&1&



#--------------------------------------------------------------------
# 4. rerun all using rep1
## 1) dhs_bins
# .input file
cp DHS_v3_100-400bp.input DHS_v3_100-400bp_oldATAC.input
awk '{FS=OFS=" "}{print $1,$2,"/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"$1"_"$2"_dhs.txt"}' \
DHS_v3_100-400bp_oldATAC.input > DHS_v3_100-400bp.input
# bed file
cut -f 1 /data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/embryonic-facial-prominence_15.5_H3K9ac_dhs.tab > tmp.txt
awk '{FS=OFS=" "}{if(NR==FNR){a[$4]=$0}else{print a[$1]}}' mm10_OCR-center_bins_v3_signal_based_space.bed tmp.txt > tmp2.txt
mv tmp2.txt mm10_OCR-center_bins_v3_signal_based_space.bed
#
nohup bash DHS_v3_100-400bp.sh > ./nohup.DHS_v3_100-400bp_Jan17.out 2>&1&


## 2) normal_bins
mkdir /data/zusers/fankaili/ideas/dhs_bins/normal_bins/
cd /data/zusers/fankaili/ideas/dhs_bins/normal_bins/
# .input file
awk '{FS=OFS=" "}{print $1,$2,"/data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/"$1"_"$2"_normal.txt"}' \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_v3_100-400bp_oldATAC.input > \
66samples_10marks_normal_bins.input
# .parafile
cp /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_v3_100-400bp.parafile 66samples_10marks_normal_bins.parafile
vim 66samples_10marks_normal_bins.parafile
# .sh file
cp /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_v3_100-400bp.sh 66samples_10marks_normal_bins.sh
vim 66samples_10marks_normal_bins.sh
#
cp /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/mm10.bed ./
cp -r /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/bin ./
cp -r /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/data ./
#
nohup bash 66samples_10marks_normal_bins.sh > ./nohup.66samples_10marks_normal_bins.out 2>&1&
