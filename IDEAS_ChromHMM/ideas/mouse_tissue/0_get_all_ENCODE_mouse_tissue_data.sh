#!/bin/bash

# -- Kaili
# This is script for getting all the 870 mouse tissue file ID from ENCODE.
# 1. try to get the expID of all mouse data
###### rep1
# 2. get all the ENCODE mouse rep1 data, for rerunning IDEAS.
## 1) ATAC-seq
## 2) make DNAme bigWig in 1bp resolution
## 3) get 8 histone mark ChIP-seq
# 3. get rep1 signal value for all bins
###### rep2
# 4. get all the ENCODE mouse rep2 data, for rerunning IDEAS.
# 5. get rep2 signal value for all bins

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/"

cd ${workDir}

# 1. try to get the expID of all mouse data
# !. need double check if you use this file.
python ${scriptDir}get_mouse_data_from_ENCODE_json.py


# 2. get all the ENCODE mouse rep1 data, for rerunning IDEAS.
## 1) ATAC-seq
python ${scriptDir}get_ENCODE_mouse_rep1_ATAC_data.py
## 2) make DNAme bigWig in 1bp resolution
# fileID are CpG bedMethyl file
python ${scriptDir}get_ENCODE_mouse_rep1_DNAme_data.py
## 3) get 8 histone mark ChIP-seq
python ${scriptDir}get_ENCODE_mouse_rep1_8HM_data.py
## 4) get CTCF data
python ${scriptDir}get_ENCODE_mouse_rep1_CTCF_data.py
##
cat ENCODE_mouse_rep1_8HM_filelist.txt ENCODE_mouse_rep1_ATAC_filelist.txt ENCODE_mouse_rep1_CTCF_filelist.txt ENCODE_mouse_rep1_DNAme_filelist.txt \
| sort -k1,1 -k2,2 > ENCODE_mouse_rep1_filelist.txt

# 3. get rep1 signal value for all bins
#
normal_bins_signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/"
dhs_bins_signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
mkdir ${normal_bins_signal_path}
mkdir ${dhs_bins_signal_path}
#
encode_data_path="/data/projects/encode/data/"
#
cat ENCODE_mouse_rep1_8HM_filelist.txt ENCODE_mouse_rep1_ATAC_filelist.txt \
ENCODE_mouse_rep1_CTCF_filelist.txt > tmp_ENCODE_mouse_rep1_filelist.txt
#
if [ -f ENCODE_rep1_signal_code.sh ]
then
    rm ENCODE_rep1_signal_code.sh
fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "cp "${encode_data_path}${expID}"/"${fileID}".bigWig /tmp/" >> ENCODE_rep1_signal_code.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/normal_bins/mm10_tab_noM.bed /tmp/"${sample}"_"${assay}"_normal.tab" >> ENCODE_rep1_signal_code.sh
    echo "awk '{print ""\$""5}' /tmp/"${sample}"_"${assay}"_normal.tab > /tmp/"${sample}"_"${assay}"_normal.txt" >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}"_normal.tab" ${normal_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}"_normal.txt" ${normal_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3.sorted.bed /tmp/"${sample}"_"${assay}"_dhs.tab" >> ENCODE_rep1_signal_code.sh
    echo "awk '{print ""\$""5}' /tmp/"${sample}"_"${assay}"_dhs.tab > /tmp/"${sample}"_"${assay}"_dhs.txt" >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}"_dhs.tab" ${dhs_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}"_dhs.txt" ${dhs_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    #
    echo "rm /tmp/"${fileID}".bigWig" >> ENCODE_rep1_signal_code.sh
done < tmp_ENCODE_mouse_rep1_filelist.txt
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "gzip -c /data/projects/encode/data/"${expID}"/"${fileID}".bed.gz > /tmp/"${fileID}".bed" >> ENCODE_rep1_signal_code.sh
    echo "bash /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/aggregation_rep1/bedMethyl2bigWig.sh /tmp/"${fileID}".bed /tmp/mm10.chrom.sizes /tmp/"${fileID}".bigWig" >> ENCODE_rep1_signal_code.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/run_ideas_p_value/mm10_tab.bed /tmp/"${sample}"_"${assay}"_normal.tab" >> ENCODE_rep1_signal_code.sh
    echo "awk '{print ""\$""6}' /tmp/"${sample}"_"${assay}"_normal.tab > /tmp/"${sample}"_"${assay}"_normal.txt" >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}"_normal.tab" ${normal_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}"_normal.txt" ${normal_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3.sorted.bed /tmp/"${sample}"_"${assay}"_dhs.tab" >> ENCODE_rep1_signal_code.sh
    echo "awk '{print ""\$""6}' /tmp/"${sample}"_"${assay}"_dhs.tab > /tmp/"${sample}"_"${assay}"_dhs.txt" >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}"_dhs.tab" ${dhs_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}"_dhs.txt" ${dhs_bins_signal_path} >> ENCODE_rep1_signal_code.sh
    #
    echo "rm /tmp/"${fileID}".bigWig /tmp/"${fileID}".bed" >> ENCODE_rep1_signal_code.sh
done < ENCODE_mouse_rep1_DNAme_filelist.txt
##### cut code into 10 samples each
mkdir /data/zusers/fankaili/ideas/code/get_rep1_signal/
#
for i in {1..68}
do
    awk -v i="$i" '{if(NR>((i-1)*100) && NR<=(100*i)){print $0}}' ENCODE_rep1_signal_code.sh > \
    /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_code_${i}.sh
done

##########################
# rerun DNAme for new bigWig

# normal_bins_signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/"
# dhs_bins_signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
# rm ENCODE_rep1_signal_code2.sh
# while read line
# do
#     sample=`awk '{print $1}' <<< ${line}`
#     assay=`awk '{print $2}' <<< ${line}`
#     expID=`awk '{print $4}' <<< ${line}`
#     fileID=`awk '{print $5}' <<< ${line}`
#     #
#     echo "bigWigAverageOverBed /data/zusers/fankaili/ideas/signal/mouse_66samples_DNAme_bigWig/"${fileID}".bigWig /data/zusers/fankaili/ideas/run_ideas_p_value/mm10_tab.bed "${normal_bins_signal_path}${sample}"_"${assay}"_normal.tab" >> ENCODE_rep1_signal_code2.sh
#     echo "awk '{print ""\$""6}' "${normal_bins_signal_path}${sample}"_"${assay}"_normal.tab > "${normal_bins_signal_path}${sample}"_"${assay}"_normal.txt" >> ENCODE_rep1_signal_code2.sh
#     #
#     echo "bigWigAverageOverBed /data/zusers/fankaili/ideas/signal/mouse_66samples_DNAme_bigWig/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3.sorted.bed "${dhs_bins_signal_path}${sample}"_"${assay}"_dhs.tab" >> ENCODE_rep1_signal_code2.sh
#     echo "awk '{print ""\$""6}' "${dhs_bins_signal_path}${sample}"_"${assay}"_dhs.tab > "${dhs_bins_signal_path}${sample}"_"${assay}"_dhs.txt" >> ENCODE_rep1_signal_code2.sh
# done < ENCODE_mouse_rep1_DNAme_filelist.txt
#
# for i in {1..11}
# do
#     awk -v i="$i" '{if(NR>((i-1)*24) && NR<=(24*i)){print $0}}' ENCODE_rep1_signal_code2.sh > \
#     /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_code2_${i}.sh
#     nohup bash /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_code2_${i}.sh \
#     > /data/zusers/fankaili/ideas/code/get_rep1_signal/nohup.ENCODE_rep1_signal_code2_${i}.out 2>&1&
# done


###########################
# rerun normal bins with new bed without chrM
# rm chrM
awk '{FS=" ";OFS="\t"}{if($1!="chrM"){print $1,$2,$3,$4}}' mm10.bed > mm10_tab_noM.bed
awk '{FS=OFS=" "}{if($1!="chrM"){print $1,$2,$3,$4}}' mm10.bed > tmp.bed
mv tmp.bed mm10.bed
#
normal_bins_signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/"
dhs_bins_signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
mkdir ${normal_bins_signal_path}
mkdir ${dhs_bins_signal_path}
#
encode_data_path="/data/projects/encode/data/"
#
if [ -f ENCODE_rep1_signal_code3.sh ]; then rm ENCODE_rep1_signal_code3.sh; fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "bigWigAverageOverBed "${encode_data_path}${expID}"/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/normal_bins/mm10_tab_noM.bed "${normal_bins_signal_path}${sample}"_"${assay}"_normal.tab" >> ENCODE_rep1_signal_code3.sh
    echo "awk '{print ""\$""5}' "${normal_bins_signal_path}${sample}"_"${assay}"_normal.tab > "${normal_bins_signal_path}${sample}"_"${assay}"_normal.txt" >> ENCODE_rep1_signal_code3.sh
done < tmp_ENCODE_mouse_rep1_filelist.txt
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "bigWigAverageOverBed /data/zusers/fankaili/ideas/signal/mouse_66samples_DNAme_bigWig/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/normal_bins/mm10_tab_noM.bed "${normal_bins_signal_path}${sample}"_"${assay}"_normal.tab" >> ENCODE_rep1_signal_code3.sh
    echo "awk '{print ""\$""6}' "${normal_bins_signal_path}${sample}"_"${assay}"_normal.tab > "${normal_bins_signal_path}${sample}"_"${assay}"_normal.txt" >> ENCODE_rep1_signal_code3.sh
done < ENCODE_mouse_rep1_DNAme_filelist.txt
###
for i in {1..14}
do
    awk -v i="$i" '{if(NR>((i-1)*100) && NR<=(100*i)){print $0}}' ENCODE_rep1_signal_code3.sh > \
    /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_code3_${i}.sh
    nohup bash /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_code3_${i}.sh > \
    /data/zusers/fankaili/ideas/code/get_rep1_signal/nohup.ENCODE_rep1_signal_code3_${i}.out 2>&1&
done


# 4. get all the ENCODE mouse rep2 data, for rerunning IDEAS.
## 1) ATAC-seq
python ${scriptDir}get_ENCODE_mouse_rep2_ATAC_data.py
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$1}else{print a[$3],$0}}' ENCODE_mouse_rep1_ATAC_filelist.txt ENCODE_mouse_rep2_ATAC_filelist.txt > tmp.ENCODE_mouse_rep2_ATAC_filelist.txt
mv tmp.ENCODE_mouse_rep2_ATAC_filelist.txt ENCODE_mouse_rep2_ATAC_filelist.txt
## 2) make DNAme bigWig in 1bp resolution
# fileID are CpG bedMethyl file
python ${scriptDir}get_ENCODE_mouse_rep2_DNAme_data.py
## 3) get 8 histone mark ChIP-seq
python ${scriptDir}get_ENCODE_mouse_rep2_8HM_data.py
## 4) get CTCF data
python ${scriptDir}get_ENCODE_mouse_rep2_CTCF_data.py
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$1}else{print a[$3],$0}}' ENCODE_mouse_rep1_CTCF_filelist.txt ENCODE_mouse_rep2_CTCF_filelist.txt > tmp.ENCODE_mouse_rep2_CTCF_filelist.txt
mv tmp.ENCODE_mouse_rep2_CTCF_filelist.txt ENCODE_mouse_rep2_CTCF_filelist.txt
##
cat ENCODE_mouse_rep2_8HM_filelist.txt ENCODE_mouse_rep2_ATAC_filelist.txt ENCODE_mouse_rep2_CTCF_filelist.txt ENCODE_mouse_rep2_DNAme_filelist.txt \
| sort -k1,1 -k2,2 > ENCODE_mouse_rep2_filelist.txt



# 5. get rep2 signal value for all bins
## 1) dhs-bins
dhs_bins_signal_path="/data/zusers/fankaili/ideas/signal/rep2_signal_dhs_bins/"
mkdir ${dhs_bins_signal_path}
#
encode_data_path="/data/projects/encode/data/"
#
cat ENCODE_mouse_rep2_8HM_filelist.txt ENCODE_mouse_rep2_ATAC_filelist.txt \
ENCODE_mouse_rep2_CTCF_filelist.txt > tmp_ENCODE_mouse_rep2_filelist.txt
#
if [ -f ENCODE_rep2_signal_code.sh ]
then
    rm ENCODE_rep2_signal_code.sh
fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "cp "${encode_data_path}${expID}"/"${fileID}".bigWig /tmp/" >> ENCODE_rep2_signal_code.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3.sorted.bed /tmp/"${sample}"_"${assay}"_dhs.tab" >> ENCODE_rep2_signal_code.sh
    echo "awk '{print ""\$""5}' /tmp/"${sample}"_"${assay}"_dhs.tab > /tmp/"${sample}"_"${assay}"_dhs.txt" >> ENCODE_rep2_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}"_dhs.tab" ${dhs_bins_signal_path} >> ENCODE_rep2_signal_code.sh
    echo "mv /tmp/"${sample}"_"${assay}"_dhs.txt" ${dhs_bins_signal_path} >> ENCODE_rep2_signal_code.sh
    #
    echo "rm /tmp/"${fileID}".bigWig" >> ENCODE_rep2_signal_code.sh
done < tmp_ENCODE_mouse_rep2_filelist.txt
#
if [ -f ENCODE_rep2_signal_code_DNAme.sh ]
then
    rm ENCODE_rep2_signal_code_DNAme.sh
fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "cat /data/projects/encode/data/"${expID}"/"${fileID}".bed.gz | gzip -d > /tmp/"${fileID}".bed" >> ENCODE_rep2_signal_code_DNAme.sh
    echo "bash /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/bedMethyl2bigWig_2.sh "${fileID}" /home/fankaili/mm10.chrom.sizes" >> ENCODE_rep2_signal_code_DNAme.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/mm10_OCR-center_bins_v3.sorted.bed /tmp/"${sample}"_"${assay}"_dhs.tab" >> ENCODE_rep2_signal_code_DNAme.sh
    echo "awk '{print ""\$""6}' /tmp/"${sample}"_"${assay}"_dhs.tab > /tmp/"${sample}"_"${assay}"_dhs.txt" >> ENCODE_rep2_signal_code_DNAme.sh
    echo "mv /tmp/"${sample}"_"${assay}"_dhs.tab" ${dhs_bins_signal_path} >> ENCODE_rep2_signal_code_DNAme.sh
    echo "mv /tmp/"${sample}"_"${assay}"_dhs.txt" ${dhs_bins_signal_path} >> ENCODE_rep2_signal_code_DNAme.sh
    #
    echo "rm /tmp/"${fileID}".bigWig /tmp/"${fileID}".bed" >> ENCODE_rep2_signal_code_DNAme.sh
done < ENCODE_mouse_rep2_DNAme_filelist.txt
##### cut code into 10 samples each
mkdir /data/zusers/fankaili/ideas/code/get_rep2_signal/
#
for i in {1..41}
do
    awk -v i="$i" '{if(NR>((i-1)*90) && NR<=(90*i)){print $0}}' ENCODE_rep2_signal_code.sh > \
    /data/zusers/fankaili/ideas/code/get_rep2_signal/ENCODE_rep2_signal_code_${i}.sh
done
#
for i in {1..11}
do
    awk -v i="$i" '{if(NR>((i-1)*42) && NR<=(42*i)){print $0}}' ENCODE_rep2_signal_code_DNAme.sh > \
    /data/zusers/fankaili/ideas/code/get_rep2_signal/ENCODE_rep2_signal_code_DNAme_${i}.sh
done
## 2) normal_bins
normal_bins_signal_path="/data/zusers/fankaili/ideas/signal/rep2_signal_normal_bins/"
mkdir ${normal_bins_signal_path}
#
encode_data_path="/data/projects/encode/data/"
#
if [ -f ENCODE_rep2_signal_code2.sh ]
then
    rm ENCODE_rep2_signal_code2.sh
fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "cp "${encode_data_path}${expID}"/"${fileID}".bigWig /tmp/" >> ENCODE_rep2_signal_code2.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/normal_bins/mm10_tab_noM.bed /tmp/"${sample}"_"${assay}"_normal.tab" >> ENCODE_rep2_signal_code2.sh
    echo "awk '{print ""\$""5}' /tmp/"${sample}"_"${assay}"_normal.tab > /tmp/"${sample}"_"${assay}"_normal.txt" >> ENCODE_rep2_signal_code2.sh
    echo "mv /tmp/"${sample}"_"${assay}"_normal.tab" ${normal_bins_signal_path} >> ENCODE_rep2_signal_code2.sh
    echo "mv /tmp/"${sample}"_"${assay}"_normal.txt" ${normal_bins_signal_path} >> ENCODE_rep2_signal_code2.sh
    #
    echo "rm /tmp/"${fileID}".bigWig" >> ENCODE_rep2_signal_code2.sh
done < tmp_ENCODE_mouse_rep2_filelist.txt
#
if [ -f ENCODE_rep2_signal_code_DNAme2.sh ]
then
    rm ENCODE_rep2_signal_code_DNAme2.sh
fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "cat /data/projects/encode/data/"${expID}"/"${fileID}".bed.gz | gzip -d > /tmp/"${fileID}".bed" >> ENCODE_rep2_signal_code_DNAme2.sh
    echo "bash /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/bedMethyl2bigWig_2.sh "${fileID}" /home/fankaili/mm10.chrom.sizes" >> ENCODE_rep2_signal_code_DNAme2.sh
    #
    echo "/home/fankaili/bigWigAverageOverBed /tmp/"${fileID}".bigWig /data/zusers/fankaili/ideas/dhs_bins/normal_bins/mm10_tab_noM.bed /tmp/"${sample}"_"${assay}"_normal.tab" >> ENCODE_rep2_signal_code_DNAme2.sh
    echo "awk '{print ""\$""6}' /tmp/"${sample}"_"${assay}"_normal.tab > /tmp/"${sample}"_"${assay}"_normal.txt" >> ENCODE_rep2_signal_code_DNAme2.sh
    echo "mv /tmp/"${sample}"_"${assay}"_normal.tab" ${normal_bins_signal_path} >> ENCODE_rep2_signal_code_DNAme2.sh
    echo "mv /tmp/"${sample}"_"${assay}"_normal.txt" ${normal_bins_signal_path} >> ENCODE_rep2_signal_code_DNAme2.sh
    #
    echo "rm /tmp/"${fileID}".bigWig /tmp/"${fileID}".bed" >> ENCODE_rep2_signal_code_DNAme2.sh
done < ENCODE_mouse_rep2_DNAme_filelist.txt
##### cut code into 10 samples each
for i in {1..41}
do
    awk -v i="$i" '{if(NR>((i-1)*90) && NR<=(90*i)){print $0}}' ENCODE_rep2_signal_code2.sh > \
    /data/zusers/fankaili/ideas/code/get_rep2_signal/ENCODE_rep2_signal_code2_${i}.sh
done
#
for i in {1..11}
do
    awk -v i="$i" '{if(NR>((i-1)*42) && NR<=(42*i)){print $0}}' ENCODE_rep2_signal_code_DNAme2.sh > \
    /data/zusers/fankaili/ideas/code/get_rep2_signal/ENCODE_rep2_signal_code_DNAme2_${i}.sh
done
