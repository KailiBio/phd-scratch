#!/bin/bash

# -- Kaili
# This script is for getting bimodal matrix of SMC3 and RAD21.
# INPUT:
# OUTPUT:

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scipts/"

#####################
## SMC3

# 1. calculate zscore
mkdir zscore_smc3
outDir="/data/zusers/fankaili/ccre/tf/zscore_smc3/"
#
for cellline in `cat /data/zusers/fankaili/ccre/tf/hg19_cellline_with_CTCF_using.txt`
do
    echo ${cellline} ;
    line=`grep "SMC3" /data/zusers/fankaili/ccre/tf/encode_tf_file_list/encode_hg19_${cellline}_tf_id_list.txt` ;
    if [ "$line" != "" ]; then
        n=`echo $line | awk '{print NF}'` ;
        nn=`expr $(($n/3))` ;
        for ((i=1;i<=$nn;i++))
        do
            a=`expr $((1+($i-1)*3))`;
            b=`expr $((2+($i-1)*3))`;
            id=`awk -v i="$a" '{print $i}' <<< $line`;
            file_id=`awk -v i="$b" '{print $i}' <<< $line`;
            bigWigAverageOverBed /data/projects/encode/data/${id}/${file_id}.bigWig \
            /data/zusers/moorej3/ENCODE-Registry/hg19/V4/hg19-rDHSs.bed ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_SMC3_signal.tab ;
            python /data/zusers/fankaili/ccre/tf/zscore-normalization.py ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_SMC3_signal.tab > ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_SMC3_signal_zscore.txt ;
        done
    fi
done

# 2. make zscore matrix
bash ${scriptDir}make_tf_zscore_matrix.sh -l /data/zusers/moorej3/ENCODE-Registry/hg19/V4/hg19-rDHSs.bed -d /data/zusers/fankaili/ccre/tf/matrix/ \
 -m hg19_rDHS_SMC3_zscore_matrix.txt -t SMC3 -f /data/zusers/fankaili/ccre/tf/zscore_smc3/



## get ubi-rDHS signal matrix
awk '{FS=OFS="\t"}{if($1=="id"){print $0}}' hg19_rDHS_CTCF_zscore_matrix.txt > hg19_ubi-rDHS_CTCF_zscore_matrix_DHSID.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt hg19_rDHS_CTCF_zscore_matrix.txt >> hg19_ubi-rDHS_CTCF_zscore_matrix_DHSID.txt
awk '{FS=OFS="\t"}{if($1=="id"){print $0}}' hg19_ubi-rDHS_CTCF_zscore_matrix_DHSID.txt > hg19_ubi-rDHS_CTCF_zscore_matrix.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_CTCF_zscore_matrix_DHSID.txt >> hg19_ubi-rDHS_CTCF_zscore_matrix.txt










#####################
## RAD21

# 1. calculate zscore
mkdir zscore_rad21
outDir="/data/zusers/fankaili/ccre/tf/zscore_rad21/"
#
for cellline in `cat /data/zusers/fankaili/ccre/tf/hg19_cellline_with_CTCF_using.txt`
do
    echo ${cellline} ;
    line=`grep "RAD21" /data/zusers/fankaili/ccre/tf/encode_tf_file_list/encode_hg19_${cellline}_tf_id_list.txt` ;
    if [ "$line" != "" ]; then
        n=`echo $line | awk '{print NF}'` ;
        nn=`expr $(($n/3))` ;
        for ((i=1;i<=$nn;i++))
        do
            a=`expr $((1+($i-1)*3))`;
            b=`expr $((2+($i-1)*3))`;
            id=`awk -v i="$a" '{print $i}' <<< $line`;
            file_id=`awk -v i="$b" '{print $i}' <<< $line`;
            bigWigAverageOverBed /data/projects/encode/data/${id}/${file_id}.bigWig \
            /data/zusers/moorej3/ENCODE-Registry/hg19/V4/hg19-rDHSs.bed ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_RAD21_signal.tab ;
            python /data/zusers/fankaili/ccre/tf/zscore-normalization.py ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_RAD21_signal.tab > ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_RAD21_signal_zscore.txt ;
        done
    fi
done
