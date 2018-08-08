#!/bin/bash

# -- Kaili
# This is the script for getting CTCF signal&zscore of ccREs in all biosamples.
# INPUT: hg19_cellline_with_CTCF.txt from Jun01.ccREs_TF.sh
#        (all celllines with CTCF data in hg19. only cell line here, tissue removed.)
# OUTPUT: 1) CTCF signal of ccREs in all celllines: /data/zusers/fankaili/ccre/tf/signal/
#         2) CTCF zscore of rDHS in all celllines: /data/zusers/fankaili/ccre/tf/zscore_ctcf/
#         3) Matrix: signal/zscore matrix for ccre/rDHS and ubi-rDHS.


scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scipts/"

# 1. calculate the CTCF signal
for cellline in `cat /data/zusers/fankaili/ccre/tf/encode_tf_file_list/encode_hg19_tf_cellline_with_cell_type_specific_list.txt`
do
    echo ${cellline} ;
    #cellline=NT2-D1
    python ${scriptDir}calculate_CTCF_signal_of_ccREs.py ${cellline} ;
done


# 2. calculate CTCF z-score
outDir="/data/zusers/fankaili/ccre/tf/zscore_ctcf/"
#
for cellline in `cat /data/zusers/fankaili/ccre/tf/encode_tf_file_list/encode_hg19_tf_cellline_with_cell_type_specific_list.txt`
do
    echo ${cellline} ;
    line=`grep "CTCF" /data/zusers/fankaili/ccre/tf/encode_tf_file_list/encode_hg19_${cellline}_tf_id_list.txt` ;
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
            /data/zusers/moorej3/Registry-of-ccREs/hg19/V4/hg19-rDHSs.bed ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_CTCF_signal.tab ;
            python /data/zusers/fankaili/ccre/tf/zscore-normalization.py ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_CTCF_signal.tab > ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_CTCF_signal_zscore.txt ;
        done
    fi
done

# 3. make the matrix
cd /data/zusers/fankaili/ccre/tf/matrix/

## 1) ubi-rDHS signal
### get ccREs signal matrix
echo "id" > hg19_ccRE_CTCF_signal_matrix.txt
awk '{print $4}' /data/zusers/fankaili/ccre/hg19-cREs.bed >> hg19_ccRE_CTCF_signal_matrix.txt
#
for file in `ls /data/zusers/fankaili/ccre/tf/signal/*_final.txt`
do
    filename0=${file%_final.txt} ;
    filename=${filename0#/data/zusers/fankaili/ccre/tf/signal/CTCF_signal_ccREs_} ;
    echo $filename ;
    echo -e "id\t"$filename > temp.txt ;
    cat $file >> temp.txt ;
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print $0,b[$1]}}}' temp.txt hg19_ccRE_CTCF_signal_matrix.txt > temp_matrix.txt ;
    mv temp_matrix.txt hg19_ccRE_CTCF_signal_matrix.txt ;
    rm temp*.txt ;
done

### get ubi-rDHS signal matrix
awk '{FS=OFS="\t"}{if($1=="id"){print $0}}' hg19_ccRE_CTCF_signal_matrix.txt > hg19_ubi-rDHS_CTCF_signal_matrix.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt hg19_ccRE_CTCF_signal_matrix.txt >> hg19_ubi-rDHS_CTCF_signal_matrix.txt


## 2) ubi-rDHS z-score
### get rDHS z-score matrix
echo "id" > hg19_rDHS_CTCF_zscore_matrix.txt
awk '{print $4}' /data/zusers/moorej3/Registry-of-ccREs/hg19/V4/hg19-rDHSs.bed >> hg19_rDHS_CTCF_zscore_matrix.txt
#
for file in `ls /data/zusers/fankaili/ccre/tf/zscore_ctcf/*_CTCF_signal_zscore.txt`
do
    filename0=${file%_CTCF_signal_zscore.txt} ;
    filename=${filename0#/data/zusers/fankaili/ccre/tf/zscore_ctcf/hg19_rDHS_} ;
    echo $filename ;
    echo -e "id\t"$filename > temp.txt ;
    awk '{FS=" \t";OFS="\t"}{print $1,$2}' $file >> temp.txt ;
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print $0,b[$1]}}}' temp.txt hg19_rDHS_CTCF_zscore_matrix.txt > temp_matrix.txt ;
    mv temp_matrix.txt hg19_rDHS_CTCF_zscore_matrix.txt ;
    rm temp*.txt;
done

### get ubi-rDHS signal matrix
awk '{FS=OFS="\t"}{if($1=="id"){print $0}}' hg19_rDHS_CTCF_zscore_matrix.txt > hg19_ubi-rDHS_CTCF_zscore_matrix_DHSID.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt hg19_rDHS_CTCF_zscore_matrix.txt >> hg19_ubi-rDHS_CTCF_zscore_matrix_DHSID.txt
awk '{FS=OFS="\t"}{if($1=="id"){print $0}}' hg19_ubi-rDHS_CTCF_zscore_matrix_DHSID.txt > hg19_ubi-rDHS_CTCF_zscore_matrix.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){b=a[$1];$1=b;print $0}}}' /data/zusers/fankaili/ccre/ccREs_ID_transfer_clean.bed hg19_ubi-rDHS_CTCF_zscore_matrix_DHSID.txt >> hg19_ubi-rDHS_CTCF_zscore_matrix.txt
