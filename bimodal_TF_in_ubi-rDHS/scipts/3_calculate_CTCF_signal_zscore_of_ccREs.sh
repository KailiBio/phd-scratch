#!/bin/bash

# -- Kaili
# This is the script for getting CTCF signal&zscore of ccREs in all biosamples.
# INPUT: hg19_cellline_with_CTCF.txt from Jun01.ccREs_TF.sh
#        (all celllines with CTCF data in hg19. only cell line here, tissue removed.)
# OUTPUT: 1) CTCF signal of ccREs in all celllines: /data/zusers/fankaili/ccre/tf/signal/
#         2) CTCF zscore of rDHS in all celllines: /data/zusers/fankaili/ccre/tf/zscore_ctcf/
#         3) Matrix

scriptDir="/data/zusers/fankaili/github/KailiBio/bimodal_TF_in_ubi-rDHS/scripts/"
outDir="/data/zusers/fankaili/ccre/tf/zscore_ctcf/"

# 1. calculate the CTCF signal
for cellline in `cat /data/zusers/fankaili/ccre/tf/hg19_cellline_with_CTCF.txt`
do
    echo ${cellline} ;
    python ${scriptDir}calculate_CTCF_signal_of_ccREs.py ${cellline} ;
done


# 2. calculate CTCF z-score
for cellline in `cat /data/zusers/fankaili/ccre/tf/encode_tf_file_list/encode_hg19_tf_cellline_with_cell_type_specific_list.txt`
do
    echo ${cellline} ;
    line=`grep "CTCF" /data/zusers/fankaili/ccre/tf/encode_tf_file_list/encode_hg19_${cellline}_tf_id_list.txt` ;
    id=`awk '{print $1}' <<< $line`;
    file_id=`awk '{print $2}' <<< $line`;
    bigWigAverageOverBed /data/projects/encode/data/${id}/${file_id}.bigWig \
    /data/zusers/moorej3/Registry-of-ccREs/hg19/V4/hg19-rDHSs.bed ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_CTCF_signal.tab ;
    python /data/zusers/fankaili/ccre/tf/zscore-normalization.py ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_CTCF_signal.tab > ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_CTCF_signal_zscore.txt ;
done


# 3. make the matrix
