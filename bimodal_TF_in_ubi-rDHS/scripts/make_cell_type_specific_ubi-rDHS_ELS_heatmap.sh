#!/bin/bash

# -- Kaili
# This script is for making heatmap for cell-type specific ubi-rDHS ELS.


# 1. get cell-type specific ubi-rDHS ELS matrix
cd /data/zusers/fankaili/ccre/tf/matrix/

echo "id" > hg19_cell_type_specific_ubi-rDHS_ELS_matrix.txt
cat /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt >> hg19_cell_type_specific_ubi-rDHS_ELS_matrix.txt

for cellline in `cat /data/zusers/fankaili/ccre/tf/hg19_cellline_with_CTCF_using.txt`
do
    if [ -s /data/zusers/fankaili/ccre/tf/cell_type_specific/ELS/ubi-rDHS_${cellline}_ELS_list.txt ]; then
        echo $cellline ;
        awk '{print $0}' /data/zusers/fankaili/ccre/tf/cell_type_specific/ELS/ubi-rDHS_${cellline}_ELS_list.txt | sed '1i'$cellline > temp_${cellline}.txt ;
        awk -v cellline="$cellline" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0,1}else{if($1=="id"){print $0,cellline}else{print $0,0}}}}' \
        temp_${cellline}.txt hg19_cell_type_specific_ubi-rDHS_ELS_matrix.txt > temp_matrix.txt ;
        mv temp_matrix.txt hg19_cell_type_specific_ubi-rDHS_ELS_matrix.txt ;
        rm temp_${cellline}.txt ;
    fi
done

## sort file
head -1 hg19_cell_type_specific_ubi-rDHS_ELS_matrix.txt > hg19_cell_type_specific_ubi-rDHS_ELS_matrix_ccreid.txt
awk 'NR>1' hg19_cell_type_specific_ubi-rDHS_ELS_matrix.txt | sort -k1n >> hg19_cell_type_specific_ubi-rDHS_ELS_matrix_ccreid.txt

# 2. make heatmap
Rscript make_cell_type_specific_ubi-rDHS_ELS_heatmap.R
