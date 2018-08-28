#!/bin/bash

# -- Kaili
# This script is for making matrix for ubi-rDHS.

cd /data/zusers/fankaili/ccre/tf/tsne/tsne_table/

##################
# 10,921 ubi-rDHS

awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5}else{print $1,a[$1]}}' /data/zusers/fankaili/ccre/tf/tsne/hg19_ubi-rDHS_distance2TSS_log10.txt \
/data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt > hg19_ubi-rDHS_matrix.txt
sed -i '1iid\tdistance' hg19_ubi-rDHS_matrix.txt
#
cut -f 3 /data/zusers/fankaili/ccre/hg19_ubi_ccRE_ccreid_agnostic.txt > temp.txt
sed -i '1iagnostic' temp.txt
paste hg19_ubi-rDHS_matrix.txt temp.txt > temp_2.txt
mv temp_2.txt hg19_ubi-rDHS_matrix.txt
#
while read line
do
    echo ${line} ;
    echo ${line}_DNase > temp_DNase.txt ;
    awk '{if(NR>1){print $3}}' /data/zusers/fankaili/ccre/tf/tsne/matrix/hg19_ubi-rDHS_${line}_five_dimension_matrix.txt >> temp_DNase.txt ;
    echo ${line}_H3K4me3 > temp_H3K4me3.txt ;
    awk '{if(NR>1){print $4}}' /data/zusers/fankaili/ccre/tf/tsne/matrix/hg19_ubi-rDHS_${line}_five_dimension_matrix.txt >> temp_H3K4me3.txt ;
    echo ${line}_H3K27ac > temp_H3K27ac.txt ;
    awk '{if(NR>1){print $5}}' /data/zusers/fankaili/ccre/tf/tsne/matrix/hg19_ubi-rDHS_${line}_five_dimension_matrix.txt >> temp_H3K27ac.txt ;
    echo ${line}_CTCF > temp_CTCF.txt ;
    awk '{if(NR>1){print $6}}' /data/zusers/fankaili/ccre/tf/tsne/matrix/hg19_ubi-rDHS_${line}_five_dimension_matrix.txt >> temp_CTCF.txt ;
    echo ${line}_group > temp_group.txt ;
    awk '{if(NR>1){print $7}}' /data/zusers/fankaili/ccre/tf/tsne/matrix/hg19_ubi-rDHS_${line}_five_dimension_matrix.txt >> temp_group.txt ;
    #
    paste hg19_ubi-rDHS_matrix.txt temp_DNase.txt temp_H3K4me3.txt temp_H3K27ac.txt temp_CTCF.txt temp_group.txt > temp_3.txt ;
    mv temp_3.txt hg19_ubi-rDHS_matrix.txt ;
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt
#
echo -e "GM_tsne_1\tGM_tsne_2" > temp_tsne.txt
cat hg19_ubi-rDHS_GM12878_tNSE.txt >> temp_tsne.txt
paste hg19_ubi-rDHS_matrix.txt temp_tsne.txt > temp_4.txt
mv temp_4.txt hg19_ubi-rDHS_matrix.txt
#
rm temp*.txt


##################
# 29,733 ubi-rDHS

awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5}else{print $1,a[$1]}}' /data/zusers/fankaili/ccre/tf/tsne/hg19_ubi-rDHS_2_distance2TSS_log10.txt \
/data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list_ccreid.txt > hg19_ubi-rDHS_2_matrix.txt
sed -i '1iid\tdistance' hg19_ubi-rDHS_2_matrix.txt
#
cut -f 3 /data/zusers/fankaili/ccre/hg19_ubi_ccRE_2_ccreid_agnostic.txt > temp.txt
sed -i '1iagnostic' temp.txt
paste hg19_ubi-rDHS_2_matrix.txt temp.txt > temp_2.txt
mv temp_2.txt hg19_ubi-rDHS_2_matrix.txt
#
while read line
do
    echo ${line} ;
    echo ${line}_DNase > temp_DNase.txt ;
    awk '{if(NR>1){print $3}}' /data/zusers/fankaili/ccre/tf/tsne/matrix_2/hg19_ubi-rDHS_2_${line}_five_dimension_matrix.txt >> temp_DNase.txt ;
    echo ${line}_H3K4me3 > temp_H3K4me3.txt ;
    awk '{if(NR>1){print $4}}' /data/zusers/fankaili/ccre/tf/tsne/matrix_2/hg19_ubi-rDHS_2_${line}_five_dimension_matrix.txt >> temp_H3K4me3.txt ;
    echo ${line}_H3K27ac > temp_H3K27ac.txt ;
    awk '{if(NR>1){print $5}}' /data/zusers/fankaili/ccre/tf/tsne/matrix_2/hg19_ubi-rDHS_2_${line}_five_dimension_matrix.txt >> temp_H3K27ac.txt ;
    echo ${line}_CTCF > temp_CTCF.txt ;
    awk '{if(NR>1){print $6}}' /data/zusers/fankaili/ccre/tf/tsne/matrix_2/hg19_ubi-rDHS_2_${line}_five_dimension_matrix.txt >> temp_CTCF.txt ;
    echo ${line}_group > temp_group.txt ;
    awk '{if(NR>1){print $7}}' /data/zusers/fankaili/ccre/tf/tsne/matrix_2/hg19_ubi-rDHS_2_${line}_five_dimension_matrix.txt >> temp_group.txt ;
    #
    paste hg19_ubi-rDHS_2_matrix.txt temp_DNase.txt temp_H3K4me3.txt temp_H3K27ac.txt temp_CTCF.txt temp_group.txt > temp_3.txt ;
    mv temp_3.txt hg19_ubi-rDHS_2_matrix.txt ;
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt
#
echo -e "GM_tsne_1\tGM_tsne_2" > temp_tsne.txt
cat hg19_ubi-rDHS_2_GM12878_tNSE.txt >> temp_tsne.txt
paste hg19_ubi-rDHS_2_matrix.txt temp_tsne.txt > temp_4.txt
mv temp_4.txt hg19_ubi-rDHS_2_matrix.txt
#
rm temp*.txt
