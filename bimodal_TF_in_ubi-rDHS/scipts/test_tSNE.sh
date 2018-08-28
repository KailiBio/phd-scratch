#!/bin/bash

# -- Kaili
# This script is for runing tSNE for ubi-rDHS/ubi-rDHS2 datasets.
# using five dimension data: log distance to TSS, DNase, H3K4me3, H3K27ac, CTCF

masterFile="/data/zusers/moorej3/ENCODE-Registry/hg19/V4/Cell-Type-Specific/Master-Cell-List.txt"
zscoreDir="/data/zusers/moorej3/ENCODE-Registry/hg19/V4/signal-output/"
workDir="/data/zusers/fankaili/ccre/tf/tsne/"

cd ${workDir}

# 1. get all the biosamples with all four dimension data
grep -v '\-\-\-'  ${masterFile}> /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_master.txt
### all 21 biosamples with total 4 dimension datas.

# 2. get distance to TSS
## 1) get ubi-rDHS bed file
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $1,$2,$3,$4}}}' /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt \
/data/zusers/fankaili/ccre/hg19-cREs.bed | sort -k1,1 -k2,2n > /data/zusers/fankaili/ccre/hg19_ubi-rDHS.bed

## 2) get distance, log10
cd ${workDir}
bedtools closest -a /data/zusers/fankaili/ccre/hg19_ubi-rDHS.bed -b \
/data/zusers/moorej3/moorej.ghpcc.project/Reference/Human/hg19/Gencode19/TSS.Filtered.bed -d -k 1 | \
awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$13}' | sort -u > hg19_ubi-rDHS_distance2TSS.txt
# get log10
awk '{FS=OFS="\t"}{print $1,$2,$3,$4,log($5+1)/log(10)}' hg19_ubi-rDHS_distance2TSS.txt > hg19_ubi-rDHS_distance2TSS_log10.txt

# 3. make matrix for each biosample (21 in total)
while read line
do
    echo ${line} ;
    dnaseExp=`awk '{FS=OFS="\t"}{print $1}' <<< ${line}` ;
    dnaseID=`awk '{FS=OFS="\t"}{print $2}' <<< ${line}` ;
    h3k4me3Exp=`awk '{FS=OFS="\t"}{print $3}' <<< ${line}` ;
    h3k4me3ID=`awk '{FS=OFS="\t"}{print $4}' <<< ${line}` ;
    h3k27acExp=`awk '{FS=OFS="\t"}{print $5}' <<< ${line}` ;
    h3k27acID=`awk '{FS=OFS="\t"}{print $6}' <<< ${line}` ;
    ctcfExp=`awk '{FS=OFS="\t"}{print $7}' <<< ${line}` ;
    ctcfID=`awk '{FS=OFS="\t"}{print $8}' <<< ${line}` ;
    cellline=`awk '{FS=OFS="\t"}{print $9}' <<< ${line}` ;
    # log10(distance)
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$4]=$5}else{if(a[$1]){print $1,b[$1]}}}' hg19_ubi-rDHS_distance2TSS_log10.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt > temp.distance.txt ;
    # DNase
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${dnaseExp}-${dnaseID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt > temp.dnase.txt ;
    # H3K4me3
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${h3k4me3Exp}-${h3k4me3ID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt > temp.h3k4me3.txt ;
    # H3K27ac
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${h3k27acExp}-${h3k27acID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt > temp.h3k27ac.txt ;
    # CTCF
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${ctcfExp}-${ctcfID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt > temp.ctcf.txt ;
    #
    paste temp.distance.txt temp.dnase.txt temp.h3k4me3.txt temp.h3k27ac.txt temp.ctcf.txt > hg19_ubi-rDHS_${cellline}_five_dimension_matrix.txt ;
    sed -i '1iid\tdistance\tdnase\th3k4me3\th3k27ac\tctcf' hg19_ubi-rDHS_${cellline}_five_dimension_matrix.txt ;
    rm temp.*.txt ;
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_master.txt

# 4. test using GM12878
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$10}else{if(a[$1]){print $0,a[$1]}}}' \
/data/zusers/moorej3/ENCODE-Registry/hg19/V4/Cell-Type-Specific/Five-Group/ENCFF131XEA_ENCFF818GNV_ENCFF180LKW_ENCFF886KRA.cREs.bed hg19_ubi-rDHS_GM12878_five_dimension_matrix.txt > hg19_ubi-rDHS_GM12878_five_dimension_matrix2.txt
sed -i '1iid\tdistance\tdnase\th3k4me3\th3k27ac\tctcf\ttype' hg19_ubi-rDHS_GM12878_five_dimension_matrix2.txt
#
Rscript test_tSNE_on_GM12878.R

##########################
## for 29,773 ubi-rDHSs

# 1. get ubi-rDHS2 list
cd /data/zusers/fankaili/ccre/
awk '{FS=OFS="\t"}{if($2>=346){print $1}}' /data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V4/Registry/V4-hg19/ccRE-DNase-Biosample-Counts.txt > ubi_rDHS_2_hg19_list.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$4}else{if(a[$1]){print a[$1]}}}' ccREs_ID_transfer_clean.bed ubi_rDHS_2_hg19_list.txt > ubi_ccREs_2_hg19_list_ccreid.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5}else{if(a[$1]){print a[$1]}}}' ccREs_ID_transfer_clean.bed ubi_ccREs_2_hg19_list_ccreid.txt > ubi_ccREs_2_hg19_list.txt

# 2. get distance to TSS
## 1) get ubi-rDHS bed file
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $1,$2,$3,$4}}}' /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list_ccreid.txt \
/data/zusers/fankaili/ccre/hg19-cREs.bed | sort -k1,1 -k2,2n > /data/zusers/fankaili/ccre/hg19_ubi-rDHS_2.bed

## 2) get distance, log10
cd ${workDir}
bedtools closest -a /data/zusers/fankaili/ccre/hg19_ubi-rDHS_2.bed -b \
/data/zusers/moorej3/moorej.ghpcc.project/Reference/Human/hg19/Gencode19/TSS.Filtered.bed -d -k 1 | \
awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$13}' | sort -u > hg19_ubi-rDHS_2_distance2TSS.txt
# get log10
awk '{FS=OFS="\t"}{print $1,$2,$3,$4,log($5+1)/log(10)}' hg19_ubi-rDHS_2_distance2TSS.txt > hg19_ubi-rDHS_2_distance2TSS_log10.txt


# 3. make matrix for each biosample (21 in total)
while read line
do
    echo ${line} ;
    dnaseExp=`awk '{FS=OFS="\t"}{print $1}' <<< ${line}` ;
    dnaseID=`awk '{FS=OFS="\t"}{print $2}' <<< ${line}` ;
    h3k4me3Exp=`awk '{FS=OFS="\t"}{print $3}' <<< ${line}` ;
    h3k4me3ID=`awk '{FS=OFS="\t"}{print $4}' <<< ${line}` ;
    h3k27acExp=`awk '{FS=OFS="\t"}{print $5}' <<< ${line}` ;
    h3k27acID=`awk '{FS=OFS="\t"}{print $6}' <<< ${line}` ;
    ctcfExp=`awk '{FS=OFS="\t"}{print $7}' <<< ${line}` ;
    ctcfID=`awk '{FS=OFS="\t"}{print $8}' <<< ${line}` ;
    cellline=`awk '{FS=OFS="\t"}{print $9}' <<< ${line}` ;
    # log10(distance)
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$4]=$5}else{if(a[$1]){print $1,b[$1]}}}' hg19_ubi-rDHS_2_distance2TSS_log10.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list_ccreid.txt > temp.distance.txt ;
    # DNase
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${dnaseExp}-${dnaseID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list.txt > temp.dnase.txt ;
    # H3K4me3
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${h3k4me3Exp}-${h3k4me3ID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list.txt> temp.h3k4me3.txt ;
    # H3K27ac
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${h3k27acExp}-${h3k27acID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list.txt > temp.h3k27ac.txt ;
    # CTCF
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${ctcfExp}-${ctcfID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list.txt > temp.ctcf.txt ;
    #
    paste temp.distance.txt temp.dnase.txt temp.h3k4me3.txt temp.h3k27ac.txt temp.ctcf.txt > hg19_ubi-rDHS_2_${cellline}_five_dimension_matrix.txt ;
    sed -i '1iid\tdistance\tdnase\th3k4me3\th3k27ac\tctcf' hg19_ubi-rDHS_2_${cellline}_five_dimension_matrix.txt ;
    rm temp.*.txt ;
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_master.txt

# 4. test using GM12878
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$10}else{if(a[$1]){print $0,a[$1]}}}' \
/data/zusers/moorej3/ENCODE-Registry/hg19/V4/Cell-Type-Specific/Five-Group/ENCFF131XEA_ENCFF818GNV_ENCFF180LKW_ENCFF886KRA.cREs.bed \
hg19_ubi-rDHS_2_GM12878_five_dimension_matrix.txt > hg19_ubi-rDHS_2_GM12878_five_dimension_matrix2.txt
sed -i '1iid\tdistance\tdnase\th3k4me3\th3k27ac\tctcf\ttype' hg19_ubi-rDHS_2_GM12878_five_dimension_matrix2.txt
#
Rscript test_tSNE_on_GM12878.R
