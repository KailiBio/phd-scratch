#!/bin/bash

# -- Kaili
# This script is for calling differential expressed gene.
# 1. get M4 tsv file list
# 2. get M4 TPM & read_counts matrix
# 3. Call DE using M4 data
# 4. comparison with Junko's results

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/"
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"

cd ${workDir}

# 1. get M4 tsv file list
python ${scriptDir}get_mm10_RNAseq_tsv_list.py
sort -u mouse_RNA_tsv_list.txt | sort -k5,5 -k3,3n > tmp.txt
mv tmp.txt mouse_RNA_tsv_list.txt

# 2. get M4 TPM & read_counts matrix
awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM4_comprehensive_gene.txt > mm10_M4_TPM_matrix.txt
awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM4_comprehensive_gene.txt > mm10_M4_counts_matrix.txt
#
while read line
do
    expID=`awk '{print $1}' <<< $line`
    fileID=`awk '{print $2}' <<< $line`
    rep=`awk '{print $3}' <<< ${line}`
    sample=`awk '{split($5,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
    echo $sample
    #
    head -1 /data/projects/encode/data/${expID}/${fileID}.tsv > tmp.gene.txt
    grep "ENSMUSG" /data/projects/encode/data/${expID}/${fileID}.tsv >> tmp.gene.txt
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $0,a[$1]}}' tmp.gene.txt mm10_M4_TPM_matrix.txt > tmp.txt
    sed "s/TPM/${sample}_${rep}/" tmp.txt > mm10_M4_TPM_matrix.txt
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $0,a[$1]}}' tmp.gene.txt mm10_M4_counts_matrix.txt > tmp.txt
    sed "s/expected_count/${sample}_${rep}/" tmp.txt > mm10_M4_counts_matrix.txt
done < mouse_RNA_tsv_list.txt

# 3. transfer TPM-log-limma-unlog back to read-counts
## get sum_RPK for each sample
if [ -f mm10_RPK_eachSample.txt ];then rm mm10_RPK_eachSample.txt; fi
#
while read line
do
    expID=`awk '{print $1}' <<< $line`
    fileID=`awk '{print $2}' <<< $line`
    rep=`awk '{print $3}' <<< ${line}`
    sample=`awk '{split($5,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
    echo $sample
    #
    grep "ENSMUSG" /data/projects/encode/data/${expID}/${fileID}.tsv > tmp.1.txt
    n=`awk '{FS=OFS="\t"}{if($6!=0){a=$5*1000/$6/$4;sum+=a;n++}}END{print sum/n}'  tmp.1.txt`
    echo -e $sample"_"$rep"\t"$n >> mm10_RPK_eachSample.txt
done < mouse_RNA_tsv_list.txt
## get gene length
grep "ENSMUSG" /data/projects/encode/data/ENCSR178GUS/ENCFF434XNZ.tsv | cut -f 1,4 > mm10_M4_gene_length.txt



# 3. Call DE using M4 data
## 1)
# Rscript normalized_M4.R
nohup Rscript ${scriptDir}call_DE.R "mm10_M4_after-limma.txt" "M4" > ./nohup/nohup.call_DE_M4.out 2>&1&
#
nohup Rscript ${scriptDir}call_DE.R "mm10_M4_read_66.txt" "M4_read" 11 > ./nohup/nohup.call_DE_M4_read.out 2>&1&
# z001 25234

# 4. comparison with Junko's results
bash ${scriptDir}compare_DEseq2.sh

# 5. call DE using betaPrior=TRUE
nohup Rscript ${scriptDir}call_DE.R "mm10_M4_read_66.txt" "M4_read_TRUE_66" 11 > ./nohup/nohup.call_DE_M4_read_TRUE.out 2>&1&
# z001 48121
nohup bash ${scriptDir}sanity_check_vM4_TRUE.sh > ./nohup/nohup.sanity_check_vM4_TRUE.out 2>&1&
# z001 655
#######
nohup Rscript ${scriptDir}call_DE.R "mm10_M4_counts_matrix_sorted.txt" "M4_read_TRUE" > ./nohup/nohup.call_DE_M4_read_TRUE_all.out 2>&1&
# z001 42510


# 6. call DE using M18 data (betaPrior=TRUE)
nohup Rscript ${scriptDir}call_DE.R "mm10_M18_readCounts_matrix.txt" "M18_readCounts" > ./nohup/nohup.call_DE_M18_readCounts.out 2>&1&
# z001 44502

# for file in `ls /data/zusers/fankaili/ccre/mm10_rnaseq/de_M18_readCounts/`
# do
#     if [ ! -f /data/zusers/fankaili/ccre/mm10_rnaseq/de_M4_read_TRUE/${file} ];then
#         echo $file
#     fi
# done

### compare DEGs: M18 vs. M4
nohup bash ${scriptDir}compare_DEGs_M4_M18.sh > ./nohup/nohup.compare_DEGs_M4_M18.out 2>&1&
# z001 18269
