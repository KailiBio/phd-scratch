#!/bin/bash

# -- Kaili
# This script is for analyzing batch effect for all 73 RNA-seq data.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/"
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"

cd ${workDir}

# 1. get expression matrix
## 1) raw signal
#
awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_exp_matrix.txt
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    fileIDs=`awk '{print $2}' <<< ${line}`
    rep=`awk '{print $3}' <<< ${line}`
    sample=`awk '{split($4,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
    #
    echo ${sample}
    #
    head -1 ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results > tmp.gene.txt
    grep "ENSMUSG" ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results >> tmp.gene.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $0,a[$1]}}' tmp.gene.txt mm10_exp_matrix.txt > tmp.txt
    sed "s/TPM/${sample}_${rep}/" tmp.txt > mm10_exp_matrix.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt
rm tmp.*txt



# awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_exp_matrix_rep1.txt
# #
# while read line
# do
#     expID=`awk '{print $1}' <<< ${line}`
#     fileIDs=`awk '{print $2}' <<< ${line}`
#     rep=`awk '{print $3}' <<< ${line}`
#     sample=`awk '{split($4,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
#     #
#     if [ "$rep" -eq "1" ];then
#         echo ${sample}
#     fi
#     #
#     awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $0,a[$1]}}' ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results mm10_exp_matrix_rep1.txt > tmp.txt
#     sed "s/TPM/${sample}_${rep}/" tmp.txt > mm10_exp_matrix_rep1.txt
# done < mouse_RNA_fastq_list_merged_SampleInfo.txt

## 2) log-transformed, z-score signal
awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_exp_zscore_matrix.txt
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    fileIDs=`awk '{print $2}' <<< ${line}`
    rep=`awk '{print $3}' <<< ${line}`
    sample=`awk '{split($4,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
    #
    echo ${sample}
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6;b[$1]=1}else{if(FNR>1 && b[$1]){print $1,a[$1]}}}' ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results mm10_exp_zscore_matrix.txt > tmp.txt
    python /data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/zscore-normalization.py tmp.txt 2 > tmp2.txt
    awk -v sample="${sample}" -v rep="${rep}" 'BEGIN{FS=OFS="\t";print "gene_id",sample"_"rep}{print $1,$2}' tmp2.txt > tmp3.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$1]}}' tmp3.txt mm10_exp_zscore_matrix.txt > tmp4.txt
    mv tmp4.txt mm10_exp_zscore_matrix.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt


# 2. detect batch effect using raw TPM
Rscript ${scriptDir}analyze_batch_effect.R

# 3. get spike-in matrix
echo "gene_id" > mm10_spikeIn_exp_matrix.txt
grep "ERCC-" gencode.vM18.basic.annotation_withSpinkIn.gtf | cut -f 1 | sort -u >> mm10_spikeIn_exp_matrix.txt
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    fileIDs=`awk '{print $2}' <<< ${line}`
    rep=`awk '{print $3}' <<< ${line}`
    sample=`awk '{split($4,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
    #
    echo ${sample}
    #
    head -1 ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results > tmp.spikein.txt
    grep "ERCC-" ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results >> tmp.spikein.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $0,a[$1]}}' tmp.spikein.txt mm10_spikeIn_exp_matrix.txt > tmp.txt
    sed "s/expected_count/${sample}_${rep}/" tmp.txt > mm10_spikeIn_exp_matrix.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt
#
rm tmp*

###
awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_exp_counts_matrix.txt
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    fileIDs=`awk '{print $2}' <<< ${line}`
    rep=`awk '{print $3}' <<< ${line}`
    sample=`awk '{split($4,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
    #
    echo ${sample}
    #
    head -1 ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results > tmp.gene.txt
    grep "ENSMUSG" ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results >> tmp.gene.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $0,a[$1]}}' tmp.gene.txt mm10_exp_counts_matrix.txt > tmp.txt
    sed "s/expected_count/${sample}_${rep}/" tmp.txt > mm10_exp_counts_matrix.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt
rm tmp.*txt


##########
# Jul 28
awk '{FS=OFS="\t"}{print $4}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_exp_matrix_spikeIn96_Jul28.txt
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    fileIDs=`awk '{print $2}' <<< ${line}`
    rep=`awk '{print $3}' <<< ${line}`
    sample=`awk '{split($4,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
    #
    echo ${sample}
    #
    head -1 ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results > tmp.gene.txt
    grep "ENSMUSG" ./rna_exp/${expID}/${fileIDs}/${fileIDs}.genes.results >> tmp.gene.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $0,a[$1]}}' tmp.gene.txt mm10_exp_matrix_spikeIn96_Jul28.txt > tmp.txt
    sed "s/TPM/${sample}_${rep}/" tmp.txt > mm10_exp_matrix_spikeIn96_Jul28.txt
done < mouse_RNA_fastq_list_merged_SampleInfo.txt
rm tmp.*txt

#################
# Aug 08
bash ${scriptDir}clustering_using_vM4_sanityCheck.sh
####
# sanity check if the matrix are correct
## 1) get gene list for M4 and M18
awk '{FS=OFS="\t"}{if(NR>1){split($4,a,".");print $4,a[1]}}' /home/fankaili/genome/mm10_vM4_comprehensive_gene.txt > mm10_M4_genelist.txt
awk '{FS=OFS="\t"}{if(NR>1){split($4,a,".");print $4,a[1]}}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt > mm10_M18_genelist.txt
cut -f 2 mm10_M4_genelist.txt > tmp.txt
cut -f 2 mm10_M18_genelist.txt >> tmp.txt
sort tmp.txt | uniq -d > mm10_M4_M18_overlap_genelist.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$1}else{print $1,a[$1]}}' mm10_M4_genelist.txt mm10_M4_M18_overlap_genelist.txt > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$1}else{print $1,$2,a[$1]}}' mm10_M18_genelist.txt tmp.txt  > mm10_M4_M18_overlap_genelist2.txt
sed -i '1i ID\tM4\tM18' mm10_M4_M18_overlap_genelist2.txt
rm tmp.txt
## 2) make scatter plot for read counts and TPM.
nohup bash ${scriptDir}compare_TPM_differentVersion.sh > ./nohup/nohup.compare_TPM_differetVersion.out 2>&1&
# z001 14419
nohup bash ${scriptDir}compare_readCounts_differetVersion.sh > ./nohup/nohup.compare_readCounts_differetVersion 2>&1&
# z001 14436
##########
# get final matrix
bash ${scriptDir}get_matrix_from_RSEM.sh
