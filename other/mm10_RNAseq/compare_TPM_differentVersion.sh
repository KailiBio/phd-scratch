#!/bin/bash

# -- Kaili
# This script is for comparing RNA-seq output in different version. (M4 vs. M18)


scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/"
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"

cd ${workDir}

while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    fileIDs=`awk '{print $2}' <<< ${line}`
    rep=`awk '{print $3}' <<< ${line}`
    sample=`awk '{split($4,a,"-");if(length(a)==5){split(a[4],b,"(");print a[2]"_"b[2]}else if(length(a)==6){split(a[5],b,"(");print a[2]"-"a[3]"_"b[2]}else if(length(a)==7){split(a[6],b,"(");print a[3]"_"b[2]}}' <<< ${line}`
    echo $sample
    sample_rep=$sample"_"$rep
    # M4
    tsvID=`grep $expID mouse_RNA_tsv_list.txt | awk -v rep="$rep" '{if($3==rep){print $2}}'`
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $0,a[$2]}}' /data/projects/encode/data/${expID}/${tsvID}.tsv mm10_M4_M18_overlap_genelist2.txt > tmp1.txt
    # spikeIn96, M18
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $0,a[$3]}}' ./rna_exp_Jul25/${expID}/${fileIDs}/${fileIDs}.genes.results tmp1.txt > tmp2.txt
    # make scatter
    cut -f 1,4,5 tmp2.txt > tmp.txt
    Rscript ${scriptDir}make_rnaseq_check_scatter.R tmp.txt ${sample_rep}
done < mouse_RNA_fastq_list_merged_SampleInfo.txt

rm tmp*.txt
