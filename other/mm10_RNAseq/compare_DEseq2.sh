#!/bin/bash

# -- Kaili
# This script is for comparing DEseq2 output between TPM & TPM-log-limma-unlog.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/"
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"

cd ${workDir}

# 0.
cut -f 5,6 mouse_RNA_fastq_list_merged_SampleInfo.txt | sort -u > mouse_sample_list.txt
#
compare_DEcall(){
    outfile=$1
    DEG_folder=$2
    #
    if [ -f ${outfile} ];then rm ${outfile}; fi
    #
    for file in `ls ./${DEG_folder}/`
    do
        sample1=`awk '{split($0,a,"_VS_");print a[1]}' <<< ${file}`
        sample2=`awk '{split($1,a,"_VS_");split(a[2],b,".txt");print b[1]}' <<< ${file}`
        echo ${sample1}" VS "${sample2}
        #
        if [[ $sample1 == "facial"* ]]; then
            sample11=${sample1/facial/embryonic_facial_prominence}
        elif [[ $sample1 == "neural.tube"* ]];then
            sample11=${sample1/neural.tube/neural_tube}
        else
            sample11=$sample1
        fi
        #
        if [[ $sample2 == *"facial"* ]]; then
            sample22=${sample2/facial/embryonic_facial_prominence}
        elif [[ $sample2 == "neural.tube"* ]];then
            sample22=${sample2/neural.tube/neural_tube}
        else
            sample22=$sample2
        fi
        #
        if [ -f /data/projects/screen/Version-4/mouse_epigenome/de_all_pairs/data/${sample11}_VS_${sample22}.txt.gz ];then
            cp /data/projects/screen/Version-4/mouse_epigenome/de_all_pairs/data/${sample11}_VS_${sample22}.txt.gz ./
            gzip -d ${sample11}_VS_${sample22}.txt.gz
            awk '{if(NR>1 && $7<0.05 && ($3>2 || $3<-2)){print $1}}' ${sample11}_VS_${sample22}.txt > tmp.1.txt
            #
            rm ${sample11}_VS_${sample22}.txt
        elif [ -f /data/projects/screen/Version-4/mouse_epigenome/de_all_pairs/data/${sample22}_VS_${sample11}.txt.gz ];then
            cp /data/projects/screen/Version-4/mouse_epigenome/de_all_pairs/data/${sample22}_VS_${sample11}.txt.gz ./
            gzip -d ${sample22}_VS_${sample11}.txt.gz
            awk '{if(NR>1 && $7<0.05 && ($3>2 || $3<-2)){print $1}}' ${sample22}_VS_${sample11}.txt > tmp.1.txt
            #
            rm ${sample22}_VS_${sample11}.txt
        fi
        #
        awk '{if(NR>1 && $7<0.05 && ($3>2 || $3<-2)){print $1}}' ./${DEG_folder}/${file} > tmp.3.txt
        #
        a=`wc -l tmp.1.txt | awk '{print $1}'`
        b=`cat tmp.1.txt tmp.3.txt | sort | uniq -d | wc -l | awk '{print $1}'`
        c=`wc -l tmp.3.txt | awk '{print $1}'`
        echo -e ${sample1}"_VS_"${sample2}"\t"$a"\t"$b"\t"$c >> ${outfile}
    done
}

# 1. get DE count matrix, Junko's vs. TPM-log-limma-unlog
outfile="DE_call_comparison_M4_TPM-log-limma-unlog.txt"
DEG_folder="de_M4"
#
compare_DEcall ${outfile} ${DEG_folder}
Rscript ${scriptDir}make_DEcall_barplot.R "/data/zusers/fankaili/ccre/mm10_rnaseq/" "DE_call_comparison_M4_TPM-log-limma-unlog.txt" "DE_comparison_TPM-log-limma-unlog.pdf"

# 2.  get DE count matrix, Junko's vs. read_counts
outfile="DE_call_comparison_M4_read-counts.txt"
DEG_folder="de_M4_read"
#
compare_DEcall ${outfile} ${DEG_folder}
Rscript ${scriptDir}make_DEcall_barplot.R "/data/zusers/fankaili/ccre/mm10_rnaseq/" "DE_call_comparison_M4_read-counts.txt" "DE_comparison_read-counts.pdf"



make_sanityCheck_volcano.sh
