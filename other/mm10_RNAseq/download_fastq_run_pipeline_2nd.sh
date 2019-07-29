#!/bin/bash

# -- Kaili
# This script is for downloading fastq file and run RNAseq pipeline.

num=$1

#######
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"
cd ${workDir}

#
expID=`awk -v num="$num" '{if(NR==num){print $1}}' mouse_RNA_fastq_list_merged.txt`
fileIDs=`awk -v num="$num" '{if(NR==num){print $2}}' mouse_RNA_fastq_list_merged.txt`

# make folder, set path
if [ ! -d ${workDir}rna_exp/${expID} ]; then mkdir ${workDir}rna_exp/${expID}; fi
cd ${workDir}rna_exp/${expID}

# get fastq: download and/or merge
#l=`echo ${fileIDs} | awk '{split($1,a,"-");print length(a)}'`
#if [[ "$l" != 1 ]]; then
#    for ((i=1; i<=$l; i++))
#    do
#        tmp_file_id=`echo ${fileIDs} | awk -v i="$i" '{split($1,a,"-");print a[i]}'`
#        wget https://www.encodeproject.org/files/${tmp_file_id}/@@download/${tmp_file_id}.fastq.gz
#        gzip -d ${tmp_file_id}.fastq.gz
#        cat ${tmp_file_id}.fastq >> ${fileIDs}.fastq
#    done
#else
#    wget https://www.encodeproject.org/files/${fileIDs}/@@download/${fileIDs}.fastq.gz
#    gzip -d ${fileIDs}.fastq.gz
#fi

# run pipeline
STARgenomeDir="/home/fankaili/STARgenome"
RSEMrefDir="/home/fankaili/RSEMgenome/RSEMref_mm10"
dataType="unstr_SE"
#
#mkdir ${workDir}rna_exp/${expID}/${fileIDs}/
echo bash /data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/run_RNAseq_pipeline_2nd.sh ${workDir}rna_exp/${expID}/${fileIDs}.fastq "" ${STARgenomeDir} ${RSEMrefDir} ${dataType} 8 8 ${workDir}rna_exp/${expID}/${fileIDs}/ ${fileIDs}
bash /data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/run_RNAseq_pipeline_2nd.sh ${workDir}rna_exp/${expID}/${fileIDs}.fastq "" ${STARgenomeDir} ${RSEMrefDir} ${dataType} 8 8 ${workDir}rna_exp/${expID}/${fileIDs}/ ${fileIDs}

# bash /data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/run_RNAseq_pipeline.sh /data/zusers/fankaili/ccre/mm10_rnaseq/ENCFF581SPK.fastq "" /home/fankaili/STARgenome /home/fankaili/RSEMgenome/RSEMref_mm10 unstr_SE 8 8 /data/zusers/fankaili/ccre/mm10_rnaseq/ss/ ss

# ----
echo "Done. Cheers!"
