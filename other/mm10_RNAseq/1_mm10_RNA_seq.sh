#!/bin/bash

# -- Kaili
# This script is for reruning mm10 RNA-seq to GENCODE vM18.
# 1. get data list
# 2. download bam files
# 3. run RSEM
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/"
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"

cd ${workDir}

# 0. STAR & RSEM preparation
## z018
cd /home/fankaili/
mkdir STARgenome
mkdir RSEMgenome
#
fastaGenome="/data/zusers/fankaili/ccre/mm10_rnaseq/mm10.fa"
fastaSpikeins="/data/tusers/yutianx/tongji2/piRNA/For_Xuan/rawdata/rnaseq/pi6_Jan29.2019/spikein.fa"
gtf="/data/zusers/fankaili/ccre/mm10_rnaseq/gencode.vM18.basic.annotation_withSpinkIn.gtf"
## 1) STAR
STARgenomeDir="/home/fankaili/STARgenome"
STARcommand="STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $STARgenomeDir --genomeFastaFiles $fastaGenome $fastaSpikeins --sjdbGTFfile $gtf --sjdbOverhang 100 --outFileNamePrefix $STARgenomeDir"
echo $STARcommand
$STARcommand
## 2) RSEM
RSEMgenomeDir="/home/fankaili/RSEMgenome"
RSEMcommand="rsem-prepare-reference --gtf $gtf $fastaGenome","$fastaSpikeins $RSEMgenomeDir/RSEMref_mm10"
echo $RSEMcommand
$RSEMcommand
# https://www.encodeproject.org/files/ENCFF931IVO/@@download/ENCFF931IVO.fastq.gz

cd ${workDir}

# 1. get data list
python ${scriptDir}get_mm10_RNAseq_fastq.py
#
sort -u mouse_RNA_fastq_list.txt | sort -k1,1 -k3,3 -k4,4 > mouse_RNA_fastq_list_sort.txt
awk 'BEGIN{FS=OFS="\t";e="a";file="a";rep="a"}{if(NR==1){e=$1;file=$2;rep=$3}else if(e==$1 && rep==$3){a=file"-"$2;file=a}else{print e,file,rep;e=$1;file=$2;rep=$3}}END{print e,file,rep}' mouse_RNA_fastq_list_sort.txt > mouse_RNA_fastq_list_merged.txt

# 2. download fastq files, STAR+RSEM
mkdir rna_exp
#
for num in {1..146}
do
    echo $num
    bash ${scriptDir}download_fastq_run_pipeline.sh $num
done

for num in {1..146}
do
    echo $num
    bash ${scriptDir}run_RSEM.sh $num
done
