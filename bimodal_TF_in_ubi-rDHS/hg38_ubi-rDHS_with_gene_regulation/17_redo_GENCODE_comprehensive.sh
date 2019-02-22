#!/bin/bash

# -- Kaili
# This script is for redo all the TSS analysis using V28, GENCODE comprehensive annotation (removing problematic transcripts).
# 0. get GENCODE v28 comprehensive annotation
# 1. process TSS file
# 2. fig1: signal
# 3. fig2: TSindex
# 4. fig3: distance

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
autoScriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/automated_pipeline/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/comprehensive_annotation/"

hg38_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed"
hg38_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed"

####################################
cd ${workDir}

##----------------------------
# 0. get GENCODE v28 comprehensive annotation
cd /home/fankaili/genome/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz
gzip -d gencode.v28.annotation.gtf.gz
bash ${dailyCodeDir}process_GENCODE_gtf.sh /home/fankaili/genome/gencode.v28.annotation.gtf hg38 \
v28 comprehensive /home/fankaili/genome/

awk '{FS=OFS="\t"}{print $1,$2-50,$3+50,$4,1,$6}' /home/fankaili/genome/hg38_v28_comprehensive_TSS.bed | sort -k1,1 -k2,2n > /home/fankaili/genome/hg38_v28_comprehensive_TSS_sorted.bed
# nohup python ${autoScriptDir}get_GRCh38_tissue_RAMPAGE_exp.py /home/fankaili/genome/hg38_v28_comprehensive_TSS_sorted.bed /data/zusers/fankaili/ccre/hg38_ubi-rDHS/tissue_rampage_v28/ > nohup.out 2>&1&

##----------------------------
# 00. process gene&TSS file
bash ${autoScriptDir}0_preprocess_gene_TSS_files.sh hg38_v28_comprehensive_TSS_filtered.bed hg38_v28_comprehensive_TSS_protein_coding.bed hg38_v28_comprehensive_gene_protein_coding.txt /data/zusers/fankaili/ccre/hg38_ubi-rDHS/comprehensive_annotation/

## parameters
DNase_signal_path="/data/zusers/moorej3/moorej.ghpcc.project/PsychENCODE/Registry/V0-hg38/signal-output/"
non_ubi_outPath="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/"
gene_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/comprehensive_annotation/hg38_v28_comprehensive_gene_labled.bed"
tss_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/comprehensive_annotation/hg38_v28_comprehensive_TSS_filtered_uniq_labeled.bed"
exp_signal_path="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/"
tss_signal_path="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/tissue_rampage_v28/"
tss_with_id_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/comprehensive_annotation/hg38_v28_comprehensive_TSS_filtered_with_uniqID.bed"

##----------------------------
# 1. fig1: signal
## get DNase-seq, RNA-seq, RAMPAGE signal for
bash ${autoScriptDir}1_compare_signal.sh ${DNase_signal_path} ${non_ubi_outPath} ${gene_labeled_file} ${tss_labeled_file} ${tss_with_id_file} ${exp_signal_path} ${tss_signal_path} ${workDir}


##----------------------------
# 2. fig2: TSindex
bash ${autoScriptDir}2_compare_TSindex.sh ${gene_labeled_file} ${tss_with_id_file} ${tss_labeled_file} ${workDir}


##----------------------------
# 3. fig3: distance
bash ${autoScriptDir}3_TSS_cluster.sh ${tss_labeled_file} ${workDir}


##----------------------------
# 4. for protein-coding only
mkdir ${workDir}protein_coding
bash ${autoScriptDir}4_repeat_for_protein_coding.sh ${non_ubi_outPath} ${gene_labeled_file} ${tss_labeled_file} ${exp_signal_path} ${tss_signal_path} ${tss_with_id_file} ${workDir}"protein_coding/"
