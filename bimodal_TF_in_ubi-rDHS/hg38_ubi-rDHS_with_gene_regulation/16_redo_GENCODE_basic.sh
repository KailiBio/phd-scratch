#!/bin/bash

# -- Kaili
# This script is for redo all the TSS analysis using V28, GENCODE basic annotation.
# 0. get GENCODE v28 basic annotation
# 00. process gene & TSS file
# 1. fig1: signal
# 2. fig2: TSindex
# 3. fig3: distance

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
autoScriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/automated_pipeline/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/"

hg38_ubi_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed"
hg38_rOCRs="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed"

####################################
cd ${workDir}

##----------------------------
# 0. get GENCODE v28 basic annotation
cd /home/fankaili/genome/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.basic.annotation.gtf.gz
gzip -d gencode.v28.basic.annotation.gtf.gz
# # get filter gene type
# bash ${dailyCodeDir}process_GENCODE_gtf.sh /home/fankaili/genome/gencode.v24.annotation.gtf hg38 \
# v24 comprehensive /home/fankaili/genome/
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0}}}' \
# /data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered.bed hg38_v24_comprehensive_transcript.txt | cut -f 10 | sort -u \
# > transcript_type_for_filter_TSS.txt
# process gtf
bash ${dailyCodeDir}process_GENCODE_gtf.sh /home/fankaili/genome/gencode.v28.basic.annotation.gtf hg38 v28 basic /home/fankaili/genome/

##----------------------------
# 00. process gene&TSS file
bash ${autoScriptDir}0_preprocess_gene_TSS_files.sh hg38_v28_basic_TSS_filtered.bed \
hg38_v28_basic_TSS_protein_coding.bed hg38_v28_basic_gene_protein_coding.txt \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

## parameters
DNase_signal_path="/data/zusers/moorej3/moorej.ghpcc.project/PsychENCODE/Registry/V0-hg38/signal-output/"
non_ubi_outPath="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/"
gene_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_gene_labled.bed"
tss_labeled_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_uniq_labeled.bed"
exp_signal_path="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/"
tss_signal_path="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/tissue_rampage_v28/"
tss_with_id_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/hg38_v28_basic_TSS_filtered_with_uniqID.bed"


##----------------------------
# 1. fig1: signal
## get DNase-seq, RNA-seq, RAMPAGE signal for
bash ${autoScriptDir}1_compare_signal.sh ${DNase_signal_path} ${non_ubi_outPath} ${gene_labeled_file} ${tss_labeled_file} ${exp_signal_path} ${tss_signal_path} ${workDir}


##----------------------------
# 2. fig2: TSindex
bash ${autoScriptDir}2_compare_TSindex.sh ${gene_labeled_file} ${tss_with_id_file} ${tss_labeled_file} ${workDir}


##----------------------------
# 3. fig3: distance
bash ${autoScriptDir}3_TSS_cluster.sh ${tss_labeled_file} ${workDir}

# ## 1) get uniq TSS genes
# cut -f 7 hg38_v28_basic_TSS_filtered_uniq_labeled.bed | sort | uniq -u > \
# hg38_v28_basic_TSS_filtered_uniq_labeled_uniqTSS.txt
# # 37,303 singular TSS genes out of 54,176 genes
# intersectBed -a hg38_v28_basic_TSS_filtered_uniq_labeled.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed \
# -wa -u > hg38_v28_basic_TSS_overlapping_ubi-rOCRs.bed
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$7]=1}else{if(a[$1]){print $0}}}' hg38_v28_basic_TSS_overlapping_ubi-rOCRs.bed \
# hg38_v28_basic_TSS_filtered_uniq_labeled_uniqTSS.txt | wc -l
# ### 1219 genes overlapping ubi-rOCRs
# ## 2) for multiple TSS genes
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]!=1){print $0}}}' \
# hg38_v28_basic_TSS_filtered_uniq_labeled_uniqTSS.txt hg38_v28_basic_TSS_filtered_uniq_labeled.bed | sort -k7,7 -k2,2n> hg38_v28_basic_TSS_filtered_uniq_labeled_multiple.bed
# awk 'BEGIN{FS=OFS="\t";gene="";tss="";end="";distance=1000000}{if(gene==""){gene=$7;tss=$4;end=$3}else{if(gene==$7){if(distance>($2-end)){print gene,tss,$2-end}else{print gene,tss,distance};tss=$4;distance=$2-end;end=$3}else{print gene,tss,distance;gene=$7;tss=$4;end=$3;distance=1000000}}}END{print gene,tss,distance}' \
# hg38_v28_basic_TSS_filtered_uniq_labeled_multiple.bed > tmp.txt
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$0}else{print a[$4],$8,$9}}' tmp.txt \
# hg38_v28_basic_TSS_filtered_uniq_labeled_multiple.bed > hg38_v28_basic_TSS_multiple_distance_sameGene.bed


##----------------------------
# 4. for protein-coding only
bash ${autoScriptDir}4_repeat_for_protein_coding.sh ${non_ubi_outPath} ${gene_labeled_file} ${tss_labeled_file} ${exp_signal_path} ${tss_signal_path} ${tss_with_id_file} ${workDir}"protein_coding/"
