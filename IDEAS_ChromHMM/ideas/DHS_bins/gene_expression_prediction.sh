#!/bin/bash

# -- Kaili
# This script is for validating the gene expression prediction for IDEAS result.
### "Accurate and reproducible functional maps in 127 human cell types via 2D genome segmentation."
### only do ±2kb here, no B splines.
# 1. get gene expression result, make matrix
# 2. get states bed file
# 3. get gene 20 bins bed file
# 4. get state proportion



scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/"

cd ${workDir}

# 1. get gene expression result, make matrix
python ${scriptDir}get_ENCODE_mouse_RNA_data_list.py
# get /data/zusers/fankaili/ideas/mouse_RNA_data_list.txt
# then get 66 samples: /data/zusers/fankaili/ideas/mouse_RNA_66_data_list.txt

# make data list
if [ -f mm10_RNA_tpm_matrix.txt ];then
    rm mm10_RNA_tpm_matrix.txt
fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    expID=`awk '{print $2}' <<< ${line}`
    fileID=`awk '{print $3}' <<< ${line}`
    #
    echo -e "geneID\t"${sample} > tmp.txt
    awk '{FS=" ";OFS="\t"}{NR==FNR="\t"}{if($1~/^ENS/){print $1,$6}}' /data/projects/encode/data/${expID}/${fileID}.tsv >> tmp.txt
    if [ -f mm10_RNA_tpm_matrix.txt ]
    then
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$1]}}' tmp.txt mm10_RNA_tpm_matrix.txt \
        > tmp2.txt
        mv tmp2.txt mm10_RNA_tpm_matrix.txt
    else
        mv tmp.txt mm10_RNA_tpm_matrix.txt
    fi
done < /data/zusers/fankaili/ideas/mouse_RNA_66_data_list.txt

### get mm10 protein coding gene
cd /home/fankaili/genome/
wget https://www.encodeproject.org/files/gencode.vM4.annotation/@@download/gencode.vM4.annotation.gtf.gz
gzip -d gencode.vM4.annotation.gtf.gz
#
awk '{FS=OFS="\t"}{if(match($9,/protein_coding/) && $3=="gene"){split($9,a,"\"");print a[2]}}' gencode.vM4.annotation.gtf \
| sort -u > mm10_vM4_protein_coding_ID.txt
# 22,032 protein coding gene: /home/fankaili/genome/mm10_vM4_protein_coding_ID.txt
awk '{FS=OFS="\t"}{if(match($9,/protein_coding/) && $3=="gene"){split($9,a,"\"");print $1,$4,$5,a[2],$6,$7}}' gencode.vM4.annotation.gtf \
| sort -u | sort -k1,1 -k2,2n > mm10_vM4_protein_coding.bed


### get protein_coding gene expression matrix
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $0}}}' \
/home/fankaili/genome/mm10_vM4_protein_coding.bed mm10_RNA_tpm_matrix.txt > \
mm10_RNA_protein-coding_tpm_matrix.txt

# 2. get states bed file
## dhs_bins
stateDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_v3_100-400bp_result/"
prefix="DHS_v3_100-400bp."
stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/"
nohup bash ${scriptDir}make_each_sample_state_bed_IDEAS.sh ${stateDir} ${prefix} ${stateBedDir} 66 \
> /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/nohup.make_each_sample_state_bed_dhs.out 2>&1&

## normal_bins
stateDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/66samples_10marks_normal_bins_result/"
prefix="66samples_10marks_normal_bins."
stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/state_bed/"
nohup bash ${scriptDir}make_each_sample_state_bed_IDEAS.sh ${stateDir} ${prefix} ${stateBedDir} 66 \
> /data/zusers/fankaili/ideas/dhs_bins/normal_bins/nohup.make_each_sample_state_bed_normal.out 2>&1&


# 3. get gene 20 bins bed file
awk '{FS=OFS="\t"}{if($6=="+"){s=$2-2000;for(i=1;i<=20;i++){print $1,s+200*(i-1),s+200*i,$4"_"i}}else{s=$3+2000;for(i=1;i<=20;i++){print $1,s-200*i,s-200*(i-1),$4"_"i}}}' \
/home/fankaili/genome/mm10_vM4_protein_coding.bed > mm10_protein_coding_promoter_bins.bed

# 4. get state proportion
count_state_proportion(){
    stateBedDir=$1
    prefix=$2
    #
    while read line
    do
        echo $line
        ### get state in each window
        intersectBed -a mm10_protein_coding_promoter_bins.bed -b \
        ${stateBedDir}${line}_state_sorted.bed -wa -wb | sort -u | \
        awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2],$5,$6,$7,$8,$9}' \
        | sort -k6,6n -k5,5 -k8,8n> ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state.txt
        ### state transfer
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$1}else{split(a[$11],b,"_");print $0,b[2]}}' \
        ${prefix}_bins_transfer.txt ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state.txt > tmp.${line}_${prefix}.txt
        ### calculate state count in 20 windows
        for i in {1..20}
        do
            awk -v n="$i" 'BEGIN{FS=OFS="\t";for(i=1;i<=20;i++){a[i]=0}}{if($6==n){a[$12]+=1}}END{for(i=1;i<=20;i++){print i,a[i]}}' \
            tmp.${line}_${prefix}.txt > mm10_gene_${line}_${prefix}_state_count_window_${n}.txt
        done
        rm tmp.${line}_${prefix}.txt
    done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt
}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/"
prefix="dhs"
count_state_proportion ${stateBedDir} ${prefix}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/state_bed/"
prefix="normal"
count_state_proportion ${stateBedDir} ${prefix}

# nohup bash ss1.sh > nohup.ss1.out 2>&1&
# nohup bash ss2.sh > nohup.ss2.out 2>&1&
