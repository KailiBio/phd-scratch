#!/bin/bash

# -- Kaili
# This script is for validating the gene expression prediction for IDEAS result (within cell types).
### "Accurate and reproducible functional maps in 127 human cell types via 2D genome segmentation."
### only do ±2kb here, no B splines.
# 1. get gene expression result, make matrix
# 2. get states bed file
# 3. get gene 20 bins bed file
# 4. get state proportion
# 5. do regression

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
head -1 mm10_RNA_tpm_matrix.txt > mm10_RNA_protein-coding_tpm_matrix.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $0}}}' \
/home/fankaili/genome/mm10_vM4_protein_coding.bed mm10_RNA_tpm_matrix.txt | sort -k1,1 >> \
mm10_RNA_protein-coding_tpm_matrix.txt

######---------------------
## try promoter with H3K4me3
# get TSS bins
## sort all bins
sort -k1,1 -k2,2n /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_signal_based.bed | awk '{FS=OFS="\t"}{print $0,NR}' > ocr_bins.bed
sort -k1,1 -k2,2n /data/zusers/fankaili/ideas/dhs_bins/normal_bins/mm10_tab_noM.bed | awk '{FS=OFS="\t"}{print $0,NR}' > normal_bins.bed
## get TSS bins
awk '{FS=OFS="\t"}{if($6=="+"){print $1,$2-1,$2,$4,$5,$6}else{print $1,$3,$3+1,$4,$5,$6}}' /home/fankaili/genome/mm10_vM4_protein_coding.bed > mm10_protein_coding_TSS.bed
intersectBed -a ocr_bins.bed -b mm10_protein_coding_TSS.bed -wa -wb | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$9}' | sort -u | sort -k1,1 -k2,2n > mm10_protein_coding_TSS_ocr_bin.txt
intersectBed -a normal_bins.bed -b mm10_protein_coding_TSS.bed -wa -wb | awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$9}' | sort -u | sort -k1,1 -k2,2n > mm10_protein_coding_TSS_normal_bin.txt
## get TSS promoter bins
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$0}else{for(i=1;i<=10;i++){x=10+i;print a[$5+i],$6"_"x;y=10-i;print a[$5-i],$6"_"(y+1)}}}' ocr_bins.bed mm10_protein_coding_TSS_ocr_bin.txt | sort -k5,5n > mm10_protein_coding_TSS_ocr_promoter_bins.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$0}else{for(i=1;i<=10;i++){x=10+i;print a[$5+i],$6"_"x;y=10-i;print a[$5-i],$6"_"(y+1)}}}' normal_bins.bed mm10_protein_coding_TSS_ocr_bin.txt | sort -k5,5n > mm10_protein_coding_TSS_normal_promoter_bins.txt

# try regression using H3K4me3
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$4]=$6}else{if(a[$1]){print b[$1],$5}}}' mm10_protein_coding_TSS_ocr_bin.txt /data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/heart_12.5_H3K4me3_dhs.tab | sort -k1,1 > try_regression_H3K4me3_heart12.5_ocr.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$4]=$6}else{if(a[$1]){print b[$1],$5}}}' mm10_protein_coding_TSS_normal_bin.txt /data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/heart_12.5_H3K4me3_normal.tab | sort -k1,1 > try_regression_H3K4me3_heart12.5_normal.txt
##
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$6}}}' try_regression_H3K4me3_heart12.5_ocr.txt /data/projects/encode/data/ENCSR150CUE/ENCFF345NEQ.tsv \
| sort -k1,1 > gene_exp_tmp_ocr_heart12.5.txt

######---------------------
# 2. get states bed file
## dhs_bins
stateDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_v3_100-400bp_result/"
prefix="DHS_v3_100-400bp."
stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/"
nohup bash ${scriptDir}make_each_sample_state_bed_IDEAS.sh ${stateDir} ${prefix} ${stateBedDir} 66 > /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/nohup.make_each_sample_state_bed_dhs.out 2>&1&

## normal_bins
stateDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/66samples_10marks_normal_bins_result/"
prefix="66samples_10marks_normal_bins."
stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/state_bed/"
nohup bash ${scriptDir}make_each_sample_state_bed_IDEAS.sh ${stateDir} ${prefix} ${stateBedDir} 66 > /data/zusers/fankaili/ideas/dhs_bins/normal_bins/nohup.make_each_sample_state_bed_normal.out 2>&1&


# 3. get gene 20 bins bed file
awk '{FS=OFS="\t"}{if($6=="+"){s=$2-2000;for(i=1;i<=20;i++){print $1,s+200*(i-1),s+200*i,$4"_"i}}else{s=$3+2000;for(i=1;i<=20;i++){print $1,s-200*i,s-200*(i-1),$4"_"i}}}' \
/home/fankaili/genome/mm10_vM4_protein_coding.bed > mm10_protein_coding_promoter_bins.bed

# 4. get state proportion
count_state_proportion(){
    stateBedDir=$1
    prefix=$2
    state_num=$3
    #
    while read line
    do
        echo $line
        # ### get state in each window
        # intersectBed -a mm10_protein_coding_promoter_bins.bed -b \
        # ${stateBedDir}${line}_state_sorted.bed -wo | sort -u | \
        # awk '{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,$4,a[1],a[2],$5,$6,$7,$8,$10,$9}' \
        # | sort -k6,6n -k5,5 -k8,8n > ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state.txt
        # ### state transfer
        # awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=$1}else{split(a[$11],b,"_");print $0,b[2]}}' \
        # ${prefix}_bins_transfer.txt ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state.txt \
        # > tmp.${line}_${prefix}.txt
        ### calculate state count in 20 windows
        for i in {1..20}
        do
            awk -v n="$i" -v state_num="$state_num" 'BEGIN{FS=OFS="\t";gene="";for(i=0;i<state_num;i++){a[i]=0}}{l=$3-$2;if($6==n){if(gene==""){gene=$5;a[$12]+=$11/l}else if(gene==$5){a[$12]+=$11/l}else{printf gene"\t"n;for(i=0;i<state_num;i++){printf "\t"a[i]};printf "\n";gene=$5;for(i=0;i<state_num;i++){a[i]=0};a[$12]+=$11/l}}}END{printf gene"\t"n;for(i=0;i<state_num;i++){printf "\t"a[i]};printf "\n"}' ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state.txt  | sort -k1,1 > ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state_count_window_${i}.txt
        done
    done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt
}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/"
prefix="dhs"
state_num=43
count_state_proportion ${stateBedDir} ${prefix} ${state_num}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/state_bed/"
prefix="normal"
state_num=44
count_state_proportion ${stateBedDir} ${prefix} ${state_num}

# nohup bash ss1.sh > nohup.ss1.out 2>&1&
# nohup bash ss2.sh > nohup.ss2.out 2>&1&

# 5. do regression
## get matched expression
head -1 mm10_RNA_protein-coding_tpm_matrix.txt > mm10_RNA_protein-coding_tpm_matrix_matched.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$0}else{if(a[$1]){print b[$1]}}}' \
mm10_RNA_protein-coding_tpm_matrix.txt ./normal_state_proportion/mm10_gene_lung_0_normal_state_count_window_9.txt >> mm10_RNA_protein-coding_tpm_matrix_matched.txt

# Rscript do_within_gene_regression.R

#####################
# try ±200bp, ±400bp, ±500bp, ±2kb
# 1) get bins
awk '{FS=OFS="\t"}{if($6=="+"){print $1,$2-200, $2,$4}else{print $1,$3, $3+200,$4}}' /home/fankaili/genome/mm10_vM4_protein_coding.bed > mm10_protein_coding_TSSup200_bins.bed
awk '{FS=OFS="\t"}{if($6=="+"){print $1,$2-200, $2+200,$4}else{print $1,$3-200, $3+200,$4}}' /home/fankaili/genome/mm10_vM4_protein_coding.bed > mm10_protein_coding_TSS200_bins.bed
awk '{FS=OFS="\t"}{if($6=="+"){print $1,$2-400, $2+400,$4}else{print $1,$3-400, $3+400,$4}}' /home/fankaili/genome/mm10_vM4_protein_coding.bed > mm10_protein_coding_TSS400_bins.bed
awk '{FS=OFS="\t"}{if($6=="+"){print $1,$2-500, $2+500,$4}else{print $1,$3-500, $3+500,$4}}' /home/fankaili/genome/mm10_vM4_protein_coding.bed > mm10_protein_coding_TSS500_bins.bed
awk '{FS=OFS="\t"}{if($6=="+"){print $1,$2-2000, $2+2000,$4}else{print $1,$3-2000, $3+2000,$4}}' /home/fankaili/genome/mm10_vM4_protein_coding.bed > mm10_protein_coding_TSS2000_bins.bed
# 2) calculate state proportion
count_state_proportion_bigbins(){
    stateBedDir=$1
    prefix=$2
    state_num=$3
    suffix=$4
    #
    while read line
    do
        echo $line
        ### get state in each window
        intersectBed -a mm10_protein_coding_TSS${suffix}_bins.bed -b \
        ${stateBedDir}${line}_state_sorted.bed -wo | sort -u | \
        awk '{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$10,$9}' \
        | sort -k4,4n > ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state_TSS${suffix}.txt
        ### calculate state count in 20 windows
        # awk -v state_num="$state_num" 'BEGIN{FS=OFS="\t";gene="";for(i=0;i<state_num;i++){a[i]=0};l=200}{if(gene==""){gene=$4;a[$10]+=$9/l}else if(gene==$4){a[$10]+=$9/l}else{printf gene;for(i=0;i<state_num;i++){printf "\t"a[i]};printf "\n";gene=$4;for(i=0;i<state_num;i++){a[i]=0};a[$10]+=$9/l}}END{printf gene;for(i=0;i<state_num;i++){printf "\t"a[i]};printf "\n"}' ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state_TSS${suffix}.txt  | sort -k1,1 > ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state_TSS${suffix}_count.txt
        awk -v state_num="$state_num" -v suffix="$suffix" 'BEGIN{FS=OFS="\t";gene="";for(i=0;i<state_num;i++){a[i]=0};l=suffix*2}{if(gene==""){gene=$4;a[$10]+=$9/l}else if(gene==$4){a[$10]+=$9/l}else{printf gene;for(i=0;i<state_num;i++){printf "\t"a[i]};printf "\n";gene=$4;for(i=0;i<state_num;i++){a[i]=0};a[$10]+=$9/l}}END{printf gene;for(i=0;i<state_num;i++){printf "\t"a[i]};printf "\n"}' ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state_TSS${suffix}.txt  | sort -k1,1 > ./${prefix}_state_proportion/mm10_gene_${line}_${prefix}_state_TSS${suffix}_count.txt
    done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt
}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/"
prefix="dhs"
state_num=43
suffix=up200
count_state_proportion_bigbins ${stateBedDir} ${prefix} ${state_num} ${suffix}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/"
prefix="dhs"
state_num=43
suffix=200
count_state_proportion_bigbins ${stateBedDir} ${prefix} ${state_num} ${suffix}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/"
prefix="dhs"
state_num=43
suffix=400
count_state_proportion_bigbins ${stateBedDir} ${prefix} ${state_num} ${suffix}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/"
prefix="dhs"
state_num=43
suffix=500
count_state_proportion_bigbins ${stateBedDir} ${prefix} ${state_num} ${suffix}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/"
prefix="dhs"
state_num=43
suffix=2000
count_state_proportion_bigbins ${stateBedDir} ${prefix} ${state_num} ${suffix}
###
stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/state_bed/"
prefix="normal"
state_num=44
suffix=up200
count_state_proportion_bigbins ${stateBedDir} ${prefix} ${state_num} ${suffix}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/state_bed/"
prefix="normal"
state_num=44
suffix=200
count_state_proportion_bigbins ${stateBedDir} ${prefix} ${state_num} ${suffix}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/state_bed/"
prefix="normal"
state_num=44
suffix=400
count_state_proportion_bigbins ${stateBedDir} ${prefix} ${state_num} ${suffix}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/state_bed/"
prefix="normal"
state_num=44
suffix=500
count_state_proportion_bigbins ${stateBedDir} ${prefix} ${state_num} ${suffix}

stateBedDir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/state_bed/"
prefix="normal"
state_num=44
suffix=2000
count_state_proportion_bigbins ${stateBedDir} ${prefix} ${state_num} ${suffix}

# ######################
# # redo regression using different state proportion
# ## get TSS promoter bins
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$0}else{for(i=1;i<=10;i++){x=10+i;print a[$5+i],$6"_"x;y=10-i;print a[$5-i],$6"_"(y+1)}}}' ocr_bins.bed mm10_protein_coding_TSS_ocr_bin.txt | sort -k5,5n > mm10_protein_coding_TSS_ocr_promoter_bins.txt
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$0}else{for(i=1;i<=10;i++){x=10+i;print a[$5+i],$6"_"x;y=10-i;print a[$5-i],$6"_"(y+1)}}}' normal_bins.bed mm10_protein_coding_TSS_ocr_bin.txt | sort -k5,5n > mm10_protein_coding_TSS_normal_promoter_bins.txt
# ## count state proportion
# count_state_proportion2(){
#     ocr_state_bed_dir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/state_bed/"
#     normal_state_bed_dir="/data/zusers/fankaili/ideas/dhs_bins/normal_bins/state_bed/"
#     #
#     while read line
#     do
#         echo $line
#         # ### get state in each window
#         # intersectBed -a mm10_protein_coding_TSS_ocr_promoter_bins.txt -b \
#         # ${ocr_state_bed_dir}embryonic-facial-prominence_11.5_state_sorted.bed -wo | sort -u | \
#         # awk '{FS=OFS="\t"}{split($6,a,"_");print $1,$2,$3,$4,$5,$6,a[1],a[2],$11}' \
#         # | sort -k7,7 -k8,8n > ./ocr_state_proportion2/mm10_gene_${line}_ocr_state.txt
#         # intersectBed -a mm10_protein_coding_TSS_normal_promoter_bins.txt -b \
#         # ${ocr_state_bed_dir}embryonic-facial-prominence_11.5_state_sorted.bed -wo | sort -u | \
#         # awk '{FS=OFS="\t"}{split($6,a,"_");print $1,$2,$3,$4,$5,$6,a[1],a[2],$11}' \
#         # | sort -k7,7 -k8,8n > ./normal_state_proportion2/mm10_gene_${line}_normal_state.txt
#         ### calculate state count in 20 windows
#         for i in {1..20}
#         do
#             awk -v n="$i" '{FS=OFS="\t"}{if($8==n){s=$9;printf $7;if(s!=0){if(s!=42){for(i=0;i<s;i++){printf "\t"0};printf "\t"1;for(i=s+1;i<=42;i++){printf "\t"0};printf "\n"}else{for(i=0;i<42;i++){printf "\t"0};printf "\t"1"\n"}}else{printf "\t"1;for(i=1;i<=42;i++){printf "\t"0};printf "\n"}}}' ./ocr_state_proportion2/mm10_gene_${line}_ocr_state.txt | sort -k1,1 > ./ocr_state_proportion2/mm10_gene_${line}_ocr_state_count_window_${i}.txt
#             #
#             awk -v n="$i" '{FS=OFS="\t"}{if($8==n){s=$9;printf $7;if(s!=0){if(s!=43){for(i=0;i<s;i++){printf "\t"0};printf "\t"1;for(i=s+1;i<=43;i++){printf "\t"0};printf "\n"}else{for(i=0;i<43;i++){printf "\t"0};printf "\t"1"\n"}}else{printf "\t"1;for(i=1;i<=43;i++){printf "\t"0};printf "\n"}}}' ./normal_state_proportion2/mm10_gene_${line}_normal_state.txt | sort -k1,1 > ./normal_state_proportion2/mm10_gene_${line}_normal_state_count_window_${i}.txt
#         done
#     done < /data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt
# }
#
# mkdir ocr_state_proportion2
# mkdir normal_state_proportion2
# count_state_proportion2


## Rscript do_within_gene_regression.R
