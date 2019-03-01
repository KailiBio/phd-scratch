#!/bin/bash

# -- Kaili
# This script is for validating the gene expression prediction for IDEAS result (across cell types).
### "Accurate and reproducible functional maps in 127 human cell types via 2D genome
### only do ±2kb here, no B splines.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/"

cd ${workDir}

# 1. get gene list with different sd
# Rscript do_across_gene_regression.R

# 2. get data matrix for each group
get_matrix_for_regression_across_celltype.py






# 1. get adjusted R-square for each gene
cut -f 1 mm10_RNA_protein-coding_tpm_matrix_matched.txt | awk 'NR>1' > gene_list.txt
## dhs bins: z008 11127
for i in {1..221}
do
    echo $i
    awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list.txt > gene_list_dhs.txt
    python ${scriptDir}get_gene_across_celltype_rsquare.py dhs /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list_dhs.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 1 > ./nohup/get_gene_across_celltype_rsquare_dhs_${i}.out
    cat gene_rsquare_matrix_dhs_1.txt >> gene_rsquare_matrix_dhs_z008.txt
done

nohup bash run_gene_regression_dhs_z008.sh > ./nohup/nohup.run_gene_regression_dhs_z008.out 2>&1&
## normal bins: z010 17416
for i in {1..221}
do
    echo $i
    awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list.txt > gene_list_normal.txt
    python ${scriptDir}get_gene_across_celltype_rsquare.py normal /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list_normal.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 1 > ./nohup/get_gene_across_celltype_rsquare_normal_${i}.out
    cat gene_rsquare_matrix_normal_1.txt >> gene_rsquare_matrix_normal_z010.txt
done

nohup bash run_gene_regression_normal_z010.sh > ./nohup/nohup.run_gene_regression_normal_z010.out 2>&1&


###########
# z008: 30367
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' gene_rsquare_matrix_dhs_z008.txt gene_list.txt > gene_list_dhs_00.txt
for i in {1..96}
do
    echo $i
    awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*110) && FNR<=(n*110)){print $0}}' gene_list_dhs_00.txt > gene_list_dhs.txt
    python ${scriptDir}get_gene_across_celltype_rsquare.py dhs /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list_dhs.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 1 > ./nohup/get_gene_across_celltype_rsquare_dhs_${i}.out
    cat gene_rsquare_matrix_dhs_1.txt >> gene_rsquare_matrix_dhs_z008.txt
done
nohup bash run_gene_regression_dhs_z008_2.sh > ./nohup/nohup.run_gene_regression_dhs_z008_2.out 2>&1&
# z010: 35197
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' gene_rsquare_matrix_normal_z010.txt gene_list.txt > gene_list_normal_00.txt
for i in {1..98}
do
    echo $i
    awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list_normal_00.txt > gene_list_normal.txt
    python ${scriptDir}get_gene_across_celltype_rsquare.py normal /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list_normal.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 1 > ./nohup/get_gene_across_celltype_rsquare_normal_${i}.out
    cat gene_rsquare_matrix_normal_1.txt >> gene_rsquare_matrix_normal_z010.txt
done
nohup bash run_gene_regression_normal_z010_2.sh > ./nohup/nohup.run_gene_regression_normal_z010_2.out 2>&1&
#
#
#
#
#
# nohup python ${scriptDir}get_gene_across_celltype_rsquare.py dhs /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ > ./nohup/get_gene_across_celltype_rsquare_dhs.out 2>&1&
# nohup python ${scriptDir}get_gene_across_celltype_rsquare.py normal /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ > ./nohup/get_gene_across_celltype_rsquare_normal.out 2>&1&
# #
# python ${scriptDir}get_gene_across_celltype_rsquare.py dhs ${workDir}ss.txt ${workDir} 0 > ${workDir}nohup/get_gene_across_celltype_rsquare_dhs_0.out 2>&1&
# #
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' gene_rsquare_matrix_dhs_0.txt gene_list.txt > gene_list2.txt
#
# #############
#
# awk '{if(NR>0 && NR<=5000){print $0}}' gene_list2.txt > gene_list2_1.txt
# awk '{if(NR>5000 && NR<=10000){print $0}}' gene_list2.txt > gene_list2_2.txt
# awk '{if(NR>10000 && NR<=15000){print $0}}' gene_list2.txt > gene_list2_3.txt
# awk '{if(NR>15000){print $0}}' gene_list2.txt > gene_list2_4.txt
#
# scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
# workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/"
#
# cd ${workDir}
# #------------------
# # normal
# # z008: 41680
# for i in {1..50}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list2_2.txt > gene_list2_z008.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py normal /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z008.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 2 > ./nohup/get_gene_across_celltype_rsquare_normal_z008_${i}.out
#     cat gene_rsquare_matrix_normal_2.txt >> gene_rsquare_matrix_normal_008.txt
# done
# # z010: 50616
# for i in {1..50}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list2_3.txt > gene_list2_z010.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py normal /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z010.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 3 > ./nohup/get_gene_across_celltype_rsquare_normal_z010_${i}.out
#     cat gene_rsquare_matrix_normal_3.txt >> gene_rsquare_matrix_normal_010.txt
# done
# # z010-2: 50657
# for i in {1..68}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list2_4.txt > gene_list2_z010_2.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py normal /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z010_2.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 4 > ./nohup/get_gene_across_celltype_rsquare_normal_z010_2_${i}.out
#     cat gene_rsquare_matrix_normal_4.txt >> gene_rsquare_matrix_normal_010_2.txt
# done
#
# cat gene_rsquare_matrix_normal_00.txt gene_rsquare_matrix_normal_001.txt gene_rsquare_matrix_normal_008.txt gene_rsquare_matrix_normal_010_2.txt gene_rsquare_matrix_normal_010.txt | sort -u > gene_rsquare_matrix_normal_000.txt
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' gene_rsquare_matrix_normal_000.txt gene_list.txt > gene_list2_05.txt
#
# # z008: 48374
# for i in {1..56}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list2_05.txt > gene_list2_z008_2.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py normal /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z008_2.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 1 > ./nohup/get_gene_across_celltype_rsquare_normal_z008_2_${i}.out
#     cat gene_rsquare_matrix_normal_1.txt >> gene_rsquare_matrix_normal_008_2.txt
# done
#
# cat gene_rsquare_matrix_normal_008_2.txt >> gene_rsquare_matrix_normal_000.txt
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' gene_rsquare_matrix_normal_000.txt gene_list.txt > gene_list2_06.txt
#
# # z008: 7161
# for i in {1..5}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list2_06.txt > gene_list2_z008_3.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py normal /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z008_3.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 1 > ./nohup/get_gene_across_celltype_rsquare_normal_z008_3_${i}.out
#     cat gene_rsquare_matrix_normal_1.txt >> gene_rsquare_matrix_normal_008_3.txt
# done
#
# cat gene_rsquare_matrix_normal_008_3.txt >> gene_rsquare_matrix_normal_000.txt
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' gene_rsquare_matrix_normal_000.txt gene_list.txt > gene_list2_07.txt
#
# #z008: 49205
# for i in {1..4}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list2_07.txt > gene_list2_z008_4.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py normal /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z008_4.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 1 > ./nohup/get_gene_across_celltype_rsquare_normal_z008_4_${i}.out
#     cat gene_rsquare_matrix_normal_1.txt >> gene_rsquare_matrix_normal_008_4.txt
# done
#
# cat gene_rsquare_matrix_normal_008_4.txt gene_rsquare_matrix_normal_000.txt | sort -k1,1 > gene_rsquare_matrix_normal.txt
# #------------------
# # dhs
# # z008: 51228
# for i in {1..50}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list2_1.txt > gene_list2_z008.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py dhs /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z008.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 2 > ./nohup/get_gene_across_celltype_rsquare_dhs_z008_${i}.out
#     cat gene_rsquare_matrix_dhs_2.txt >> gene_rsquare_matrix_dhs_008.txt
# done
#
# # z010: 16065
# for i in {1..50}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list2_2.txt > gene_list2_z010.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py dhs /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z010.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 3 > ./nohup/get_gene_across_celltype_rsquare_dhs_z010_${i}.out
#     cat gene_rsquare_matrix_dhs_3.txt >> gene_rsquare_matrix_dhs_010.txt
# done
#
# # z010-2: 16212
# for i in {1..68}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list2_4.txt > gene_list2_z010_2.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py dhs /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z010_2.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 4 > ./nohup/get_gene_across_celltype_rsquare_dhs_z010_2_${i}.out
#     cat gene_rsquare_matrix_dhs_4.txt >> gene_rsquare_matrix_dhs_010_2.txt
# done
#
# cat gene_rsquare_matrix_dhs_008.txt gene_rsquare_matrix_dhs_010.txt gene_rsquare_matrix_dhs_00.txt | sort -u > gene_rsquare_matrix_dhs_000.txt
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' gene_rsquare_matrix_dhs_000.txt gene_list.txt > gene_list_dhs_00.txt
# awk '{FS=OFS="\t"}{if(NR<=5500){print $0}}' gene_list_dhs_00.txt > gene_list_dhs_1.txt
# awk '{FS=OFS="\t"}{if(NR>5500){print $0}}' gene_list_dhs_00.txt > gene_list_dhs_2.txt
#
# # z008: 7016
# for i in {1..55}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list_dhs_1.txt > gene_list2_z008.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py dhs /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z008.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 2 > ./nohup/get_gene_across_celltype_rsquare_dhs_z008_2_${i}.out
#     cat gene_rsquare_matrix_dhs_2.txt >> gene_rsquare_matrix_dhs_008_2.txt
# done
#
# # z010: 57068
# for i in {1..63}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*100) && FNR<=(n*100)){print $0}}' gene_list_dhs_2.txt > gene_list2_z010.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py dhs /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z010.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 3 > ./nohup/get_gene_across_celltype_rsquare_dhs_z010_3_${i}.out
#     cat gene_rsquare_matrix_dhs_3.txt >> gene_rsquare_matrix_dhs_010_3.txt
# done
#
# cat gene_rsquare_matrix_dhs_008_2.txt gene_rsquare_matrix_dhs_010_3.txt >> gene_rsquare_matrix_dhs_000.txt
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' gene_rsquare_matrix_dhs_000.txt gene_list.txt > gene_list_dhs_01.txt
#
# # z010: 13268
# for i in {1..3}
# do
#     echo $i
#     awk -v n="$i" '{FS=OFS="\t"}{if(NR>((n-1)*110) && FNR<=(n*110)){print $0}}' gene_list_dhs_01.txt > gene_list2_z010.txt
#     python ${scriptDir}get_gene_across_celltype_rsquare.py dhs /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list2_z010.txt /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ 4 > ./nohup/get_gene_across_celltype_rsquare_dhs_z010_4_${i}.out
#     cat gene_rsquare_matrix_dhs_4.txt >> gene_rsquare_matrix_dhs_010_4.txt
# done
#
# cat gene_rsquare_matrix_dhs_010_4.txt gene_rsquare_matrix_dhs_000.txt | sort -k1,1 > gene_rsquare_matrix_dhs.txt


# 2. calculate gene expression sd, make figures

# Rscript do_across_gene_regression.R
