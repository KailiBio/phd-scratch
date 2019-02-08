#!/bin/bash

# -- Kaili
# This script is for validating the gene expression prediction for IDEAS result (across cell types).
### "Accurate and reproducible functional maps in 127 human cell types via 2D genome
### only do ±2kb here, no B splines.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/"

cd ${workDir}

# 1. get adjusted R-square for each gene
for i in {1..20}
do
    awk -v i="$i" '{if(NR>(1100*(i-1)) && NR<=(1100*i)){print $1}}' /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/dhs_state_proportion/mm10_gene_kidney_15.5_dhs_state_count_window_1.txt > gene_${i}.txt
    nohup bash ${scriptDir}get_gene_across_celltype_rsquare.sh ${i} > ./nohup/nohup.gene_across_${i}.out 2>&1&
done

cat gene_rsquare_matrix_dhs_*.txt > gene_rsquare_matrix_dhs.txt
cat gene_rsquare_matrix_normal_*.txt > gene_rsquare_matrix_normal.txt
rm gene_rsquare_matrix_dhs_*.txt gene_rsquare_matrix_normal_*.txt

nohup python ${scriptDir}get_gene_across_celltype_rsquare.py dhs /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ > ./nohup/nohup.get_gene_across_celltype_rsquare_dhs.out 2>&1&
nohup python ${scriptDir}get_gene_across_celltype_rsquare.py normal /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/ > ./nohup/nohup.get_gene_across_celltype_rsquare_normal.out 2>&1&
rm tmp_rsquare_dhs.txt tmp_rsquare_normal.txt

# 2. calculate gene expression sd, make figures
