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
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs 0-0.5 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs 0.5-1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs 1-2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs 2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
#
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal 0-0.5 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal 0.5-1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal 1-2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal 2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/

### no log sd
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs 2_0 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs 2_0-1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs 2_1-2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs 2_2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/

#
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal 2_0 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal 2_0-1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal 2_1-2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal 2_2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/

# mean
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs mean_1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs mean_1-10 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs mean_10-50 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs mean_50 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
#
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal mean_1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal mean_1-10 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal mean_10-50 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal mean_50 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/

# median
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs median_1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs median_1-10 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs median_10-50 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs median_50 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
#
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal median_1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal median_1-10 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal median_10-50 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal median_50 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/

# p0.01
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p0.01_0-0.2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p0.01_0.2-0.5 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p0.01_0.5-1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p0.01_1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
#
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p0.01_0-0.2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p0.01_0.2-0.5 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p0.01_0.5-1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p0.01_1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/

# p0.1
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p0.1_0-0.2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p0.1_0.2-0.5 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p0.1_0.5-1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p0.1_1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
#
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p0.1_0-0.2 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p0.1_0.2-0.5 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p0.1_0.5-1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p0.1_1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/


# p1
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p1_0-0.1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p1_0.1-0.3 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p1_0.3-0.5 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py dhs p1_0.5 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
#
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p1_0-0.1 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p1_0.1-0.3 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p1_0.3-0.5 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/
python ${scriptDir}get_matrix_for_regression_across_celltype.py normal p1_0.5 \
/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/

# 3. do regression for each matrix
rm ./across_sample_matrix/r-square_across_celltype_*.txt
for i in {1..20}
do
    echo ${i}
    for range in p1_0.5
    do
        ## dhs
        Rscript ${scriptDir}calculate_rsqaure.R /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/regression_across_celltype_data_dhs_${range}_window${i}.txt dhs | awk '{print $2}' >> ./across_sample_matrix/r-square_across_celltype_dhs_${range}.txt
        ## normal
        Rscript ${scriptDir}calculate_rsqaure.R /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/regression_across_celltype_data_normal_${range}_window${i}.txt normal | awk '{print $2}' >> ./across_sample_matrix/r-square_across_celltype_normal_${range}.txt
    done
done

0-0.5 0.5-1 1-2 2
2_0 2_0-1 2_1-2 2_2
mean_1 mean_1-10 mean_10-50 mean_50
median_1 median_1-10 median_10-50 median_50
p0.01_0-0.2 p0.01_0.2-0.5 p0.01_0.5-1 p0.01_1
p0.1_0-0.2 p0.1_0.2-0.5 p0.1_0.5-1 p0.1_1
p1_0-0.1 p1_0.1-0.3 p1_0.3-0.5 p1_0.5

# 4. make figures
# Rscript do_across_gene_regression.R
