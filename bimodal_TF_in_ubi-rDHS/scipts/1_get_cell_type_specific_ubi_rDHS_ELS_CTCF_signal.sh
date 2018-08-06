#！/bin/bash

# -- Kaili
# This is script for getting the cell-type-specific ubi-rDHS PLS/ESL/CTCF-only/DNase/inactive.
# INPUT: ubi-rDHS: /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt
# OUTPUT: cell-type specific ubi-rDHS: /data/zusers/fankaili/ccre/tf/cell_type_specific/


scriptDir="/data/zusers/fankaili/github/KailiBio/bimodal_TF_in_ubi-rDHS/scripts/"

# 1. get cell type specific definition file list
python ${script}get_cell_type_specific_definition_file_list.py

# 2. divide ubi-rDHS
python ${script}divide_ubi_rDHS_by_cell_type_specific.py.py
#
nohup python /data/zusers/fankaili/github/KailiBio/bimodal_TF_in_ubi-rDHS/scripts/divide_ubi_rDHS_by_cell_type_specific.py > /data/zusers/fankaili/ccre/tf/cell_type_specific/nohup.divide_ubi_rDHS_by_cell_type_specific.py.out 2>&1&
