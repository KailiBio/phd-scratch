#!/bin/bash

# -- Kaili
# This script is for getting TF file list.
# INPUT: /data/zusers/fankaili/ccre/tf/cell_type_specific/hg19_ccREs_cell_type_specific_definition_file_list.txt
# OUTPUT: file ID list in: /data/zusers/fankaili/ccre/tf/encode_tf_file_list/

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scipts/"
hg19_cell_type_specific_file="/data/zusers/fankaili/ccre/tf/cell_type_specific/hg19_ccREs_cell_type_specific_definition_file_list.txt"
outDir="/data/zusers/fankaili/ccre/tf/encode_tf_file_list/"
outName="encode_hg19_tf_cellline_with_cell_type_specific"


cd /data/zusers/fankaili/ccre/tf/encode_tf_file_list/

# 1. get cell lines
## 1-1  all celllines in ccRE cell-type specific file
awk '{FS=OFS="\t"}{print $1}' ${hg19_cell_type_specific_file} > hg19_ccRE_biosample.txt
awk '{FS=OFS="\t"}{if($1!~/_/){print $0}}' hg19_ccRE_biosample.txt > hg19_ccRE_cellline.txt

## 1-2 all cellines with hg19 TF ChIP-seq data in ENCODE
python ${scriptDir}get_biosample_list.py
sort -u encode_hg19_tf_cellline_list0.txt > encode_hg19_tf_cellline_list.txt

## 1-3 get repeated celllines
cat hg19_ccRE_cellline.txt encode_hg19_tf_cellline_list.txt | sort | uniq -d > encode_hg19_tf_cellline_with_cell_type_specific_list.txt

## remove astrocyte,H1-hESC,myotube,neutrophil,osteoblast from /data/zusers/fankaili/ccre/tf/hg19_cellline_with_CTCF.txt

# 2. get URL
awk '{FS=OFS="\t"}{print $1,"https://www.encodeproject.org/search/?type=Experiment&assay_title=ChIP-seq&target.investigated_as=transcription+factor&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assembly=hg19&biosample_term_name="$1"&assay_title=ChIP-seq&limit=all&format=json"}' \
${outName}_list.txt > ${outName}_URL.txt
sed -i 's/NT2\/D1\\t/NT2-D1\\t/g' ${outName}_URL.txt

# 3. get ENCODE hg19 TF list for each cellline
python ${scriptDir}get_encode_tf_file_list.py ${outDir}${outName}_URL.txt
