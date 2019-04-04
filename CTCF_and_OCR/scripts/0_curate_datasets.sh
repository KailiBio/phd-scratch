#!/bin/bash

# -- Kaili
# This script is for curating datasets.
# 1. ENCODE CTCF data list
# 2. ENCODE DNase/ATAC data list
# 3. Cistrome CTCF data list
# 4. Cistrome DNase/ATAC data list
# 5. ENCODE 3D data list
# 6. match all the data

workDir="/data/zusers/fankaili/CTCF_and_OCR/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/CTCF_and_OCR/scripts/"

cd ${workDir}

# 1. ENCODE CTCF data list
## 1) human
ENCODE_human_CTCF_url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=ChIP-seq&target.label=CTCF&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&limit=all&format=json"
python ${scriptDir}get_ENCODE_expID_info_from_URL.py ${ENCODE_human_CTCF_url} ${workDir}"encode_human_CTCF_explist.txt"
## 2) mouse
ENCODE_mouse_CTCF_url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=ChIP-seq&target.label=CTCF&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&limit=all&format=json"
python ${scriptDir}get_ENCODE_expID_info_from_URL.py ${ENCODE_mouse_CTCF_url} ${workDir}"encode_mouse_CTCF_explist.txt"

# 2. ENCODE DNase/ATAC data list
## 1) human
ENCODE_human_OCR_url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_slims=DNA+accessibility&assay_title=DNase-seq&assay_title=ATAC-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&limit=all&format=json"
python ${scriptDir}get_ENCODE_expID_info_from_URL.py ${ENCODE_human_OCR_url} ${workDir}"encode_human_OCR_explist.txt"
## 2) mouse
ENCODE_mouse_OCR_url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_slims=DNA+accessibility&assay_title=DNase-seq&assay_title=ATAC-seq&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&limit=all&format=json"
python ${scriptDir}get_ENCODE_expID_info_from_URL.py ${ENCODE_mouse_OCR_url} ${workDir}"encode_mouse_OCR_explist.txt"

# 3. Cistrome CTCF data list
## 1) human
awk '{OFS="\t"}{if($8=="CTCF"){print $2,$5,$6,$7,$1,$9}}' /data/projects/cistrome/metadata/TF_human_data_information.txt > Cistrome_human_CTCF_datalist0.txt
awk '{FS=OFS="\t"}{if($2!="None"){print $1,"ChIP-seq","CTCF",$2}else{print $1,"ChIP-seq","CTCF",$3}}' Cistrome_human_CTCF_datalist0.txt | sort -k4 > Cistrome_human_CTCF_datalist.txt
## 2) mouse
awk '{OFS="\t"}{if($8=="CTCF"){print $2,$5,$6,$7,$1,$9}}' /data/projects/cistrome/metadata/TF_mouse_data_information.txt > Cistrome_mouse_CTCF_datalist0.txt
awk '{FS=OFS="\t"}{if($2!="None"){print $1,"ChIP-seq","CTCF",$2}else{print $1,"ChIP-seq","CTCF",$3}}' Cistrome_mouse_CTCF_datalist0.txt | sort -k4 > Cistrome_mouse_CTCF_datalist.txt

# 4. Cistrome DNase/ATAC data list
## 1) human
cut -f 1-2,4-8 /data/projects/cistrome/metadata/ca_human_data_information.txt > Cistrome_human_OCR_datalist0.txt
awk '{FS=OFS="\t"}{if($3!="None"){print $2,$6,"-",$3}else{print $2,$6,"-",$4}}' Cistrome_human_OCR_datalist0.txt | sort -k4 > Cistrome_human_OCR_datalist.txt
## 2) mouse
cut -f 1-2,4-8 /data/projects/cistrome/metadata/ca_mouse_data_information.txt > Cistrome_mouse_OCR_datalist0.txt
awk '{FS=OFS="\t"}{if($3!="None"){print $2,$6,"-",$3}else{print $2,$6,"-",$4}}' Cistrome_mouse_OCR_datalist0.txt | sort -k4 > Cistrome_mouse_OCR_datalist.txt

# 5. ENCODE 3D data list
## 1) human
ENCODE_human_3D_url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_slims=3D+chromatin+structure&assay_title=Hi-C&assay_title=ChIA-PET&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&limit=all&format=json"
python ${scriptDir}get_ENCODE_expID_info_from_URL.py ${ENCODE_human_3D_url} ${workDir}"encode_human_3D_explist.txt"
## 2) mouse
ENCODE_mouse_3D_url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_slims=3D+chromatin+structure&assay_title=Hi-C&assay_title=ChIA-PET&limit=all&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&format=json"
python ${scriptDir}get_ENCODE_expID_info_from_URL.py ${ENCODE_mouse_3D_url} ${workDir}"encode_mouse_3D_explist.txt"

# 6. match all the data
python ${scriptDir}match_data_by_biosample.py
