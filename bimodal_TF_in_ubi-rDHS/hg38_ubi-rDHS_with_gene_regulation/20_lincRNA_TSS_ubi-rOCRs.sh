#!/bin/bash

# -- Kaili
# This script is for analyzing lincRNA that whose TSS overlap ubi-rOCRs.
# 0. count rOCRs overlapping situation with genes


cd /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/

# 0. count rOCRs overlapping situation with genes
###################
# Apr 18
# count gene that whose TSS overlap rOCRs
###################
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$6}else{print $0,a[$7]}}' /home/fankaili/genome/hg38_v28_basic_gene_filtered.txt hg38_v28_basic_TSS_filtered_uniq.bed | sort -u > hg38_v28_basic_TSS_filtered_uniq_withGeneType.bed
cut -f 4 hg38_v28_basic_TSS_filtered_uniq_withGeneType.bed | sort | uniq -d > TSS_be_shared.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' TSS_be_shared.txt hg38_v28_basic_TSS_filtered_uniq_withGeneType.bed
####
# protein-coding TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $1,$2,$3,$8,$5,$6}}}' /home/fankaili/genome/hg38_v28_basic_TSS_protein_coding.bed hg38_v28_basic_TSS_filtered_with_uniqID.bed | sort -u > ./rOCRs_TSSs/TSS_protein-coding.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$8]){print $0}}}' ./rOCRs_TSSs/TSS_protein-coding.bed hg38_v28_basic_TSS_filtered_with_uniqID.bed | sort -u > ./rOCRs_TSSs/TSS_protein-coding_withGENCODEid.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{print $0,a[$4]}}' /home/fankaili/genome/hg38_v28_basic_transcript_filtered.txt ./rOCRs_TSSs/TSS_protein-coding_withGENCODEid.bed | sort -u | sort -k1,1 -k2,2n > ./rOCRs_TSSs/TSS_protein-coding_withGENCODEid_withGeneType.bed
# non-protein-coding TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $1,$2,$3,$4,$5,$6}}}' ./rOCRs_TSSs/TSS_protein-coding.bed hg38_v28_basic_TSS_filtered_uniq.bed | sort -u > ./rOCRs_TSSs/TSS_non-protein-coding.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$8]){print $0}}}' ./rOCRs_TSSs/TSS_non-protein-coding.bed hg38_v28_basic_TSS_filtered_with_uniqID.bed | sort -u > ./rOCRs_TSSs/TSS_non-protein-coding_withGENCODEid.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{print $0,a[$4]}}' /home/fankaili/genome/hg38_v28_basic_transcript_filtered.txt ./rOCRs_TSSs/TSS_non-protein-coding_withGENCODEid.bed | sort -u | sort -k1,1 -k2,2n > ./rOCRs_TSSs/TSS_non-protein-coding_withGENCODEid_withGeneType.bed
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' GRCh38_ubi-rOCRs_EDGEid.bed GRCh38-rOCRs.bed | sort -u > GRCh38-non-ubi-rOCRs.bed
#
grep "protein_coding" TSS_protein-coding_withGENCODEid_withGeneType.bed | cut -f 7 | sort -u > protein-coding_list.txt
cut -f 7 TSS_protein-coding_withGENCODEid_withGeneType.bed | sort -u > tmp.txt
cut -f 7 TSS_non-protein-coding_withGENCODEid_withGeneType.bed | sort -u > non-protein-coding_list.txt
#########
cd /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/rOCRs_TSSs/
#
intersectBed -a TSS_protein-coding_withGENCODEid_withGeneType.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb > protein-coding-TSS_ubi-rOCRs.bed
cut -f 8 protein-coding-TSS_ubi-rOCRs.bed | sort -u | wc -l
grep "protein_coding" protein-coding-TSS_ubi-rOCRs.bed | cut -f 7 | sort -u | wc -l
grep "protein_coding" protein-coding-TSS_ubi-rOCRs.bed | cut -f 7 | sort -u > protein_overlap_ubi-rOCRs.txt
#
intersectBed -a TSS_protein-coding_withGENCODEid_withGeneType.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed -wa -wb > protein-coding-TSS_non-ubi-rOCRs.bed
cut -f 8 protein-coding-TSS_non-ubi-rOCRs.bed | sort -u | wc -l
grep "protein_coding" protein-coding-TSS_non-ubi-rOCRs.bed | cut -f 7 | sort -u > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' protein_overlap_ubi-rOCRs.txt tmp.txt > protein_overlap_non-ubi-rOCRs.txt
#
intersectBed -a TSS_non-protein-coding_withGENCODEid_withGeneType.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb > non-protein-coding-TSS_ubi-rOCRs.bed
cut -f 8 non-protein-coding-TSS_ubi-rOCRs.bed | sort -u | wc -l
cut -f 7 non-protein-coding-TSS_ubi-rOCRs.bed | sort -u | wc -l
cut -f 7 non-protein-coding-TSS_ubi-rOCRs.bed | sort -u | awk '{if($1!="ENSG00000234912.11" && $1!="ENSG00000281398.3"){print $0}}'> non-protein-coding_overlap_ubi-rOCRs.txt
cut -f 7,9 non-protein-coding-TSS_ubi-rOCRs.bed | sort -u | cut -f 2 | sort | uniq -c
#
intersectBed -a TSS_non-protein-coding_withGENCODEid_withGeneType.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed -wa -wb > non-protein-coding-TSS_non-ubi-rOCRs.bed
cut -f 8 non-protein-coding-TSS_non-ubi-rOCRs.bed | sort -u | wc -l
cut -f 7 non-protein-coding-TSS_non-ubi-rOCRs.bed | sort -u > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' non-protein-coding_overlap_ubi-rOCRs.txt tmp.txt | awk '{if($1!="ENSG00000250802.7" && $1!="ENSG00000281005.1"){print $0}}'> non-protein_overlap_non-ubi-rOCRs.txt
##
intersectBed -a TSS_protein-coding_withGENCODEid_withGeneType.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -v | sort -u > tmp.bed
intersectBed -a tmp.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed -v | sort -u > tmp2.bed
cut -f 8 tmp2.bed | sort -u | wc -l
grep "protein_coding" tmp2.bed | cut -f 7 | sort -u > tmp2.txt
cat protein_overlap_ubi-rOCRs.txt tmp.txt protein_overlap_non-ubi-rOCRs.txt > protein_overlap_rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' protein_overlap_rOCRs.txt tmp2.txt | sort -u > protein_non-overlap_rOCRs.txt
#
intersectBed -a TSS_non-protein-coding_withGENCODEid_withGeneType.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -v | sort -u > tmp.bed
intersectBed -a tmp.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed -v | sort -u > tmp2.bed
cut -f 8 tmp2.bed | sort -u | wc -l
cut -f 7 tmp2.bed | sort -u > tmp2.txt
cat non-protein-coding_overlap_ubi-rOCRs.txt non-protein_overlap_non-ubi-rOCRs.txt > non-protein_overlap_rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $0}}}' non-protein_overlap_rOCRs.txt tmp2.txt | sort -u | awk '{if($1!="ENSG00000282057.1"){print $0}}'> non-protein_non-overlap_rOCRs.txt

###################
# Apr 07
# analysis on lincRNA that overlap ubi-rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$6}else{print $0,a[$1]}}' /home/fankaili/genome/hg38_v28_basic_gene_filtered.txt sample_RNA_signal_comparison.txt > tmp.txt
awk '{FS=OFS="\t"}{if($3=="genes_whose_TSSs_overlap_ubi-rOCRs" && $6=="protein_coding"){print $0,"ubi-PC"}else if($3=="genes_whose_TSSs_overlap_ubi-rOCRs" && $6=="lincRNA"){print $0,"ubi-lincRNA"}else if($3=="genes_whose_TSSs_overlap_other_active_rOCRs" && $6=="protein_coding"){print $0,"other-PC"}else if($3=="genes_whose_TSSs_overlap_other_active_rOCRs" && $6=="lincRNA"){print $0,"other-lincRNA"}}' tmp.txt > sample_RNA_signal_comparison_PC_linc.txt
rm tmp.txt
###########
# Apr 18
# lincRNA transcript level
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0,"protein_coding"}else{print $0,"non-coding"}}}' ./rOCRs_TSSs/protein-coding_list.txt sample_RNA_signal_comparison.txt sample_RNA_signal_comparison.txt > tmp.txt
grep "lincRNA" ./rOCRs_TSSs/TSS_non-protein-coding_withGENCODEid_withGeneType.bed | cut -f 7,9 > tmp2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if($6=="protein_coding"){print $0}else{if(a[$1]){print $1,$2,$3,$4,$5,"lincRNA"}}}}' tmp2.txt tmp.txt | sort -u > tmp3.txt
awk '{FS=OFS="\t"}{if($3=="genes_whose_TSSs_overlap_ubi-rOCRs" && $6=="protein_coding"){print $0,"ubi-PC"}else if($3=="genes_whose_TSSs_overlap_ubi-rOCRs" && $6=="lincRNA"){print $0,"ubi-lincRNA"}else if($3=="genes_whose_TSSs_overlap_other_active_rOCRs" && $6=="protein_coding"){print $0,"other-PC"}else if($3=="genes_whose_TSSs_overlap_other_active_rOCRs" && $6=="lincRNA"){print $0,"other-lincRNA"}}' tmp3.txt > sample_RNA_signal_comparison_PC_linc_new.txt
# lincRNA TSS level
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $0,"protein_coding"}else{print $0,"non-coding"}}}' ./rOCRs_TSSs/TSS_protein-coding.bed sample_RAMPAGE_signal_comparison.txt > tmp.txt
grep "lincRNA" ./rOCRs_TSSs/TSS_non-protein-coding_withGENCODEid_withGeneType.bed | cut -f 8,9 > tmp2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if($6=="protein_coding"){print $0}else{if(a[$1]){print $1,$2,$3,$4,$5,"lincRNA"}}}}' tmp2.txt tmp.txt | sort -u > tmp3.txt
awk '{FS=OFS="\t"}{if($3=="TSS_overlap_ubi-rOCRs" && $6=="protein_coding"){print $0,"ubi-PC"}else if($3=="TSS_overlap_ubi-rOCRs" && $6=="lincRNA"){print $0,"ubi-lincRNA"}else if($3=="TSS_overlap_non-ubi_active-rOCRs" && $6=="protein_coding"){print $0,"other-PC"}else if($3=="TSS_overlap_non-ubi_active-rOCRs" && $6=="lincRNA"){print $0,"other-lincRNA"}}' tmp3.txt > sample_RAMPAGE_signal_comparison_PC_linc_new.txt

# Rscript fig1.basic_v28.R

# 1. TSS-linc-ubi-rOCRs to TSS-PC-ubi-rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6;b[$1]=$7}else{if(a[$4]){print $0,a[$4],b[$4]}}}' sample_RAMPAGE_signal_comparison_PC_linc_new.txt hg38_v28_basic_TSS_filtered_uniq.bed > ./rOCRs_TSSs/PC_lincRNA_TSS.bed
grep "ubi-lincRNA" ./rOCRs_TSSs/PC_lincRNA_TSS.bed | sort -k1,1 -k2,2n > ./rOCRs_TSSs/ubi-lincRNA_TSS.bed
grep "ubi-PC" ./rOCRs_TSSs/PC_lincRNA_TSS.bed | sort -k1,1 -k2,2n > ./rOCRs_TSSs/ubi-PC_TSS.bed
#
grep "ubi-lincRNA" ./rOCRs_TSSs/PC_lincRNA_TSS.bed > tmp.bed
grep "ubi-PC" ./rOCRs_TSSs/PC_lincRNA_TSS.bed >> tmp.bed
sort -k1,1 -k2,2n tmp.bed > ./rOCRs_TSSs/ubi-PC_lincRNA_TSS.bed
intersectBed -a ./rOCRs_TSSs/ubi-PC_lincRNA_TSS.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb | sort -u | cut -f 1-4,9,13 | sort -k6,6 -k2,2n > ./rOCRs_TSSs/ubi-rOCRs_PC_lincRNA.txt
cut -f 6 ./rOCRs_TSSs/ubi-rOCRs_PC_lincRNA.txt | sort | uniq -c | head
#
intersectBed -a ./rOCRs_TSSs/ubi-lincRNA_TSS.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb | cut -f 1-4,10-13 | sort -u > ./rOCRs_TSSs/ubi-lincRNA_TSS_ubi-rOCRs.bed
intersectBed -a ./rOCRs_TSSs/PC_lincRNA_TSS.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb | cut -f 1-4,13 | sort -u
# group ubi-rOCRs with lincRNA
cut -f 8 ./rOCRs_TSSs/ubi-lincRNA_TSS_ubi-rOCRs.bed | sort -u > ./rOCRs_TSSs/ubi-rOCRs_with_lincRNA.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$6]){print $0}}}' ./rOCRs_TSSs/ubi-rOCRs_with_lincRNA.txt ./rOCRs_TSSs/ubi-rOCRs_PC_lincRNA.txt > ./rOCRs_TSSs/gene_in_ubi-rOCRs.txt
grep "ubi-PC" ./rOCRs_TSSs/gene_in_ubi-rOCRs.txt | cut -f 6 | sort -u > ./rOCRs_TSSs/ubi-rOCRs_with_both.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print $1}}}' ./rOCRs_TSSs/ubi-rOCRs_with_both.txt ./rOCRs_TSSs/ubi-rOCRs_with_lincRNA.txt > ./rOCRs_TSSs/ubi-rOCRs_only_lincRNA.txt
# count
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$6]){print $0}}}' ./rOCRs_TSSs/ubi-rOCRs_only_lincRNA.txt ./rOCRs_TSSs/gene_in_ubi-rOCRs.txt > ./rOCRs_TSSs/lincRNA_in_ubi-rOCRs_withoutPC.bed
cut -f 4 ./rOCRs_TSSs/lincRNA_in_ubi-rOCRs_withoutPC.bed | sort -u > lincRNA_in-linc-only_ubi-rOCRs.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$6]){print $0}}}' ./rOCRs_TSSs/ubi-rOCRs_with_both.txt ./rOCRs_TSSs/gene_in_ubi-rOCRs.txt > ./rOCRs_TSSs/lincRNA_in_ubi-rOCRs_withPC.bed
grep "ubi-lincRNA" ./rOCRs_TSSs/lincRNA_in_ubi-rOCRs_withPC.bed | cut -f 4 | sort -u > lincRNA_in-linc-PC_ubi-rOCRs.txt
# RAMPAGE signal
awk '{FS=FOS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0,"in_linc_only_ubi-rOCRs"}}}' lincRNA_in-linc-only_ubi-rOCRs.txt sample_RAMPAGE_signal_comparison_PC_linc_new.txt > sample_RAMPAGE_signal_comparison_linc_in_ubi-rOCRs.txt
awk '{FS=FOS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0,"in_linc_PC_ubi-rOCRs"}}}' lincRNA_in-linc-PC_ubi-rOCRs.txt sample_RAMPAGE_signal_comparison_PC_linc_new.txt >> sample_RAMPAGE_signal_comparison_linc_in_ubi-rOCRs.txt

#Rscript make_figures_lincRNA-ubi.R

# 2.

# ####################
# # Mar 25th
# # TSSs overlap with ubi-rOCRs and rOCRs
# # TSSs
# intersectBed -a hg38_v28_basic_TSS_protein_coding.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -u | cut -f 4 | sort -u | wc -l
# #
# intersectBed -a hg38_v28_basic_TSS_protein_coding.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -wa -u | cut -f 4 | sort -u | wc -l
# #
# intersectBed -a hg38_v28_basic_TSS_protein_coding.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -v | cut -f 4 | sort -u | wc -l
# ###
# intersectBed -a hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -u | cut -f 4 | sort -u | wc -l
# #
# intersectBed -a hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -wa -u | cut -f 4 | sort -u | wc -l
# #
# intersectBed -a hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -v | cut -f 4 | sort -u | wc -l
# # genes
# intersectBed -a hg38_v28_basic_TSS_protein_coding.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -u | cut -f 7 | sort -u | wc -l
# #
# intersectBed -a hg38_v28_basic_TSS_protein_coding.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -wa -u | cut -f 7 | sort -u | wc -l
# #
# intersectBed -a hg38_v28_basic_TSS_protein_coding.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -v | cut -f 7 | sort -u | wc -l
# ###
# intersectBed -a hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -u | cut -f 7 | sort -u | wc -l
# #
# intersectBed -a hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -wa -u | cut -f 7 | sort -u | wc -l
# #
# intersectBed -a hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -v | cut -f 7 | sort -u | wc -l
#
#
# intersectBed -a hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -u | cut -f 7 | sort -u > ss.txt
#
# awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' ss.txt hg38_v28_basic_gene_filtered.txt | head



##############
# Jun16
# re-calculate gene overlapping ubi-rOCRs, remove ubi-rOCRs that overlapping PC while calculating for non-ubi-rOCRs.
cd /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/ubi-rOCRs_overlap_geneType/
#
intersectBed -a /home/fankaili/genome/hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs_withLable.bed -wa -wb > all_TSS_rOCRs.txt
## ubi-rOCRs
awk '{if($12=="ubi-rOCR"){print $0}}' all_TSS_rOCRs.txt > all_TSS_ubi-rOCRs.txt
# 1. for PC_ubi
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0}}}' /home/fankaili/genome/hg38_v28_basic_TSS_protein_coding.bed all_TSS_ubi-rOCRs.txt > PC-TSS_ubi-rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' /home/fankaili/genome/hg38_v28_basic_TSS_protein_coding.bed all_TSS_ubi-rOCRs.txt > tmp1.txt
cut -f 7 PC-TSS_ubi-rOCRs.txt | sort -u > PC_ubi-rOCRs_genelist.txt
cut -f 11 PC-TSS_ubi-rOCRs.txt | sort -u > PC_ubi-rOCRs_rOCRlist.txt
# 2. for non-PC_ubi
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]!=1){print $0}}}' PC_ubi-rOCRs_genelist.txt tmp1.txt > tmp1.1.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$11]!=1){print $0}}}' PC_ubi-rOCRs_rOCRlist.txt tmp1.1.txt > nonPC-TSS_ubi-rOCRs.txt
cut -f 7 nonPC-TSS_ubi-rOCRs.txt | sort -u > nonPC_ubi-rOCRs_genelist.txt
cut -f 11 nonPC-TSS_ubi-rOCRs.txt | sort -u > nonPC_ubi-rOCRs_rOCRlist.txt
## non-ubi-rOCRs
awk '{if($12=="non-ubi-rOCR"){print $0}}' all_TSS_rOCRs.txt > tmp.all_TSS_non-ubi-rOCRs.txt
cat PC_ubi-rOCRs_genelist.txt nonPC_ubi-rOCRs_genelist.txt > tmp.genelist.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]!=1){print $0}}}' tmp.genelist.txt tmp.all_TSS_non-ubi-rOCRs.txt > all_TSS_non-ubi-rOCRs.txt
# 3. for PC_non-ubi
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0}}}' /home/fankaili/genome/hg38_v28_basic_TSS_protein_coding.bed all_TSS_non-ubi-rOCRs.txt > PC-TSS_non-ubi-rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1){print $0}}}' /home/fankaili/genome/hg38_v28_basic_TSS_protein_coding.bed all_TSS_non-ubi-rOCRs.txt > tmp2.txt
cut -f 7 PC-TSS_non-ubi-rOCRs.txt | sort -u > PC_non-ubi-rOCRs_genelist.txt
cut -f 11 PC-TSS_non-ubi-rOCRs.txt | sort -u > PC_non-ubi-rOCRs_rOCRlist.txt
# 4. for non-PC_non-ubi
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]!=1){print $0}}}' PC_non-ubi-rOCRs_genelist.txt tmp2.txt > tmp2.1.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$11]!=1){print $0}}}' PC_non-ubi-rOCRs_rOCRlist.txt tmp2.1.txt > nonPC-TSS_non-ubi-rOCRs.txt
cut -f 7 nonPC-TSS_non-ubi-rOCRs.txt | sort -u > nonPC_non-ubi-rOCRs_genelist.txt
cut -f 11 nonPC-TSS_non-ubi-rOCRs.txt | sort -u > nonPC_non-ubi-rOCRs_rOCRlist.txt
# gene without rOCRs
intersectBed -a /home/fankaili/genome/hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs_withLable.bed -v > tmp.TSS_no-rOCRs.txt
cat PC_ubi-rOCRs_genelist.txt nonPC_ubi-rOCRs_genelist.txt PC_non-ubi-rOCRs_genelist.txt nonPC_non-ubi-rOCRs_genelist.txt > tmp.genelist2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]!=1){print $0}}}' tmp.genelist2.txt tmp.TSS_no-rOCRs.txt > TSS_no-rOCRs.txt
##
awk '{FS=OFS="\t"}{if(NR==FNR){a[$7]=1}else{if(a[$7]){print $0}}}' /home/fankaili/genome/hg38_v28_basic_TSS_protein_coding.bed TSS_no-rOCRs.txt > PC-TSS_no-rOCRs.txt
cut -f 7 PC-TSS_no-rOCRs.txt | sort -u > PC-TSS_no-rOCRs_genelist.txt
##
awk '{FS=OFS="\t"}{if(NR==FNR){a[$7]=1}else{if(a[$7]!=1){print $0}}}' /home/fankaili/genome/hg38_v28_basic_TSS_protein_coding.bed TSS_no-rOCRs.txt > non-PC-TSS_no-rOCRs.txt
cut -f 7 non-PC-TSS_no-rOCRs.txt | sort -u > non-PC-TSS_no-rOCRs_genelist.txt




### PCgene
# PCgene & ubi-rOCRs
intersectBed -a /home/fankaili/genome/hg38_v28_basic_TSS_protein_coding.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb > PC-TSS_ubi-rOCRs.txt
cut -f 7 PC-TSS_ubi-rOCRs.txt | sort -u > PCgene_overlapping_ubi-rOCRs.txt
cut -f 11 PC-TSS_ubi-rOCRs.txt | sort -u > ubi-rOCRs_overlapping_PCgene-TSS.txt
# PCgene & non-ubi-rOCRs
intersectBed -a /home/fankaili/genome/hg38_v28_basic_TSS_protein_coding.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed -wa -wb > tmp.PC-TSS_non-ubi-rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]!=1){print $7}}}' PCgene_overlapping_ubi-rOCRs.txt tmp.PC-TSS_non-ubi-rOCRs.txt | sort -u > PC-TSS_non-ubi-rOCRs.txt
cut -f 11 tmp.PC-TSS_non-ubi-rOCRs.txt | sort -u > non-ubi-rOCRs_overlapping_PCgene-TSS.txt
# PCgene no rOCRs
intersectBed -a /home/fankaili/genome/hg38_v28_basic_TSS_protein_coding.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -v > tmp.PC-TSS_no-rOCRs.txt
cat PCgene_overlapping_ubi-rOCRs.txt PC-TSS_non-ubi-rOCRs.txt > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$7]!=1){print $7}}}' tmp.txt tmp.PC-TSS_no-rOCRs.txt | sort -u > PC-TSS_no-rOCRs.txt

### non-PCgene
# non-PCgene & ubi-rOCRs
intersectBed -a /home/fankaili/genome/hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb > tmp.non-PC-TSS_ubi-rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$7]!=1){print $0}}}' /home/fankaili/genome/hg38_v28_basic_gene_protein_coding.txt tmp.non-PC-TSS_ubi-rOCRs.txt > tmp.non-PC-TSS_ubi-rOCRs2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$11]!=1){print $0}}}' ubi-rOCRs_overlapping_PCgene-TSS.txt tmp.non-PC-TSS_ubi-rOCRs2.txt > non-PC-TSS_ubi-rOCRs.txt
cut -f 7 non-PC-TSS_ubi-rOCRs.txt | sort -u > non-PCgene_overlapping_ubi-rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' non-PCgene_overlapping_ubi-rOCRs.txt /home/fankaili/genome/hg38_v28_basic_gene_filtered.txt > non-PCgene_overlapping_ubi-rOCRs_withType.txt
# non-PCgene & non-ubi-rOCRs
intersectBed -a /home/fankaili/genome/hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-non-ubi-rOCRs.bed -wa -wb > tmp.non-PC-TSS_non-ubi-rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$7]!=1){print $0}}}' /home/fankaili/genome/hg38_v28_basic_gene_protein_coding.txt tmp.non-PC-TSS_non-ubi-rOCRs.txt > tmp.non-PC-TSS_non-ubi-rOCRs2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$11]!=1){print $0}}}' non-ubi-rOCRs_overlapping_PCgene-TSS.txt tmp.non-PC-TSS_non-ubi-rOCRs2.txt > non-PC-TSS_non-ubi-rOCRs.txt
cut -f 7 non-PC-TSS_non-ubi-rOCRs.txt | sort -u > non-PCgene_overlapping_non-ubi-rOCRs.txt
# non-PC no rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$7]!=1){print $0}}}' /home/fankaili/genome/hg38_v28_basic_gene_protein_coding.txt /home/fankaili/genome/hg38_v28_basic_TSS_filtered.bed > tmp.noPC.bed
intersectBed -a tmp.noPC.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -v > non-PC-TSS_no-rOCRs.txt
cut -f 7 non-PC-TSS_no-rOCRs.txt | sort -u > nonPC-TSS_no-rOCRs.txt



intersectBed -a /home/fankaili/genome/hg38_v28_basic_TSS_filtered.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -v | cut -f 7 | sort -u | wc -l
