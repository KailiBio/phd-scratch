#!/bin/bash

# -- Kaili
# This script is for generating table of high enhancer density region beyond gene.

workDir="/data/zusers/fankaili/regions/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/regions_for_testing/script/"

rna_filelist="/data/zusers/fankaili/regions/mouse_e11.5_RNA-seq_list.txt"

cd ${workDir}


# 1. all gene gene
bash ${scriptDir}get_high_density_table_with_geneList.sh "/home/fankaili/genome/mm10_vM18_basic_TSS_filtered.bed" "all_gene"

# 2. all protein-coding gene
bash ${scriptDir}get_high_density_table_with_geneList.sh "/home/fankaili/genome/mm10_vM18_basic_TSS_protein_coding.bed" "PC_gene"

# 3. e11.5 active gene
## 1) get e11.5 active gene - TPM>1
if [ -f ./e11.5_expressed_genes/e11.5_expressed_gene_TPM1_TSS.bed ];then rm ./e11.5_expressed_genes/e11.5_expressed_gene_TPM1_TSS.bed;fi
#
while read line
do
    sample=`awk '{print $1}' <<< $line`
    expID=`awk '{print $2}' <<< $line`
    geneTsv=`awk '{print $3}' <<< $line`
    echo ${sample}
    # get active genes
    grep "ENSMUSG" /data/projects/encode/data/${expID}/${geneTsv}.tsv | awk '{FS=OFS="\t"}{if($6>1){print $1,$2,$6}}' > ./e11.5_expressed_genes/e11.5_${sample}_expressed_gene_TPM1.txt
    # get TSS loci of active genes
    awk '{FS=OFS="\t"}{if(NR==FNR){split($1,c,".");a[c[1]]=1;b[c[1]]=$3}else{split($7,d,".");if(a[d[1]]){print $0,b[d[1]]}}}' ./e11.5_expressed_genes/e11.5_${sample}_expressed_gene_TPM1.txt /home/fankaili/genome/mm10_vM18_basic_TSS_filtered.bed | sort -k1,1 -k2,2n | awk '{if(NR>1){print $0}}' > ./e11.5_expressed_genes/e11.5_${sample}_expressed_gene_TPM1_TSS.txt
    cat ./e11.5_expressed_genes/e11.5_${sample}_expressed_gene_TPM1_TSS.txt >> ./e11.5_expressed_genes/e11.5_expressed_gene_TPM1_TSS.bed
done < ${rna_filelist}
#
cut -f 1-3 ./e11.5_expressed_genes/e11.5_expressed_gene_TPM1_TSS.bed | sort -u | sort -k1,1 -k2,2n > ./e11.5_expressed_genes/e11.5_expressed_gene_TPM1_TSS_sorted.bed
## 2) run pipeline
bash ${scriptDir}get_high_density_table_with_geneList.sh "/data/zusers/fankaili/regions/e11.5_expressed_genes/e11.5_expressed_gene_TPM1_TSS_sorted.bed" "e11.5_active_gene"
## 3) number of inactive genes
intersectBed -a ./high_enhancer_density_without_e11.5_active_gene/e11.5_active_gene_desert_regions_dELS_UCSCLink.bed -b /home/fankaili/genome/mm10_vM18_basic_TSS_filtered.bed -wa -wb | cut -f 1-4,22 | sort -u | cut -f 1-4 | sort | uniq -c | awk '{OFS="\t"}{print $5,$1}' > ./high_enhancer_density_without_e11.5_active_gene/e11.5_active_gene_desert_regions_inactiveGene_count.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$6}else{print $0,a[$4]}}' ./high_enhancer_density_without_e11.5_active_gene/e11.5_active_gene_desert_regions_NumOfGeneBody.bed ./high_enhancer_density_without_e11.5_active_gene/e11.5_active_gene_desert_regions_dELS_UCSCLink.bed > ./high_enhancer_density_without_e11.5_active_gene/tmp_all.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{if(a[$4]){print $0,a[$4]}else{print $0,0}}}' ./high_enhancer_density_without_e11.5_active_gene/e11.5_active_gene_desert_regions_inactiveGene_count.txt ./high_enhancer_density_without_e11.5_active_gene/tmp_all.txt > ./high_enhancer_density_without_e11.5_active_gene/e11.5_active_gene_desert_regions_Matrix.txt
#
awk '{FS=OFS="\t"}{if(NR>1){loci=$1":"$2"-"$3;link="=HYPERLINK(\""$15"\",\""loci"\")";if($6>0){print link,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16,$17}}}' ./high_enhancer_density_without_e11.5_active_gene/e11.5_active_gene_desert_regions_Matrix.txt | sort -k4nr > ./high_enhancer_density_without_e11.5_active_gene/e11.5_active_gene_desert_regions_table2.txt
