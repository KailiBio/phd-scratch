#!/bin/bash

# -- Kaili
# This script is for expression analysis of ubi-rDHS related genes.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scripts/"
workDir="/data/zusers/fankaili/ccre/tf/closest_gene/"

cd ${workDir}

# 1. get gene exp from ENCODE tsv file
## 183 RNA-seq in total
### get all 57820 genes
python ${scriptDir}get_RNA_exp.py
#
echo "gene_id" > all_gene_exp_matrix.txt
cut -f 1 ./all_gene_exp/ENCSR023ZXN.txt >> all_gene_exp_matrix.txt
for file in `ls /data/zusers/fankaili/ccre/tf/closest_gene/all_gene_exp/`
do
    echo $file;
    echo ${file%.txt} > temp.txt;
    awk '{print $2}' ./all_gene_exp/${file} >> temp.txt ;
    paste all_gene_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt all_gene_exp_matrix.txt;
done
rm temp.txt
### get 11188 closest genes
head -1 all_gene_exp_matrix.txt > closest_gene_exp_matrix.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$0}else{print a[$1]}}' all_gene_exp_matrix.txt hg19_ubi-rDHS_closest_gene_list0.txt >> closest_gene_exp_matrix.txt


# 2. gtf file processing
awk '{FS=" ";OFS="\t"}{if($3=="gene"){split($10,a,"\"");split($14,b,"\"");split($18,c,"\"");print a[2],b[2],c[2]}}' \
gencode.v19.annotation.gtf > hg19_geneID_geneType_geneSymbol.txt

## get protein-coding gene list
awk '{if($2=="protein_coding"){print $0}}' hg19_geneID_geneType_geneSymbol.txt > hg19_proteinCoding_geneID_geneType_geneSymbol.txt


Rscript make_closest_gene_exp_analysis_figs.R



#------------------------------------------------------------------------------------------------------------------------

# 3. get hg38 gene expression matrix
cd /data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/
python ${scriptDir}get_GRCh38_RNA_exp.py
#
echo "gene_id" > hg38_all_gene_exp_matrix.txt
cut -f 1 ./all_gene_exp/ENCSR023ZXN.txt >> hg38_all_gene_exp_matrix.txt
for file in `ls /data/zusers/fankaili/ccre/hg38_ubi-rDHS/all_gene_exp/`
do
    echo $file;
    echo ${file%.txt} > temp.txt;
    awk '{print $2}' ./all_gene_exp/${file} >> temp.txt ;
    paste hg38_all_gene_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt hg38_all_gene_exp_matrix.txt;
done
rm temp.txt

### get ref -- using GENCODE.v24
cp /data/zusers/moorej3/moorej.ghpcc.project/Reference/Human/hg38/GENCODE24/gencode.v24.annotation.gtf /home/fankaili/genome/
awk '{FS=" ";OFS="\t"}{if($3=="gene"){split($10,a,"\"");split($12,b,"\"");split($16,c,"\"");print a[2],b[2],c[2]}}' \
/home/fankaili/genome/gencode.v24.annotation.gtf > hg38_geneID_geneType_geneSymbol.txt
#
awk '{if($2=="protein_coding"){print $0}}' hg38_geneID_geneType_geneSymbol.txt > hg38_proteinCoding_geneID_geneType_geneSymbol.txt

# 4. get housekeeping genes and tissue-specific genes

## get housekeeping gene from paper
cd /home/fankaili/genome/
wget https://www.tau.ac.il/~elieis/HKG/HK_genes.txt
sed -i 's/ \t/\t/g' HK_genes.txt
# transfer RefSeq ID to GENCODE
# awk '{FS=OFS="\t"}{if($3=="mRNA"){split($9,a,";");split(a[4],b,"=");split(a[6],c,"=");split(b[2],d,".");print d[1],c[2]}}' \
# GRCh38_latest_genomic.gff | grep "NM_" > RefSeq_hg38_ID_symbol.txt
#
awk '{FS=OFS="\t"}{if($3=="mRNA"){print $0}}' GRCh38_latest_genomic.gff | cut -f 9 | \
awk '{FS=OFS="\t"}{split($0,a,";");for(i=1;i<=length(a);i++){split(a[i],b,"=");if(b[1]=="Name"){split(b[2],c,".");printf c[1]}else if(b[1]=="gene"){printf "\t"b[2]"\n"}}}' \
| grep "NM_" > RefSeq_hg38_ID_symbol.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$3]){print $0}}}' HK_genes.txt hg38_geneID_geneType_geneSymbol.txt \
> hg38_housekeeping_geneID_geneType_geneSymbol.txt
# rm "ENSG00000265264.5	protein_coding	TIMM10B"
# rm "ENSG00000282031.1	protein_coding	TMBIM4"
# rm "ENSG00000282795.1	pseudogene	UBE2NL"
cut -f 3 hg38_housekeeping_geneID_geneType_geneSymbol.txt > temp.txt
cut -f 1 HK_genes.txt >> temp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]!=1){print 0}}}' temp.txt temp1.txt | wc -l
sort temp.txt | uniq -u > temp2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' temp2.txt HK_genes.txt > temp3.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$2]){print $0,a[$2]}}}' RefSeq_hg38_ID_symbol.txt temp3.txt > temp4.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$3]=1}else{if(a[$3]){print $0}}}' temp4.txt hg38_geneID_geneType_geneSymbol.txt \
>> hg38_housekeeping_geneID_geneType_geneSymbol.txt
grep "ARRDC1-AS1" hg38_geneID_geneType_geneSymbol.txt >> hg38_housekeeping_geneID_geneType_geneSymbol.txt
grep "PPP4R3A" hg38_geneID_geneType_geneSymbol.txt >> hg38_housekeeping_geneID_geneType_geneSymbol.txt
grep "TMEM261" hg38_geneID_geneType_geneSymbol.txt >> hg38_housekeeping_geneID_geneType_geneSymbol.txt
grep "MACF1" hg38_geneID_geneType_geneSymbol.txt >> hg38_housekeeping_geneID_geneType_geneSymbol.txt
grep "PRKRIP1" hg38_geneID_geneType_geneSymbol.txt >> hg38_housekeeping_geneID_geneType_geneSymbol.txt
grep "RP11-74E24.2" hg38_geneID_geneType_geneSymbol.txt >> hg38_housekeeping_geneID_geneType_geneSymbol.txt
grep "SCNM1" hg38_geneID_geneType_geneSymbol.txt >> hg38_housekeeping_geneID_geneType_geneSymbol.txt
grep "AC073869.1" hg38_geneID_geneType_geneSymbol.txt >> hg38_housekeeping_geneID_geneType_geneSymbol.txt
#
sort -u hg38_housekeeping_geneID_geneType_geneSymbol.txt > temp5.txt
mv temp5.txt hg38_housekeeping_geneID_geneType_geneSymbol.txt
rm temp*.txt

## get tissue specificity score
### 1) get tissues expression matrix
# 46 tissues in total 113 samples
python ${scriptDir}get_GRCh38_tissue_RNA_exp.py
#
echo "gene_id" > hg38_tissue_gene_exp_matrix.txt
cut -f 1 ./all_gene_exp/ENCSR023ZXN.txt >> hg38_tissue_gene_exp_matrix.txt
for file in `ls /data/zusers/fankaili/ccre/hg38_ubi-rDHS/tissue_gene_exp/`
do
    echo $file;
    echo ${file%.txt} > temp.txt;
    awk '{print $2}' ./all_gene_exp/${file} >> temp.txt ;
    paste hg38_tissue_gene_exp_matrix.txt temp.txt > temp2.txt;
    mv temp2.txt hg38_tissue_gene_exp_matrix.txt;
done
rm temp.txt

### 2) calculate tissue specificity index
cp /data/zusers/zhangx/seq/mouse_ccRE/ts_all.py /data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scripts/
#
awk '{if(NR>1){print $0}}' hg38_tissue_gene_exp_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_tissue_gene_exp_TSscore.txt
rm tmp.txt
#
awk '{if(NR>1){print $0}}' hg38_all_gene_exp_matrix.txt > tmp.txt
python ${scriptDir}ts_all.py tmp.txt hg38_all_gene_exp_TSscore.txt
rm tmp.txt
