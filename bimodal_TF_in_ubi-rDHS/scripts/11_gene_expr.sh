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

# 3. get RNAseq matched DNase sample
python ${scriptDir}get_RNAseq_match_DNase.py
sort -k3,3 hg19_RNA_matched_DNase_file_list.txt > hg19_RNA_matched_DNase_file_list_sorted.txt

### get not ubi-rDHS overlapped TSS file
awk '{FS=OFS="\t"}{if($13==0){print $5,$6,$7,$8,$9,$10,$11}}' hg19_ubi-rDHS_closest_gene.bed | sort -u | cut -f 4 \
> hg19_ubi-rDHS_overlapped_TSS_list.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]!=1){print $0}}}' hg19_ubi-rDHS_overlapped_TSS_list.txt TSS.Filtered.bed \
> hg19_ubi-rDHS_non_overlapped_TSS_list.txt

### get non ubi-rDHSs
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]!=1){print $0}}}' /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt \
/data/zusers/moorej3/ENCODE-Registry/hg19/V4/hg19-rDHSs.bed > /data/zusers/fankaili/ccre/hg19_non_ubi-rDHS.bed

#
while read line
do
    echo ${line};
    rna_exp_id=`awk '{FS=OFS="\t"}{print $1}' <<< ${line}` ;
    sample=`awk '{FS=OFS="\t"}{print $3}' <<< ${line}` ;
    dnase_exp_id=`awk '{FS=OFS="\t"}{print $5}' <<< ${line}` ;
    dnase_file_id=`awk '{FS=OFS="\t"}{print $6}' <<< ${line}` ;
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){if($2>1.64){a[$1]=1}}else{if(a[$4]){print $0}}}' \
    /data/zusers/moorej3/ENCODE-Registry/hg19/V4/signal-output/${dnase_exp_id}"-"${dnase_file_id}.txt \
    /data/zusers/fankaili/ccre/hg19_non_ubi-rDHS.bed > ./non_ubi_active_DHS/${dnase_exp_id}_DHS.bed ;
    #
    bedtools intersect -a ./non_ubi_active_DHS/${dnase_exp_id}_DHS.bed -b hg19_ubi-rDHS_non_overlapped_TSS_list.txt -wa -wb \
    | cut -f 11 | sort -u > temp.txt ;
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"ubi-rDHS_overlapped"}}}' \
    hg19_ubi-rDHS_overlapped_gene_list0.txt ./all_gene_exp/${rna_exp_id}.txt > ./gene_exp_comparison_file/${sample}.txt ;
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]==1){print $1,$2,"active-rDHS_overlapped"}}}' \
    temp.txt ./all_gene_exp/${rna_exp_id}.txt >> ./gene_exp_comparison_file/${sample}.txt ;
    #
    Rscript /data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scripts/make_comparison_barplot.R ${sample} ;
done < hg19_RNA_matched_DNase_file_list_sorted.txt

while read line
do
    sample=`awk '{FS=OFS="\t"}{print $3}' <<< ${line}` ;
    Rscript /data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scripts/make_comparison_barplot.R ${sample} ;
done < hg19_RNA_matched_DNase_file_list_sorted.txt


# 4. RAMGAE signal

###### RAMPAGE
head -1 hg19-tss-rampage-matrix.txt | awk '{FS=OFS="\t"}{for(i=1;i<=NF;i++){if($i=="ENCFF198YEH"){print i}}}'
# 239 +
head -1 hg19-tss-rampage-matrix.txt | awk '{FS=OFS="\t"}{for(i=1;i<=NF;i++){if($i=="ENCFF707TAV"){print i}}}'
# 251 -

bedtools intersect -a ./non_ubi_active_DHS/ENCSR000EMT_DHS.bed -b hg19_ubi-rDHS_non_overlapped_TSS_list.txt -wa -wb \
| cut -f 5,6,7,8,9,10 | sort -u > temp.txt ;

awk '{FS=OFS="\t"}{print $4,0,$6,"active-rDHS_overlapped"}' temp.txt > GM12878_rampage_file.txt
awk '{FS=OFS="\t"}{if($13==0){print $8,0,$10,"ubi-rDHS_overlapped"}}' hg19_ubi-rDHS_closest_gene.bed >> GM12878_rampage_file.txt

awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$239}else{if($3=="+" && a[$1]==1){print $1,b[$1],$3,$4}else{print $0}}}' \
hg19-tss-rampage-matrix.txt GM12878_rampage_file.txt > temp2.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$251}else{if($3=="-" && a[$1]==1){print $1,b[$1],$3,$4}else{print $0}}}' \
hg19-tss-rampage-matrix.txt temp2.txt > GM12878_rampage_file.txt

#------------------------------------------------------------------------------------------------------------------------

# 1. get housekeeping genes and tissue-specific genes

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
