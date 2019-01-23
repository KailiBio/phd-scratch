#!/bin/bash

# -- Kaili
# This script is for proving ubi-rOCRs overlapped TSSs are more likely to form clusters.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

# 1. histogram: number of TSSs/ubi-rOCRs overlapped TSSs in gene
## number of TSSs in each gene
cut -f 7-8 TSS.Filtered.uniqID.bed | sort -u | cut -f 1 | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' \
> GRCh38_gene_allTSS_count.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$2}}}' \
/home/fankaili/genome/hg38_proteinCoding_geneID_geneType_geneSymbol.txt GRCh38_gene_allTSS_count.txt \
> GRCh38_coding-gene_allTSS_count.txt
## number of ubi-rOCRs overlapped gene in each gene
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$8]){print $7,$8}}}' \
GRCh38_ubi-rOCR_overlapped_TSS_uniqID.txt TSS.Filtered.uniqID.bed | sort -u | cut -f 1 | sort | uniq -c | \
awk '{OFS="\t"}{print $2,$1}' > GRCh38_gene_ubi-rOCR_TSS_count.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$2}}}' \
/home/fankaili/genome/hg38_proteinCoding_geneID_geneType_geneSymbol.txt GRCh38_gene_ubi-rOCR_TSS_count.txt \
> GRCh38_coding-gene_ubi-rOCR_TSS_count.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$1]){print $1,$2,a[$1],a[$1]/$2}else{print $1,$2,0,0}}}' \
GRCh38_coding-gene_ubi-rOCR_TSS_count.txt GRCh38_coding-gene_allTSS_count.txt > \
GRCh38_coding-gene_ubi-rOCR_TSS_percentage.txt
#Rscript ubi-rOCRs_overlapped_TSSs_to_clusters.R


# 2. histogram: number of TSSs overlapped in each rOCRs/ubi-rOCRs
intersectBed -a GRCh38-rOCRs.bed -b TSS.uniq.bed -wa -c > GRCh38_rOCRs_TSS_count.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5,"ubi-rOCRs"}else{print $4,$5,"non_ubi-rOCRs"}}}' \
GRCh38_ubi-rOCRs_EDGEid.bed GRCh38_rOCRs_TSS_count.txt > GRCh38_rOCRs_TSS_count_annotated.txt

# 3. number of nearby TSS for each TSS
# expand 50bp
awk '{FS=OFS="\t"}{print $1,$2-50,$3+50,$8,$5,$6,$7}' TSS.Filtered.uniq.bed | sort -k1,1 -k2,2n > \
TSS_up_down_50bp.bed
# count neary TSSs
intersectBed -a TSS_up_down_50bp.bed -b TSS.Filtered.uniq.bed -wa -c > TSS_neary_TSS_count.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{$8-=1;if(a[$4]){print $0,"ubi-rOCRs_overlapped_TSS"}else{print $0,"remaining_TSS"}}}' \
GRCh38_ubi-rOCR_overlapped_TSS_uniqID.txt TSS_neary_TSS_count.txt > TSS_neary_TSS_count_annotated.txt
####### for ubi-TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$8]){print $0}}}' \
GRCh38_ubi-rOCR_overlapped_TSS_uniqID.txt TSS.Filtered.uniq.bed > TSS_ubi-rOCRs_overlapped_uniq.bed
intersectBed -a TSS_up_down_50bp.bed -b TSS_ubi-rOCRs_overlapped_uniq.bed -wa -c > TSS_neary_ubi-TSS_count.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){if($8>0){$8-=1};print $0,"ubi-rOCRs_overlapped_TSS"}else{print $0,"remaining_TSS"}}}' \
GRCh38_ubi-rOCR_overlapped_TSS_uniqID.txt TSS_neary_ubi-TSS_count.txt > TSS_neary_ubi-TSS_count_annotated.txt
###### merge
cut -f 4,9,8 TSS_neary_TSS_count_annotated.txt | sort -u | sort -k1,1 > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$8}else{print $1,$3,$2,a[$1]}}' TSS_neary_ubi-TSS_count_annotated.txt \
tmp.txt > TSS_neary_TSS_ubi-TSS_count_merge.txt
