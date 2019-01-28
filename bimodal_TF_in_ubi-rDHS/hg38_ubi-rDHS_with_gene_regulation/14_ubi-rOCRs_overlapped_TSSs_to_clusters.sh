#!/bin/bash

# -- Kaili
# This script is for proving ubi-rOCRs overlapped TSSs are more likely to form clusters.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

#Rscript ubi-rOCRs_overlapped_TSSs_to_clusters.R
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
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$1]){print $1,$2,a[$1],a[$1]/$2}else{print $1,$2,0,0}}}' \
GRCh38_gene_ubi-rOCR_TSS_count.txt GRCh38_gene_allTSS_count.txt > \
GRCh38_gene_ubi-rOCR_TSS_percentage.txt

cut -f 2-3 GRCh38_coding-gene_ubi-rOCR_TSS_percentage.txt | sort | uniq -c | awk '{OFS="\t"}{print $2,$3,$1}' |
sort -k3,3nr > GRCh38_coding-gene_ubi-rOCR_TSS_percentage_sort_count.txt

# 2. histogram: number of TSSs overlapped in each rOCRs/ubi-rOCRs
intersectBed -a GRCh38-rOCRs.bed -b TSS.uniq.bed -wa -c > GRCh38_rOCRs_TSS_count.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $4,$5,"ubi-rOCRs"}else{print $4,$5,"non_ubi-rOCRs"}}}' \
GRCh38_ubi-rOCRs_EDGEid.bed GRCh38_rOCRs_TSS_count.txt > GRCh38_rOCRs_TSS_count_annotated.txt
## for clusters: ubi-rOCRs overlapped TSSs percentage in each gene
cut -f 7 ./merged-TSS/GRCh38_ubi-rOCR_overlapped_merged-TSS.bed | sort | uniq -c | \
awk '{FS=" ";OFS="\t"}{print $2,$1}' > ./closest_gene/hg38_gene_overlapped_TSS_count_merged2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{if(a[$1]){print $0,a[$1]}else{print $0,0}}}' \
./closest_gene/hg38_gene_overlapped_TSS_count_merged2.txt GRCh38_coding-gene_ubi-rOCR_TSS_percentage.txt \
> GRCh38_coding-gene_ubi-rOCR_TSS_percentage_with_cluster.txt
# num of TSS in a gene vs singular&cluster
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$7]){if($3==$2){print $7,$4,a[$7],"singular"}else{print $7,$4,a[$7],"cluster"}}}}' \
GRCh38_coding-gene_ubi-rOCR_TSS_percentage.txt ./merged-TSS/GRCh38_ubi-rOCR_overlapped_merged-TSS.bed \
> GRCh38_coding-gene_numOfTSS_overlappedTSSLength.txt
# gene that with only one TSS overlapped ubi-rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR){if($5==1){a[$1]=1}}else{if(a[$1]){print $0}}}' \
GRCh38_coding-gene_ubi-rOCR_TSS_percentage_with_cluster.txt GRCh38_coding-gene_numOfTSS_overlappedTSSLength.txt \
> GRCh38_coding-gene_numOfTSS1_overlappedTSSLength.txt
# gene that with only two TSS overlapped ubi-rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR){if($5==2){a[$1]=1}}else{if(a[$1]){print $0}}}' \
GRCh38_coding-gene_ubi-rOCR_TSS_percentage_with_cluster.txt GRCh38_coding-gene_numOfTSS_overlappedTSSLength.txt \
| sort -k1,1 > GRCh38_coding-gene_numOfTSS2_overlappedTSSLength.txt
awk 'BEGIN{FS=OFS="\t";gene="";type="";num=""}{if(gene!=$1){gene=$1;num=$3;type=$4}else{if(type==$4){print gene,num,type}else{print gene,num,"half"}}}' \
GRCh38_coding-gene_numOfTSS2_overlappedTSSLength.txt > GRCh38_coding-gene_numOfTSS2_overlappedTSSLength_new.txt

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


# 4. distance to nearby TSSs
bedtools closest -a TSS.Filtered.uniq_sorted.bed -b TSS.Filtered.uniq_sorted.bed -io -d | \
awk '{FS=OFS="\t"}{print $8,$7,$16,$15,$17}' > TSS_nearest_TSS.bed
cut -f 1,5 TSS_nearest_TSS.bed | sort -u > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,$2,"ubi-rOCRs_overlapped_TSS"}else{print $1,$2,"remaining_TSS"}}}' \
GRCh38_ubi-rOCR_overlapped_TSS_uniqID.txt tmp.txt > TSS_nearest_TSS_forHist_annotated.bed


## 5. distance to nearby TSS in the same gene
sort -k7,7 -k2,2n TSS.Filtered.uniq.bed > TSS.Filtered.uniq_sortedByGene.bed
awk 'BEGIN{FS=OFS="\t";gene="";tss="";end="";distance=1000000}{if(gene==""){gene=$7;tss=$8;end=$3}else{if(gene==$7){if(distance>($2-end)){print gene,tss,$2-end}else{print gene,tss,distance};tss=$8;distance=$2-end;end=$3}else{print gene,tss,distance;gene=$7;tss=$8;end=$3;distance=1000000}}}END{print gene,tss,distance}' \
TSS.Filtered.uniq_sortedByGene.bed > TSS_nearest_TSS_distance_sameGene.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$2]){print $0,"ubi-rOCRs_overlapped_TSS"}else{print $0,"remaining_TSS"}}}' \
GRCh38_ubi-rOCR_overlapped_TSS_uniqID.txt TSS_nearest_TSS_distance_sameGene.bed > \
TSS_nearest_TSS_distance_sameGene_annotated.bed
