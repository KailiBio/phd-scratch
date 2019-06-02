#!/bin/bash

# -- Kaili
# This script is for generating table for merged regions of all tissues together.

workDir="/data/zusers/fankaili/regions/all_tissue/"
cd ${workDir}
rna_filelist="/data/zusers/fankaili/regions/mouse_e11.5_RNA-seq_list.txt"

# 1. get merged region across all 8 tissues
if [ -f e11.5_all_tissue_merged_region.bed ];then rm e11.5_all_tissue_merged_region.bed;fi
#
while read line
do
    sample=`awk '{print $1}' <<< $line`
    echo $sample
    #
    cat ../${sample}/e11.5_${sample}_merged_regions.bed >> e11.5_all_tissue_merged_region.bed
done < ${rna_filelist}
#
sort -k1,1 -k2,2n e11.5_all_tissue_merged_region.bed > tmp.bed
bedtools merge -i tmp.bed | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{print $0,"merged_region_"NR,$3-$2}' > e11.5_all_tissue_merged_region.bed
rm tmp.bed

# 2. get number of dELS in these regions in all tissues, as well as non-redundant ELSs
awk 'BEGIN{FS=OFS="\t";print "chr\ts\te\tid\tlength"}{print $0}' e11.5_all_tissue_merged_region.bed > e11.5_all_tissue_merged_regions_dELS.bed
#
awk -v sample="total_ELS" '{if(NR==1){print $0,sample}}' e11.5_all_tissue_merged_regions_dELS.bed > tmp.bed
intersectBed -a e11.5_all_tissue_merged_regions_dELS.bed -b ../e11.5_dELS/e11.5_dELS.bed -c >> tmp.bed
mv tmp.bed e11.5_all_tissue_merged_regions_dELS.bed
#
while read line
do
    sample=`awk '{split($9,a,"_");print a[2]}' <<< $line`
    if [ "$sample" = "embryonic" ]; then sample="facial"; elif [ "$sample" = "neural" ]; then sample="neural-tube"; fi
    echo $sample
    #
    awk -v sample="$sample" '{if(NR==1){print $0,sample}}' e11.5_all_tissue_merged_regions_dELS.bed > tmp.bed
    intersectBed -a e11.5_all_tissue_merged_regions_dELS.bed -b ../e11.5_dELS/e11.5_${sample}_dELS.bed -c >> tmp.bed
    mv tmp.bed e11.5_all_tissue_merged_regions_dELS.bed
done < ../e11.5_master_list.txt
###-----
echo "**get number of dELS in the region.**"

# 3. num&name of other active genes in the region
if [ -f ../e11.5_expressed_genes/e11.5_all_tissue_expressedGene_withPLS.txt ];then rm ../e11.5_expressed_genes/e11.5_all_tissue_expressedGene_withPLS.txt; fi
#
while read line
do
    sample=`awk '{print $1}' <<< $line`
    echo $sample
    #
    cat  ../e11.5_expressed_genes/e11.5_${sample}_expressedGene_withPLS.txt >> ../e11.5_expressed_genes/e11.5_all_tissue_expressedGene_withPLS.txt
done < ${rna_filelist}
#
cut -f 1-7 ../e11.5_expressed_genes/e11.5_all_tissue_expressedGene_withPLS.txt | sort -u > ../e11.5_expressed_genes/e11.5_all_tissue_expressedGene_withPLS_uniq.txt
### num
intersectBed -a e11.5_all_tissue_merged_region.bed -b ../e11.5_expressed_genes/e11.5_all_tissue_expressedGene_withPLS_uniq.txt -wa -wb | awk '{FS=OFS="\t"}{print $4,$12}' | sort -u | awk 'BEGIN{FS=OFS="\t";id="";num=0}{if(NR==1){id=$1;num+=1}else{if(id==$1){num+=1}else{print id,num;id=$1;num=1}}}END{print id,num}' > e11.5_all_tissue_merged_regions_otherActiveGenesNum.txt
### name
intersectBed -a e11.5_all_tissue_merged_region.bed -b ../e11.5_expressed_genes/e11.5_all_tissue_expressedGene_withPLS_uniq.txt -wa -wb | awk '{FS=OFS="\t"}{print $4,$12}' | sort -u | awk 'BEGIN{FS=OFS="\t";id="";gene=""}{if(NR==1){id=$1;gene=$2}else{if(id==$1){a=gene";"$2;gene=a}else{print id,gene;id=$1;gene=$2}}}END{print id,gene}' > e11.5_all_tissue_merged_regions_otherActiveGenesID.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){b[$4]=$7}else{split($2,a,";");printf $1"\t"b[a[1]];if(length(a)>1){for(i=2;i<=length(a);i++){printf ";"b[a[i]]}};printf "\n"}}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt e11.5_all_tissue_merged_regions_otherActiveGenesID.txt > e11.5_all_tissue_merged_regions_otherActiveGenesSymbol.txt
###-----
echo "**get number and symbol of active gene in the region.**"

# 4. get UCSC link
## 1) make merged ccre tracks
if [ -f /data/public_html_users/fankaili/ccre/e11.5_ccre_all.bed ];then rm /data/public_html_users/fankaili/ccre/e11.5_ccre_all.bed;fi
#
while read line
do
    sample=`awk '{print $1}' <<< $line`
    echo $sample
    #
    #bigBedToBed /data/public_html_users/fankaili/ccre/${sample}.7group.bigBed /data/public_html_users/fankaili/ccre/${sample}.7group.bed
    cat /data/public_html_users/fankaili/ccre/${sample}.7group.bed >> /data/public_html_users/fankaili/ccre/e11.5_ccre_all.bed
done < ${rna_filelist}
#
cd /data/public_html_users/fankaili/ccre/
sort -u e11.5_ccre_all.bed > tmp1.bed
cut -f 4 tmp1.bed | sort | uniq -c > tmp2.bed
awk '{if($1==1){print $2}}' tmp2.bed > tmp3.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' tmp3.bed tmp1.bed > e11.5_ccre_uniq.bed
awk '{if($1>1){print $2}}' tmp2.bed > tmp4.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' tmp4.bed tmp1.bed > tmp5.bed
## for ambigous one, pELS,dELS > PLS > DNase-H3K4me3 > DNase-only > Low-DNase
awk '{FS=OFS="\t"}{if($10=="pELS" || $10=="dELS"){print $0}}' tmp5.bed >> e11.5_ccre_uniq.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1 && $10=="PLS"){print $0}}}' e11.5_ccre_uniq.bed tmp5.bed > tmp6.bed
cat tmp6.bed >> e11.5_ccre_uniq.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1 && $10=="DNase-H3K4me3"){print $0}}}' e11.5_ccre_uniq.bed tmp5.bed > tmp7.bed
cat tmp7.bed >> e11.5_ccre_uniq.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]!=1 && $10=="DNase-only"){print $0}}}' e11.5_ccre_uniq.bed tmp5.bed > tmp8.bed
cat tmp8.bed >> e11.5_ccre_uniq.bed
#
sort -k1,1 -k2,2n e11.5_ccre_uniq.bed | cut -f 1-9 > e11.5_ccre_uniq_sorted.bed
bedToBigBed e11.5_ccre_uniq_sorted.bed /home/fankaili/genome/mm10.chrom.sizes.clean e11.5_ccre_uniq.bigBed
rm tmp*.bed
## 2) get links
cd ${workDir}
awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=725069115_b5ud9s9radKDU8MlXIiwttqpAaiw\n"}' e11.5_all_tissue_merged_region.bed > e11.5_all_tissue_merged_regions_UCSCLink.bed
###-----
echo "**get UCSC links.**"

# 5. get table
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{if(FNR>1){print $1,$2,$3,$4,$5,a[$4],$6,$7,$8,$9,$10,$11,$12,$13,$14}}}' e11.5_all_tissue_merged_regions_otherActiveGenesNum.txt e11.5_all_tissue_merged_regions_dELS.bed > tmp1.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$3,$4,$5,$6,a[$4],$7,$8,$9,$10,$11,$12,$13,$14,$15}}' e11.5_all_tissue_merged_regions_otherActiveGenesSymbol.txt tmp1.bed > tmp2.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$6}else{print $0,a[$4]}}' e11.5_all_tissue_merged_regions_UCSCLink.bed tmp2.bed > e11.5_all_tissue_merged_regions_Matrix.txt
#
sort -k 8nr e11.5_all_tissue_merged_regions_Matrix.txt | awk '{FS=OFS="\t"}{loci=$1":"$2"-"$3;a="=HYPERLINK(\""$17"\",\""loci"\")";print a,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > e11.5_all_tissue_merged_regions_Table.txt
###-----
echo "**get Matrix and tables! Cheers!**"
rm tmp*.bed
