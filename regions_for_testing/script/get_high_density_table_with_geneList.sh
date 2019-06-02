#!/bin/bash

# -- Kaili
# This script is for generating table of high enhancer density region with given gene file.
# 1. get gene desert region, with length
# 2. get number of dELSs in these regions
# 3. get UCSC links
# 4. make table

gene_file=$1
prefix=$2

mkdir /data/zusers/fankaili/regions/high_enhancer_density_without_${prefix}
cd /data/zusers/fankaili/regions/high_enhancer_density_without_${prefix}

# 1. get gene desert region, with length
## 1) get TSS±5kb regions.
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{if(a[$1]){s=int($2-5000);e=int($3+5000);if(s<0){s=0};if(e>a[$1]){e=a[$1]};printf $1"\t"s"\t"e"\n"}}}' /home/fankaili/genome/mm10.chrom.sizes.clean ${gene_file} | sort -k1,1 -k2,2n > tmp.bed
## 2) merged gene regions
bedtools merge -i tmp.bed > tmp2.bed
## 3) get gene desert regions
bedtools subtract -a /home/fankaili/genome/mm10.chrom.sizes.clean.bed -b tmp2.bed > tmp3.bed
bedtools subtract -a tmp3.bed -b /home/fankaili/genome/mm10.blacklist.bed | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{if(NR==FNR){print $0,"region_"NR, $3-$2}}' > ${prefix}_desert_regions.bed
###-----
echo "**get gene desert regions.**"

# 2. get number of dELSs in these regions
awk 'BEGIN{FS=OFS="\t";print "chr\ts\te\tid\tlength"}{print $0}' ${prefix}_desert_regions.bed > ${prefix}_desert_regions_dELS.bed
#
awk -v sample="total_dELS" '{if(NR==1){print $0,sample}}' ${prefix}_desert_regions_dELS.bed > tmp.bed
intersectBed -a ${prefix}_desert_regions_dELS.bed -b ../e11.5_dELS/e11.5_dELS.bed -c >> tmp.bed
mv tmp.bed ${prefix}_desert_regions_dELS.bed
#
while read line
do
    sample=`awk '{split($9,a,"_");print a[2]}' <<< $line`
    if [ "$sample" = "embryonic" ]; then sample="facial"; elif [ "$sample" = "neural" ]; then sample="neural-tube"; fi
    echo $sample
    #
    awk -v sample="$sample" '{if(NR==1){print $0,sample}}' ${prefix}_desert_regions_dELS.bed > tmp.bed
    intersectBed -a ${prefix}_desert_regions_dELS.bed -b ../e11.5_dELS/e11.5_${sample}_dELS.bed -c >> tmp.bed
    mv tmp.bed ${prefix}_desert_regions_dELS.bed
done < ../e11.5_master_list.txt
###-----
echo "**get number of dELS in the region.**"

# 3. get number of gene body
intersectBed -a ${prefix}_desert_regions.bed -b /home/fankaili/genome/mm10_vM18_basic_gene_filtered.bed -c > ${prefix}_desert_regions_NumOfGeneBody.bed

# 4. get UCSC links
awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=725069115_b5ud9s9radKDU8MlXIiwttqpAaiw\n"}' ${prefix}_desert_regions_dELS.bed > ${prefix}_desert_regions_dELS_UCSCLink.bed
###-----
echo "**get UCSC link.**"

# 5. make table
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$6}else{print $0,a[$4]}}' ${prefix}_desert_regions_NumOfGeneBody.bed ${prefix}_desert_regions_dELS_UCSCLink.bed > tmp_all.txt
awk '{FS=OFS="\t"}{if(NR>1){loci=$1":"$2"-"$3;link="=HYPERLINK(\""$15"\",\""loci"\")";if($6>0){print link,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16}}}' tmp_all.txt | sort -k4nr > ${prefix}_desert_regions_table.txt
###-----
echo "**get final table. Cheers!**"

rm tmp*
