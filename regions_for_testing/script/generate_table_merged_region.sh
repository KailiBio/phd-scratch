#!/bin/bash

# -- Kaili
# This script is for generating table for merged regions.

tissue=$1

workDir="/data/zusers/fankaili/regions/"
cd ${workDir}${tissue}


# 1. get merged region
sort -k1,1 -k2,2n e11.5_${tissue}_center_regions_Matrix.txt > tmp.bed
bedtools merge -i tmp.bed | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{print $0,"merged_region_"NR,$3-$2}' > e11.5_${tissue}_merged_regions.bed
###-----
echo "**get merged regions.**"

# 2. get number of dELS in these regions in all tissues, as well as non-redundant ELSs
awk 'BEGIN{FS=OFS="\t";print "chr\ts\te\tid\tlength"}{print $0}' e11.5_${tissue}_merged_regions.bed > e11.5_${tissue}_merged_regions_dELS.bed
#
awk -v sample="total_ELS" '{if(NR==1){print $0,sample}}' e11.5_${tissue}_merged_regions_dELS.bed > tmp.bed
intersectBed -a e11.5_${tissue}_merged_regions_dELS.bed -b ../e11.5_ELS/e11.5_ELS.bed -c >> tmp.bed
mv tmp.bed e11.5_${tissue}_merged_regions_dELS.bed
#
while read line
do
    sample=`awk '{split($9,a,"_");print a[2]}' <<< $line`
    if [ "$sample" = "embryonic" ]; then sample="facial"; elif [ "$sample" = "neural" ]; then sample="neural-tube"; fi
    echo $sample
    #
    awk -v sample="$sample" '{if(NR==1){print $0,sample}}' e11.5_${tissue}_merged_regions_dELS.bed > tmp.bed
    intersectBed -a e11.5_${tissue}_merged_regions_dELS.bed -b ../e11.5_ELS/e11.5_${sample}_ELS.bed -c >> tmp.bed
    mv tmp.bed e11.5_${tissue}_merged_regions_dELS.bed
done < ../e11.5_master_list.txt
###-----
echo "**get number of dELS in the region.**"

# 3. num&name of other active genes in the region
### num
intersectBed -a e11.5_${tissue}_merged_regions.bed -b ../e11.5_expressed_genes/e11.5_${tissue}_expressedGene_withPLS.txt -wa -wb | awk '{FS=OFS="\t"}{print $4,$12}' | sort -u | awk 'BEGIN{FS=OFS="\t";id="";num=0}{if(NR==1){id=$1;num+=1}else{if(id==$1){num+=1}else{print id,num;id=$1;num=1}}}END{print id,num}' > e11.5_${tissue}_merged_regions_otherActiveGenesNum.txt
### name
intersectBed -a e11.5_${tissue}_merged_regions.bed -b ../e11.5_expressed_genes/e11.5_${tissue}_expressedGene_withPLS.txt -wa -wb | awk '{FS=OFS="\t"}{print $4,$12}' | sort -u | awk 'BEGIN{FS=OFS="\t";id="";gene=""}{if(NR==1){id=$1;gene=$2}else{if(id==$1){a=gene";"$2;gene=a}else{print id,gene;id=$1;gene=$2}}}END{print id,gene}' > e11.5_${tissue}_merged_regions_otherActiveGenesID.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){b[$4]=$7}else{split($2,a,";");printf $1"\t"b[a[1]];if(length(a)>1){for(i=2;i<=length(a);i++){printf ";"b[a[i]]}};printf "\n"}}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt e11.5_${tissue}_merged_regions_otherActiveGenesID.txt > e11.5_${tissue}_merged_regions_otherActiveGenesSymbol.txt
###-----
echo "**get number and symbol of active gene in the region.**"

# 4. get UCSC link
if [ "$tissue" = "forebrain" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716617_WQaSVVUPPzYl5Scn7RpUGenLeQa8\n"}' e11.5_${tissue}_merged_regions.bed > e11.5_${tissue}_merged_regions_UCSCLink.bed
elif [ "$tissue" = "midbrain" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716047_C2IX7f2FAiUWvqN5N3kbKAlRFyNA\n"}' e11.5_${tissue}_merged_regions.bed > e11.5_${tissue}_merged_regions_UCSCLink.bed
elif [ "$tissue" = "hindbrain" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716067_2GcUjeKAeMp2kyYEmg9MzRiYRijZ\n"}' e11.5_${tissue}_merged_regions.bed > e11.5_${tissue}_merged_regions_UCSCLink.bed
elif [ "$tissue" = "neural-tube" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716085_FHZsFaKGlycmtG38h9cEM1EKZJR5\n"}' e11.5_${tissue}_merged_regions.bed > e11.5_${tissue}_merged_regions_UCSCLink.bed
elif [ "$tissue" = "heart" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716133_vgDF7r2M7wGaQqB4qDTouLZhxtHl\n"}' e11.5_${tissue}_merged_regions.bed > e11.5_${tissue}_merged_regions_UCSCLink.bed
elif [ "$tissue" = "facial" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716167_WvQpIgjp82QSEGndXZMMwA93oRhG\n"}' e11.5_${tissue}_merged_regions.bed > e11.5_${tissue}_merged_regions_UCSCLink.bed
elif [ "$tissue" = "limb" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716199_e645Pyy6KNhxjj3jn4OQpsDXRsCi\n"}' e11.5_${tissue}_merged_regions.bed > e11.5_${tissue}_merged_regions_UCSCLink.bed
elif [ "$tissue" = "liver" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716225_8UFE2jh6bMpDsIHotyyK40G1dyE5\n"}' e11.5_${tissue}_merged_regions.bed > e11.5_${tissue}_merged_regions_UCSCLink.bed
fi
###-----
echo "**get UCSC links.**"

# 5. get table
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{if(FNR>1){print $1,$2,$3,$4,$5,a[$4],$6,$7,$8,$9,$10,$11,$12,$13,$14}}}' e11.5_${tissue}_merged_regions_otherActiveGenesNum.txt e11.5_${tissue}_merged_regions_dELS.bed > tmp1.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $1,$2,$3,$4,$5,$6,a[$4],$7,$8,$9,$10,$11,$12,$13,$14,$15}}' e11.5_${tissue}_merged_regions_otherActiveGenesSymbol.txt tmp1.bed > tmp2.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$6}else{print $0,a[$4]}}' e11.5_${tissue}_merged_regions_UCSCLink.bed tmp2.bed > e11.5_${tissue}_merged_regions_Matrix.txt
#
sort -k 8nr e11.5_${tissue}_merged_regions_Matrix.txt | awk '{FS=OFS="\t"}{loci=$1":"$2"-"$3;a="=HYPERLINK(\""$17"\",\""loci"\")";print a,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > e11.5_${tissue}_merged_regions_Table.txt
###-----
echo "**get Matrix and tables! Cheers!**"
rm tmp*.bed
