#!/bin/bash

# -- Kaili
# This script is for generating a table for given tissue.
# 1. get ±150bp window
# 2. get number of dELS in these regions in all tissues, as well as non-redundant ELSs
# 3. num&name of other active genes in the region
# 4. get gene expression in all e11.5 tissues
# 5. get UCSC link
# 6. get pubMed citation
# 7. get table

tissue=$1

workDir="/data/zusers/fankaili/regions/"
cd ${workDir}

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/regions_for_testing/script/"
ccreDir="/data/projects/encode/Registry/V2/mm10/"
signalDir="/data/projects/encode/Registry/V2/mm10/Signal-Files/"
rna_filelist="/data/zusers/fankaili/regions/mouse_e11.5_RNA-seq_list.txt"
expDir="/data/zusers/fankaili/regions/e11.5_expressed_genes/"

if [ ! -d /data/zusers/fankaili/regions/${tissue}/ ];then mkdir /data/zusers/fankaili/regions/${tissue}/;fi

# 1. get ±150bp window
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{s=$2-150000;e=$3+150000;if(s<0){s=0};if(e>a[$1]){e=a[$1]};print $1,s,e,"id",$5,$6,$7,$8}}' /home/fankaili/genome/mm10.chrom.sizes.clean ./e11.5_expressed_genes/e11.5_${tissue}_expressedGene_withPLS.txt | sort -u | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{print $1,$2,$3,"region_"NR,$5,$6,$7,$8}' > ./${tissue}/e11.5_${tissue}_center_regions.bed
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7}else{print $0,a[$7]}}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_withGeneSymbol.bed
cut -f 7,9 ./${tissue}/e11.5_${tissue}_center_regions_withGeneSymbol.bed | sort -u > ./${tissue}/e11.5_${tissue}_center_regions_geneList.txt
###-----
echo "**get 300kb region.**"

# 2. get number of dELS in these regions in all tissues, as well as non-redundant ELSs
awk 'BEGIN{FS=OFS="\t";print "chr\ts\te\tid\tscore\tstrand\tgene\texp"}{print $0}' ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_dELS.bed
#
awk -v sample="total_dELS" '{if(NR==1){print $0,sample}}' ./${tissue}/e11.5_${tissue}_center_regions_dELS.bed > ./${tissue}/tmp.bed
intersectBed -a ./${tissue}/e11.5_${tissue}_center_regions_dELS.bed -b ./e11.5_dELS/e11.5_dELS.bed -c >> ./${tissue}/tmp.bed
mv ./${tissue}/tmp.bed ./${tissue}/e11.5_${tissue}_center_regions_dELS.bed
#
while read line
do
    sample=`awk '{split($9,a,"_");print a[2]}' <<< $line`
    if [ "$sample" = "embryonic" ]; then sample="facial"; elif [ "$sample" = "neural" ]; then sample="neural-tube"; fi
    echo $sample
    #
    awk -v sample="$sample" '{if(NR==1){print $0,sample}}' ./${tissue}/e11.5_${tissue}_center_regions_dELS.bed > ./${tissue}/tmp.bed
    intersectBed -a ./${tissue}/e11.5_${tissue}_center_regions_dELS.bed -b ./e11.5_dELS/e11.5_${sample}_dELS.bed -c >> ./${tissue}/tmp.bed
    mv ./${tissue}/tmp.bed ./${tissue}/e11.5_${tissue}_center_regions_dELS.bed
done < e11.5_master_list.txt
###-----
echo "**get number of dELS in the region.**"

# 3. num&name of other active genes in the region
### num
intersectBed -a ./${tissue}/e11.5_${tissue}_center_regions.bed -b ./e11.5_expressed_genes/e11.5_${tissue}_expressedGene_withPLS.txt -wa -wb | awk '{FS=OFS="\t"}{if($7!=$15){print $4,$7,$15}}' | sort -u | awk 'BEGIN{FS=OFS="\t";id="";num=0}{if(NR==1){id=$1;num+=1}else{if(id==$1){num+=1}else{print id,num;id=$1;num=1}}}END{print id,num}' > ./${tissue}/e11.5_${tissue}_center_regions_otherActiveGenesNum.txt
### name
intersectBed -a ./${tissue}/e11.5_${tissue}_center_regions.bed -b ./e11.5_expressed_genes/e11.5_${tissue}_expressedGene_withPLS.txt -wa -wb | awk '{FS=OFS="\t"}{if($7!=$15){print $4,$7,$15}}' | sort -u | awk 'BEGIN{FS=OFS="\t";id="";gene=""}{if(NR==1){id=$1;gene=$3}else{if(id==$1){a=gene";"$3;gene=a}else{print id,gene;id=$1;gene=$3}}}END{print id,gene}' > ./${tissue}/e11.5_${tissue}_center_regions_otherActiveGenesID.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){b[$4]=$7}else{split($2,a,";");printf $1"\t"b[a[1]];if(length(a)>1){for(i=2;i<=length(a);i++){printf ";"b[a[i]]}};printf "\n"}}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt ./${tissue}/e11.5_${tissue}_center_regions_otherActiveGenesID.txt > ./${tissue}/e11.5_${tissue}_center_regions_otherActiveGenesSymbol.txt
###-----
echo "**get number and symbol of active gene in the region.**"

# 4. get gene expression in all e11.5 tissues
awk 'BEGIN{FS=OFS="\t";print "id","geneID"}{print $4,$7}' ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_exp.txt
#
while read line
do
    sample=`awk '{print $1}' <<< $line`
    expID=`awk '{print $2}' <<< $line`
    geneTsv=`awk '{print $3}' <<< $line`
    #
    awk -v sample="$sample" '{FS=OFS="\t"}{if(NR==FNR){split($1,c,".");a[c[1]]=1;b[c[1]]=$6}else{if(FNR==1){print $0,sample}else{split($2,d,".");if(a[d[1]]){print $0,b[d[1]]}}}}' /data/projects/encode/data/${expID}/${geneTsv}.tsv ./${tissue}/e11.5_${tissue}_center_regions_exp.txt > ./${tissue}/tmp.txt
    mv ./${tissue}/tmp.txt ./${tissue}/e11.5_${tissue}_center_regions_exp.txt
done < ${rna_filelist}
###-----
echo "**get the expression of center gene in all e11.5 tissues.**"

# 5. get UCSC link
if [ "$tissue" = "forebrain" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716617_WQaSVVUPPzYl5Scn7RpUGenLeQa8\n"}' ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_UCSCLink.bed
elif [ "$tissue" = "midbrain" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716047_C2IX7f2FAiUWvqN5N3kbKAlRFyNA\n"}' ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_UCSCLink.bed
elif [ "$tissue" = "hindbrain" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716067_2GcUjeKAeMp2kyYEmg9MzRiYRijZ\n"}' ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_UCSCLink.bed
elif [ "$tissue" = "neural-tube" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716085_FHZsFaKGlycmtG38h9cEM1EKZJR5\n"}' ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_UCSCLink.bed
elif [ "$tissue" = "heart" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716133_vgDF7r2M7wGaQqB4qDTouLZhxtHl\n"}' ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_UCSCLink.bed
elif [ "$tissue" = "facial" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716167_WvQpIgjp82QSEGndXZMMwA93oRhG\n"}' ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_UCSCLink.bed
elif [ "$tissue" = "limb" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716199_e645Pyy6KNhxjj3jn4OQpsDXRsCi\n"}' ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_UCSCLink.bed
elif [ "$tissue" = "liver" ];then
    awk '{FS=OFS="\t"}{printf $0;printf "\thttps://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A";printf("%d", $2);printf "-";printf("%d", $3);printf "&hgsid=724716225_8UFE2jh6bMpDsIHotyyK40G1dyE5\n"}' ./${tissue}/e11.5_${tissue}_center_regions.bed > ./${tissue}/e11.5_${tissue}_center_regions_UCSCLink.bed
fi
###-----
echo "**get UCSC links.**"


# 6. get pubMed citation
# bash /data/zusers/fankaili/github/weng-lab/Kaili/regions_for_testing/script/get_pubMed_citation.sh ./${tissue}/e11.5_${tissue}_center_regions_geneList.txt ./${tissue}/e11.5_${tissue}_center_regions_pubMed_citation.txt
# echo "**get pubMed citation.**"

# 7. get table
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$4]){print $1,$2,$3,$4,$7,$9,a[$4]}else{print $1,$2,$3,$4,$7,$9,0}}}' ./${tissue}/e11.5_${tissue}_center_regions_otherActiveGenesNum.txt ./${tissue}/e11.5_${tissue}_center_regions_withGeneSymbol.bed > ./${tissue}/tmp1.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$4]){print $0,a[$4]}else{print $0,"-"}}}' ./${tissue}/e11.5_${tissue}_center_regions_otherActiveGenesSymbol.txt ./${tissue}/tmp1.bed > ./${tissue}/tmp2.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}else{print $0,a[$4]}}' ./${tissue}/e11.5_${tissue}_center_regions_dELS.bed ./${tissue}/tmp2.bed > ./${tissue}/tmp3.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}else{print $0,a[$4]}}' ./${tissue}/e11.5_${tissue}_center_regions_exp.txt ./${tissue}/tmp3.bed > ./${tissue}/tmp4.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{print $0,a[$4]}}' ./${tissue}/e11.5_${tissue}_center_regions_UCSCLink.bed ./${tissue}/tmp4.bed > ./${tissue}/tmp5.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$3]=$4;b[$3]=$5;c[$3]=1}else{if(c[$5]){print $0,a[$5],b[$5]}else{print $0,"-","-"}}}' ./${tissue}/e11.5_${tissue}_center_regions_pubMed_citation.txt ./${tissue}/tmp5.bed | sort -k9nr > ./${tissue}/e11.5_${tissue}_center_regions_Matrix.txt
awk '{FS=OFS="\t"}{loci=$1":"$2"-"$3;a="=HYPERLINK(\""$26"\",\""loci"\")";if($27!="-"){b="=HYPERLINK(\""$28"\",\""$27"\")"}else{b="-"};print a,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,b}' ./${tissue}/e11.5_${tissue}_center_regions_Matrix.txt > ./${tissue}/e11.5_${tissue}_center_regions_Table1.txt
awk '{FS=OFS="\t"}{print $1,"300,000",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$23}' ./${tissue}/e11.5_${tissue}_center_regions_Table1.txt > ./${tissue}/e11.5_${tissue}_center_regions_Table2.txt
###-----
echo "**get Matrix and tables! Cheers!**"
rm ./${tissue}/tmp*.bed
