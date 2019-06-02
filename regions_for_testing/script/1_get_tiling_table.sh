#!/bin/bash

# -- Kaili
# This script is for getting tiling table. (e11.5)

# 1. get e11.5 expressed gene in all 8 tissues
# 2. for forebrain, get ± 150 bp window
# 3. get number of dELS in these regions in all tissues
# 4. count other active genes in these regions
# 5. get expression in other 7 e11.5 tissues
# 6. get UCSC session
# 7. merge all, make tables

workDir="/data/zusers/fankaili/regions/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/regions_for_testing/script/"

ccreDir="/data/projects/encode/Registry/V2/mm10/"
signalDir="/data/projects/encode/Registry/V2/mm10/Signal-Files/"
seven_groupDir="/data/projects/encode/Registry/V2/mm10/Cell-Type-Specific/Seven-Group/"
rna_filelist="/data/zusers/fankaili/regions/mouse_e11.5_RNA-seq_list.txt"

cd ${workDir}

# 0. data preparation
## get master file in e11.5
grep "C57BL/6" ${ccreDir}Cell-Type-Specific/Master-Cell-List.txt | grep "embryo_11.5_days" > e11.5_master_list.txt
## mouse GENCODE M4 GTF
bash /data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/process_GENCODE_gtf.sh /home/fankaili/genome/gencode.vM4.annotation.gtf mm10 vM4 comprehensive
## mouse GENCODE M18 GTF
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.basic.annotation.gtf.gz
#
bash /data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/process_GENCODE_gtf.sh /home/fankaili/genome/gencode.vM18.basic.annotation.gtf mm10 vM18 basic

# 1. get e11.5 expressed gene with PLS in all 8 tissues
mkdir e11.5_expressed_genes
## 1) filter gene by TPM>5, get TSS loci
while read line
do
    sample=`awk '{print $1}' <<< $line`
    expID=`awk '{print $2}' <<< $line`
    geneTsv=`awk '{print $3}' <<< $line`
    # get active genes
    grep "ENSMUSG" /data/projects/encode/data/${expID}/${geneTsv}.tsv | awk '{FS=OFS="\t"}{if($6>5){print $1,$2,$6}}' > ./e11.5_expressed_genes/e11.5_${sample}_expressed_gene.txt
    # get TSS loci of active genes
    awk '{FS=OFS="\t"}{if(NR==FNR){split($1,c,".");a[c[1]]=1;b[c[1]]=$3}else{split($7,d,".");if(a[d[1]]){print $0,b[d[1]]}}}' ./e11.5_expressed_genes/e11.5_${sample}_expressed_gene.txt /home/fankaili/genome/mm10_vM18_basic_TSS_filtered.bed | sort -k1,1 -k2,2n | awk '{if(NR>1){print $0}}' > ./e11.5_expressed_genes/e11.5_${sample}_expressed_gene_TSS.txt
done < ${rna_filelist}
## 2) get PLS in all tissues
mkdir e11.5_PLS
# while read line
# do
#     dnase_exp=`awk '{print $1}' <<< $line`
#     dnase_file=`awk '{print $2}' <<< $line`
#     h3k4me3_exp=`awk '{print $3}' <<< $line`
#     h3k4me3_file=`awk '{print $4}' <<< $line`
#     h3k27ac_exp=`awk '{print $5}' <<< $line`
#     h3k27ac_file=`awk '{print $6}' <<< $line`
#     sample=`awk '{split($9,a,"_");print a[2]}' <<< $line`
#     if [ "$sample" = "embryonic" ]; then sample="facial"; elif [ "$sample" = "neural" ]; then sample="neural-tube"; fi
#     echo $sample
#     #
#     awk '{if(NR==FNR){if($2>1.64){a[$1]=1}}else{if($2>1.64 && a[$1]){print $1}}}' ${signalDir}${dnase_exp}-${dnase_file}.txt ${signalDir}${h3k4me3_exp}-${h3k4me3_file}.txt | sort -u > ./e11.5_PLS/e11.5_${sample}_PLS.txt
#     awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' ./e11.5_PLS/e11.5_${sample}_PLS.txt ${ccreDir}mm10-rDHSs.bed | sort -k1,1 -k2,2n > ./e11.5_PLS/e11.5_${sample}_PLS.bed
# done < e11.5_master_list.txt
while read line
do
    dnase_exp=`awk '{print $1}' <<< $line`
    dnase_file=`awk '{print $2}' <<< $line`
    h3k4me3_exp=`awk '{print $3}' <<< $line`
    h3k4me3_file=`awk '{print $4}' <<< $line`
    h3k27ac_exp=`awk '{print $5}' <<< $line`
    h3k27ac_file=`awk '{print $6}' <<< $line`
    sample=`awk '{split($9,a,"_");print a[2]}' <<< $line`
    if [ "$sample" = "embryonic" ]; then sample="facial"; elif [ "$sample" = "neural" ]; then sample="neural-tube"; fi
    echo $sample
    #
    awk '{FS=OFS="\t"}{if($10=="PLS"){print $1,$2,$3,$4}}' ${seven_groupDir}${dnase_file}_${h3k4me3_file}_${h3k27ac_file}.7group.bed > ./e11.5_PLS/e11.5_${sample}_PLS.bed
done < e11.5_master_list.txt

## 3) filter active genes by PLS
while read line
do
    sample=`awk '{print $1}' <<< $line`
    expID=`awk '{print $2}' <<< $line`
    geneTsv=`awk '{print $3}' <<< $line`
    #
    intersectBed -a ./e11.5_expressed_genes/e11.5_${sample}_expressed_gene_TSS.txt -b ./e11.5_PLS/e11.5_${sample}_PLS.bed -wa -wb > ./e11.5_expressed_genes/e11.5_${sample}_expressedGene_withPLS.txt
done < ${rna_filelist}


# 2. for forebrain, get ± 150 bp window
awk '{FS=OFS="\t"}{print $1,$2-150000,$3+150000,"id",$5,$6,$7,$8}' ./e11.5_expressed_genes/e11.5_forebrain_expressedGene_withPLS.txt | sort -u | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{print $1,$2,$3,"region_"NR,$5,$6,$7,$8}' > e11.5_forebrain_center_regions.bed

# 3. get number of dELS in these regions in all tissues
mkdir e11.5_ELS
#
awk 'BEGIN{FS=OFS="\t";print "chr\ts\te\tid\tscore\tstrand\tgene\texp"}{print $0}' e11.5_forebrain_center_regions.bed > e11.5_forebrain_center_regions_dELS.bed
#
# while read line
# do
#     dnase_exp=`awk '{print $1}' <<< $line`
#     dnase_file=`awk '{print $2}' <<< $line`
#     h3k4me3_exp=`awk '{print $3}' <<< $line`
#     h3k4me3_file=`awk '{print $4}' <<< $line`
#     h3k27ac_exp=`awk '{print $5}' <<< $line`
#     h3k27ac_file=`awk '{print $6}' <<< $line`
#     sample=`awk '{split($9,a,"_");print a[2]}' <<< $line`
#     if [ "$sample" = "embryonic" ]; then sample="facial"; elif [ "$sample" = "neural" ]; then sample="neural-tube"; fi
#     echo $sample
#     #
#     awk '{if(NR==FNR){if($2>1.64){a[$1]=1}}else{if($2>1.64 && a[$1]){print $1}}}' ${signalDir}${dnase_exp}-${dnase_file}.txt ${signalDir}${h3k27ac_exp}-${h3k27ac_file}.txt | sort -u > ./e11.5_ELS/e11.5_${sample}_ELS.txt
#     awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' ./e11.5_ELS/e11.5_${sample}_ELS.txt ${ccreDir}mm10-rDHSs.bed | sort -k1,1 -k2,2n > ./e11.5_ELS/e11.5_${sample}_ELS.bed
#     #
#     awk -v sample="$sample" '{if(NR==1){print $0,sample}}' e11.5_forebrain_center_regions_dELS.bed > tmp.bed
#     intersectBed -a e11.5_forebrain_center_regions_dELS.bed -b ./e11.5_ELS/e11.5_${sample}_ELS.bed -c >> tmp.bed
#     mv tmp.bed e11.5_forebrain_center_regions_dELS.bed
# done < e11.5_master_list.txt
mkdir e11.5_dELS
while read line
do
    dnase_exp=`awk '{print $1}' <<< $line`
    dnase_file=`awk '{print $2}' <<< $line`
    h3k4me3_exp=`awk '{print $3}' <<< $line`
    h3k4me3_file=`awk '{print $4}' <<< $line`
    h3k27ac_exp=`awk '{print $5}' <<< $line`
    h3k27ac_file=`awk '{print $6}' <<< $line`
    sample=`awk '{split($9,a,"_");print a[2]}' <<< $line`
    if [ "$sample" = "embryonic" ]; then sample="facial"; elif [ "$sample" = "neural" ]; then sample="neural-tube"; fi
    echo $sample
    #
    awk '{FS=OFS="\t"}{if($10=="dELS"){print $1,$2,$3,$4}}' ${seven_groupDir}${dnase_file}_${h3k4me3_file}_${h3k27ac_file}.7group.bed > ./e11.5_dELS/e11.5_${sample}_dELS.bed
done < e11.5_master_list.txt

# 4. count other active genes in these regions
intersectBed -a e11.5_forebrain_center_regions.bed -b ./e11.5_expressed_genes/e11.5_forebrain_expressedGene_withPLS.txt -wa -wb | awk '{FS=OFS="\t"}{if($7!=$15){print $4,$7,$15}}' | sort -u | awk 'BEGIN{FS=OFS="\t";id="";num=0}{if(NR==1){id=$1;num+=1}else{if(id==$1){num+=1}else{print id,num;id=$1;num=1}}}END{print id,num}' > e11.5_forebrain_center_regions_otherActiveGenesNum.txt

# 5. get expression in other 7 e11.5 tissues
awk 'BEGIN{FS=OFS="\t";print "loci","geneID"}{print $1":"$2"-"$3,$7}' e11.5_forebrain_center_regions.bed > e11.5_forebrain_center_regions_matrix.txt
while read line
do
    sample=`awk '{print $1}' <<< $line`
    expID=`awk '{print $2}' <<< $line`
    geneTsv=`awk '{print $3}' <<< $line`
    #
    awk -v sample="$sample" '{FS=OFS="\t"}{if(NR==FNR){split($1,c,".");a[c[1]]=1;b[c[1]]=$6}else{if(FNR==1){print $0,sample}else{split($2,d,".");if(a[d[1]]){print $0,b[d[1]]}}}}' /data/projects/encode/data/${expID}/${geneTsv}.tsv e11.5_forebrain_center_regions_matrix.txt > tmp.txt
    mv tmp.txt e11.5_forebrain_center_regions_matrix.txt
done < ${rna_filelist}

# 6. get UCSC link
## 1) ccRE tracks
# signal track
while read line
do
    dnase_exp=`awk '{print $1}' <<< $line`
    dnase_file=`awk '{print $2}' <<< $line`
    h3k4me3_exp=`awk '{print $3}' <<< $line`
    h3k4me3_file=`awk '{print $4}' <<< $line`
    h3k27ac_exp=`awk '{print $5}' <<< $line`
    h3k27ac_file=`awk '{print $6}' <<< $line`
    sample=`awk '{split($9,a,"_");print a[2]}' <<< $line`
    echo $sample
    #
    cp /data/projects/encode/data/${dnase_exp}/${dnase_file}.bigWig /data/public_html_users/fankaili/ccre/e11.5_${sample}_DNase.bigWig
    cp /data/projects/encode/data/${h3k4me3_exp}/${h3k4me3_file}.bigWig /data/public_html_users/fankaili/ccre/e11.5_${sample}_H3K4me3.bigWig
    cp /data/projects/encode/data/${h3k27ac_exp}/${h3k27ac_file}.bigWig /data/public_html_users/fankaili/ccre/e11.5_${sample}_H3K27ac.bigWig
done < e11.5_master_list.txt
### agnostic tracks
awk '{FS=OFS="\t"}{if($6=="PLS" || $6=="PLS,CTCF-bound"){print $1,$2,$3,$4,1,".",$2,$3,"255,0,0"}else if($6=="pELS" || $6=="pELS,CTCF-bound"){print $1,$2,$3,$4,1,".",$2,$3,"255,167,0"}else if($6=="dELS" || $6=="dELS,CTCF-bound"){print $1,$2,$3,$4,1,".",$2,$3,"255,205,0"}else if($6=="DNase-H3K4me3" || $6=="DNase-H3K4me3,CTCF-bound"){print $1,$2,$3,$4,1,".",$2,$3,"255,170,170"}else{print $1,$2,$3,$4,1,".",$2,$3,"0,176,240"}}' ${ccreDir}mm10-ccREs.bed | sort -k1,1 -k2,2n > /data/public_html_users/fankaili/ccre/mm10-ccREs_bed9.bed
bedToBigBed /data/public_html_users/fankaili/ccre/mm10-ccREs_bed9.bed /home/fankaili/genome/mm10.chrom.sizes.clean /data/public_html_users/fankaili/ccre/mm10-ccREs_bed9.bb
# e11.5 forebrain ccRE
awk '{FS=OFS="\t"}{if($2>1.64){print $1}}' ${signalDir}ENCSR014SFF-ENCFF009ALP.txt > /data/public_html_users/fankaili/ccre/tmp_active_DNase.txt
awk '{FS=OFS="\t"}{if($2>1.64){print $1}}' ${signalDir}ENCSR739DVM-ENCFF880ZZM.txt > /data/public_html_users/fankaili/ccre/tmp_active_H3K4me3.txt
awk '{FS=OFS="\t"}{if($2>1.64){print $1}}' ${signalDir}ENCSR275KPI-ENCFF827WXG.txt > /data/public_html_users/fankaili/ccre/tmp_active_H3K27ac.txt
cat /data/public_html_users/fankaili/ccre/tmp_active_H3K4me3.txt /data/public_html_users/fankaili/ccre/tmp_active_H3K27ac.txt | sort -u > /data/public_html_users/fankaili/ccre/tmp.active.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $1,"e11.5_forebrain_ccRE"}}}' /data/public_html_users/fankaili/ccre/tmp_active_DNase.txt /data/public_html_users/fankaili/ccre/tmp.active.txt > /data/public_html_users/fankaili/ccre/e11.5_forebrain_ccRE.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' /data/public_html_users/fankaili/ccre/e11.5_forebrain_ccRE.txt ${ccreDir}mm10-rDHSs.bed | sort -k1,1 -k2,2n > /data/public_html_users/fankaili/ccre/e11.5_forebrain_ccRE.bed
bedtools closest -a /data/public_html_users/fankaili/ccre/e11.5_forebrain_ccRE.bed -b /home/fankaili/genome/mm10_vM18_basic_TSS_filtered.bed -d | awk '{FS=OFS="\t"}{if($12>2000){print $1,$2,$3,$4,"distal"}else if($12<200){print $1,$2,$3,$4,"TSS"}else{print $1,$2,$3,$4,"proximal"}}' | sort -u > /data/public_html_users/fankaili/ccre/e11.5_forebrain_ccRE_distance.bed
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4] && $3=="TSS"){print $0,"PLS"}else{print $0,"ccRE"}}}' /data/public_html_users/fankaili/ccre/tmp_active_H3K4me3.txt /data/public_html_users/fankaili/ccre/e11.5_forebrain_ccRE_distance.bed > /data/public_html_users/fankaili/ccre/tmp1.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4] && $5=="proximal"){print $1,$2,$3,$4,$5,"pELS"}else if(a[$4] && $5=="TSS" && $6=="ccRE"){print $1,$2,$3,$4,$5,"pELS"}else if(a[$4] && $5=="distal"){print $1,$2,$3,$4,$5,"dELS"}else{print $0}}}' /data/public_html_users/fankaili/ccre/tmp_active_H3K27ac.txt /data/public_html_users/fankaili/ccre/tmp1.txt > /data/public_html_users/fankaili/ccre/tmp2.txt
awk '{FS=OFS="\t"}{if($6=="PLS"){print $1,$2,$3,$4,1,".",$2,$3,"255,0,0"}else if($6="pELS"){print $1,$2,$3,$4,1,".",$2,$3,"255,167,0"}else if($6=="dELS"){print $1,$2,$3,$4,1,".",$2,$3,"255,205,0"}else{print $1,$2,$3,$4,1,".",$2,$3,"255,170,170"}}' /data/public_html_users/fankaili/ccre/tmp2.txt | sort -k1,1 -k2,2n > /data/public_html_users/fankaili/ccre/e11.5_forebrain_ccREs_bed9.bed
bedToBigBed /data/public_html_users/fankaili/ccre/e11.5_forebrain_ccREs_bed9.bed /home/fankaili/genome/mm10.chrom.sizes.clean /data/public_html_users/fankaili/ccre/e11.5_forebrain_ccREs_bed9.bb
#
## 2) get links
# https://genome.ucsc.edu/s/Kaili/e11.5_regions
awk '{FS=OFS="\t"}{l=0.2*($3-$2);print $0,"https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A"$2-l"-"$3+l"&hgsid=722611183_YgcSzAt5gcOgBtdKaPnvzd4PYh4J"}' e11.5_forebrain_center_regions.bed > e11.5_forebrain_center_regions_withLink.bed


# 7. merge all, make tables
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$4]){print $1":"$2"-"$3,$7,b[$4]}else{print $1":"$2"-"$3,$7,0}}}' e11.5_forebrain_center_regions_otherActiveGenesNum.txt e11.5_forebrain_center_regions.bed > tmp1.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7}else{print $1,$2,a[$2],$3}}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt tmp1.bed > tmp1-1.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1":"$2"-"$3]=$15"\t"$13"\t"$14"\t"$11"\t"$12"\t"$9"\t"$16"\t"$10}else{print $0,a[$1]}}' e11.5_forebrain_center_regions_dELS.bed tmp1-1.bed > tmp2.bed
awk '{FS=OFS="\t"}{if(NR==FNR){if(NR>1){a[$1]=$0}}else{print $0,a[$1]}}' e11.5_forebrain_center_regions_matrix.txt tmp2.bed | cut -f 1-12,15-22 > tmp3.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1":"$2"-"$3]=$9}else{print $0,a[$1]}}' e11.5_forebrain_center_regions_withLink.bed tmp3.bed | sort -k5,5nr > e11.5_forebrain_center_regions_TABLE.txt
rm tmp*.bed



###########
# add another sheet for expression of other genes in the region.
mkdir other_expressed_gene
intersectBed -a e11.5_forebrain_center_regions.bed -b ./e11.5_expressed_genes/e11.5_forebrain_expressedGene_withPLS.txt -wa -wb | awk '{FS=OFS="\t"}{if($7!=$15){print $4,$15}}' | sort -u > ./other_expressed_gene/e11.5_forebrain_other_expressed_gene_list.txt
grep "ENSMUSG" /data/projects/encode/data/ENCSR160IIN/ENCFF465SNB.tsv | awk '{FS=OFS="\t"}{print $1,$6}' > tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){split($1,c,".");a[c[1]]=$2}else{split($2,d,".");print $0,a[d[1]]}}' tmp.txt ./other_expressed_gene/e11.5_forebrain_other_expressed_gene_list.txt > tmp2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7}else{print $1,$2,a[$2],$3}}' /home/fankaili/genome/mm10_vM18_basic_gene_filtered.txt tmp2.txt > tmp3.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$1":"$2"-"$3}else{print a[$1],$2,$3,$4}}' e11.5_forebrain_center_regions.bed tmp3.txt > tmp4.txt
awk 'BEGIN{FS=OFS="\t";id=""}{if(NR==1){id=$1;printf id"\t"$2"\t"$3"\t"$4}else if(id==$1){printf "\t"$2"\t"$3"\t"$4}else{id=$1;printf "\n"id"\t"$2"\t"$3"\t"$4}}END{printf "\n"}' tmp4.txt > tmp5.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$0;b[$1]=1}else{if(b[$1]){print $2,$3,$4,a[$1]}else{print $2,$3,$4,$1}}}' tmp5.txt e11.5_forebrain_center_regions_TABLE.txt > ./other_expressed_gene/e11.5_forebrain_other_expressed_gene_TABLE.txt
