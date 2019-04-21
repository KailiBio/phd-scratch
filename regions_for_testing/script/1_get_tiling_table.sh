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
rna_filelist="/data/zusers/fankaili/regions/mouse_e11.5_RNA-seq_list.txt"

cd ${workDir}

# 0. data preparation
## get master file in e11.5
grep "C57BL/6" ${ccreDir}Cell-Type-Specific/Master-Cell-List.txt | grep "embryo_11.5_days" > e11.5_master_list.txt
## mouse GENCODE M4 GTF
bash /data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/process_GENCODE_gtf.sh /home/fankaili/genome/gencode.vM4.annotation.gtf mm10 vM4 comprehensive


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
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$3}else{if(a[$7]){print $0,b[$7]}}}' ./e11.5_expressed_genes/e11.5_${sample}_expressed_gene.txt /home/fankaili/genome/mm10_vM4_comprehensive_TSS_filtered.bed | sort -k1,1 -k2,2n | awk '{if(NR>1){print $0}}' > ./e11.5_expressed_genes/e11.5_${sample}_expressed_gene_TSS.txt
done < ${rna_filelist}
## 2) get PLS in all tissues
mkdir e11.5_PLS
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
    awk '{if(NR==FNR){if($2>1.64){a[$1]=1}}else{if($2>1.64 && a[$1]){print $1}}}' ${signalDir}${dnase_exp}-${dnase_file}.txt ${signalDir}${h3k4me3_exp}-${h3k4me3_file}.txt | sort -u > ./e11.5_PLS/e11.5_${sample}_PLS.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' ./e11.5_PLS/e11.5_${sample}_PLS.txt ${ccreDir}mm10-rDHSs.bed | sort -k1,1 -k2,2n > ./e11.5_PLS/e11.5_${sample}_PLS.bed
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
    awk '{if(NR==FNR){if($2>1.64){a[$1]=1}}else{if($2>1.64 && a[$1]){print $1}}}' ${signalDir}${dnase_exp}-${dnase_file}.txt ${signalDir}${h3k27ac_exp}-${h3k27ac_file}.txt | sort -u > ./e11.5_ELS/e11.5_${sample}_ELS.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(a[$4]){print $0}}}' ./e11.5_ELS/e11.5_${sample}_ELS.txt ${ccreDir}mm10-rDHSs.bed | sort -k1,1 -k2,2n > ./e11.5_ELS/e11.5_${sample}_ELS.bed
    #
    awk -v sample="$sample" '{if(NR==1){print $0,sample}}' e11.5_forebrain_center_regions_dELS.bed > tmp.bed
    intersectBed -a e11.5_forebrain_center_regions_dELS.bed -b ./e11.5_ELS/e11.5_${sample}_ELS.bed -c >> tmp.bed
    mv tmp.bed e11.5_forebrain_center_regions_dELS.bed
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
    awk -v sample="$sample" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$6}else{if(FNR==1){print $0,sample}else{if(a[$2]){print $0,b[$2]}}}}' /data/projects/encode/data/${expID}/${geneTsv}.tsv e11.5_forebrain_center_regions_matrix.txt > tmp.txt
    mv tmp.txt e11.5_forebrain_center_regions_matrix.txt
done < ${rna_filelist}

# 6. get UCSC session
awk '{FS=OFS="\t"}{if($6=="dELS" || $6=="dELS,CTCF-bound"){print $1,$2,$3,$4,1,".",$2,$3,"255,121,3"}else if($6=="pELS" || $6=="pELS,CTCF-bound"){print $1,$2,$3,$4,1,".",$2,$3,"255,205,0"}else if($6=="PLS" || $6=="PLS,CTCF-bound"){print $1,$2,$3,$4,1,".",$2,$3,"255,0,0"}else{print $1,$2,$3,$4,1,".",$2,$3,"140,140,140"}}' ${ccreDir}mm10-ccREs.bed | sort -k1,1 -k2,2n > /data/public_html_users/fankaili/ccre/mm10-ccREs_bed9.bed
bedToBigBed /data/public_html_users/fankaili/ccre/mm10-ccREs_bed9.bed /home/fankaili/genome/mm10.chrom.sizes.clean /data/public_html_users/fankaili/ccre/mm10-ccREs_bed9.bb
#
# https://genome.ucsc.edu/s/Kaili/e11.5_regions
awk '{FS=OFS="\t"}{l=0.2*($3-$2);print $0,"https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="$1"%3A"$2-l"-"$3+l"&hgsid=722068975_Z6MlV9ysQMRMkTiBT6lgo6uLr5GJ"}' e11.5_forebrain_center_regions.bed > e11.5_forebrain_center_regions_withLink.bed

# 7. merge all, make tables
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$4]){print $1":"$2"-"$3,$7,b[$4]}else{print $1":"$2"-"$3,$7,0}}}' e11.5_forebrain_center_regions_otherActiveGenesNum.txt e11.5_forebrain_center_regions.bed > tmp1.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1":"$2"-"$3]=$15"\t"$13"\t"$14"\t"$11"\t"$12"\t"$9"\t"$16"\t"$10}else{print $0,a[$1]}}' e11.5_forebrain_center_regions_dELS.bed tmp1.bed > tmp2.bed
awk '{FS=OFS="\t"}{if(NR==FNR){if(NR>1){a[$1]=$0}}else{print $0,a[$1]}}' e11.5_forebrain_center_regions_matrix.txt tmp2.bed | cut -f 1-11,14-21 > tmp3.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1":"$2"-"$3]=$9}else{print $0,a[$1]}}' e11.5_forebrain_center_regions_withLink.bed tmp3.bed | sort -k4,4nr > e11.5_forebrain_center_regions_TABLE.txt
# get gene symbol
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7}else{print a[$2]}}' /home/fankaili/genome/mm10_vM4_comprehensive_gene_filtered.txt e11.5_forebrain_center_regions_TABLE.txt > e11.5_forebrain_center_regions_geneSymbol.txt
#
rm tmp*.bed
