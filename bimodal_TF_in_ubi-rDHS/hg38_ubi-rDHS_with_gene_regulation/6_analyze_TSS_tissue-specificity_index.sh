#！/bin/bash

# -- Kaili
# This script is for analyzing TSS tissue-specificity index.
# 1. RAMPAGE signal between biosamples
# 2. RAMPAGE signal between ubi-rDHS overlapped and not overlapped TSSs.
# 3. Scatter plot for TS index of all TSS in a gene and overlapped TSS.
# 4. Scatter plot for TS index of overlapped TSS in a gene and rest TSS.
# 5. maximum TSS TS index in gene

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/"

cd ${workDir}

# 0. sort TSS id
awk '{FS=OFS="\t"}{print $1,$2,$3,$6,$1"_"$2"_"$3"_"$6}' TSS.Filtered.bed | sort -u | sort -k1,1 -k2n | awk '{FS=OFS="\t"}{print $0,"TSS"NR}' \
> TSS.uniq.bed
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$6}else{if(a[$1"_"$2"_"$3"_"$6]){print $0,a[$1"_"$2"_"$3"_"$6]}}}' TSS.uniq.bed TSS.Filtered.bed \
> TSS.Filtered.uniqID.bed
#
sort -k8 -u TSS.Filtered.uniqID.bed > TSS.Filtered.uniq.bed
awk '{FS=OFS="\t"}{if($6=="+"){print $0}}' TSS.Filtered.uniq.bed > TSS.Filtered.uniq_plus.bed
awk '{FS=OFS="\t"}{if($6=="-"){print $0}}' TSS.Filtered.uniq.bed > TSS.Filtered.uniq_minus.bed

# get overlapped TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$11]=1}else{if(a[$4]){print $0}}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed TSS.Filtered.uniqID.bed \
> GRCh38_ubi-rOCR_overlapped_TSS.bed



# 1. RAMPAGE signal between biosamples
## function for getting RAMPAGE signal
getRampageSignal(){
    # setting path
    rampage_signal_dir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage/"
    tss_plus="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered.uniq_plus.bed"
    tss_minus="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered.uniq_minus.bed"
    output_dir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/rampage_tissue/"
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$4}else{print $8,a[$4]}}' ${rampage_signal_dir}$1.tab  ${tss_plus}> ${output_dir}$3_rampage.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$4}else{print $8,a[$4]}}' ${rampage_signal_dir}$2.tab  ${tss_minus} >> ${output_dir}$3_rampage.txt
}

## calculate signal
## choose 8 biosample here
mkdir rampage_tissue
getRampageSignal ENCFF612UBB ENCFF302QEN liver_32_year
getRampageSignal ENCFF039WHT ENCFF143FSY GM12878
getRampageSignal ENCFF438VNC ENCFF618GJP H7-hESC
getRampageSignal ENCFF391ZRW ENCFF730ESX stomach_40_week
getRampageSignal ENCFF224DNG ENCFF984PKS lung_24_week
getRampageSignal ENCFF417AXI ENCFF355LSJ spleen_53_year
getRampageSignal ENCFF043KVG ENCFF058KLC A172
getRampageSignal ENCFF783EAC ENCFF518WII K562

while read line
do
    id=`awk '{print $1}' <<< ${line}`
    plus=`awk '{print $2}' <<< ${line}`
    minus=`awk '{print $3}' <<< ${line}`
    biosample=`awk '{print $4}' <<< ${line}`
    #
    getRampageSignal ${plus} ${minus} ${biosample}
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt

## get ubi-rOCR overlapped gene with TSS file
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$2]=1}else{if(a[$7]){if(b[$8]){print $7,$8,"overlapped"}else{print $7,$8,"non_overlapped"}}}}' \
./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_TSS_uniqID.txt TSS.Filtered.uniqID.bed | sort -u > hg38_ubi-rOCR_overlapped_gene_all_TSS_uniqID.txt

# scatter plot
# Rscript make_rampage_signal_scatter_between_tissue.R



# 2. RAMPAGE signal between ubi-rDHS overlapped and not overlapped TSSs.
# Rscript make_rampage_signal_boxplot_between_ubi-oOCR_overlapped_nonoverlapped.R



# 3. Scatter plot for TS index of all TSS in a gene and overlapped TSS.
## get overlapped TSS TS index for each gene
awk '{FS=OFS="\t"}{if(NR==FNR){if($15=="0"){a[$11]=$14}}else{if(a[$4]){print a[$4],$8}}}' ./closest_gene/GRCh38_ubi-rOCR_closest_gene.bed \
TSS.Filtered.uniqID.bed | sort -u > ./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_TSS_uniqID.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{print $1,$2,b[$2]}}' hg38_tissue_TSS_exp_TSscore_uniqID.txt \
./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_TSS_uniqID.txt | sort -k1,1 > hg38_ubi-rOCR_overlapped_gene_TSS_TSindex.txt
# calculate average
awk '{FS=OFS="\t"}{if(NR==1){id=$1;sum=$3;n=1}else{if(id!=$1){print id,sum/n;id=$1;sum=$3;n=1}else{sum+=$3;n+=1}}}END{print id,sum/n}' \
hg38_ubi-rOCR_overlapped_gene_TSS_TSindex.txt > hg38_ubi-rOCR_overlapped_gene_TSS_TSindex_average.txt

## for all TSS of each gene
awk '{FS=OFS="\t"}{print $7,$8}' TSS.Filtered.uniqID.bed | sort -u > hg38_gene_TSS_uniqID.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{print $1,$2,b[$2]}}' hg38_tissue_TSS_exp_TSscore_uniqID.txt \
hg38_gene_TSS_uniqID.txt | sort -k1,1 > hg38_gene_TSS_TSindex.txt
# calculate average
awk '{FS=OFS="\t"}{if(NR==1){id=$1;sum=$3;n=1}else{if(id!=$1){print id,sum/n;id=$1;sum=$3;n=1}else{sum+=$3;n+=1}}}END{print id,sum/n}' \
hg38_gene_TSS_TSindex.txt > hg38_gene_TSS_TSindex_average.txt

## plot
# Rscript compare_TSindex.R


# 4. Scatter plot for TS index of overlapped TSS in a gene and rest TSS.
# for rest TSS of each gene
awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=1}else{if(a[$2]!=1){print $0}}}' ./closest_gene/GRCh38_ubi-rOCR_overlapped_gene_TSS_uniqID.txt \
hg38_gene_TSS_TSindex.txt > hg38_non_overlapped_gene_TSS_TSindex.txt
# calculate average
awk '{FS=OFS="\t"}{if(NR==1){id=$1;sum=$3;n=1}else{if(id!=$1){print id,sum/n;id=$1;sum=$3;n=1}else{sum+=$3;n+=1}}}END{print id,sum/n}' \
hg38_non_overlapped_gene_TSS_TSindex.txt > hg38_non_overlapped_gene_TSS_TSindex_average.txt

## plot
# Rscript compare_TSindex.R



# 5. maximum TSS TS index in gene
awk '{FS=OFS="\t"}{if(NR==1){id=$1;max=$3}else{if(id!=$1){print id,max;id=$1;max=$3}else{if(max<$3){max=$3}}}}END{print id,max}' \
hg38_gene_TSS_TSindex.txt > hg38_gene_TSS_TSindex_max.txt

## plot
# Rscript compare_TSindex.R
