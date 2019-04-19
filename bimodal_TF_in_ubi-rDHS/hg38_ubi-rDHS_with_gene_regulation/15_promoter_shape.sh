#!/bin/bash

# -- Kaili
# This script is for analyzing ubi-rOCRs & promoter shapes.
# 1. GC content: ubi-rOCRs vs. rOCRs with TSSs overlapping
# 2. CpG islands
# 3. pre-process

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/"
promoterShapeResult="/data/zusers/zhangx/projects/rampage/0_rampage_peak/entropy/promoter_type/"

cd ${workDir}

# Rscript basic_figures_for_promoter_shape.R

# 1. GC content: ubi-rOCRs vs. rOCRs with TSSs overlapping
# get TSS overlapped rOCRs
intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/TSS.Filtered.uniq.bed -wa |sort -u | sort -k1,1 -k2,2n > TSS_overlapped_rOCRs.bed
# get fasta
bedtools getfasta -fi /home/fankaili/genome/hg38.fa -bed TSS_overlapped_rOCRs.bed > TSS_overlapped_rOCRs.fa
python ${dailyCodeDir}count_GC_content.py /data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/TSS_overlapped_rOCRs.fa /data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/TSS_overlapped_rOCRs_GCcontent.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{id=$1":"$2"-"$3;if(b[id]){print $0,a[id]}}}' TSS_overlapped_rOCRs_GCcontent.txt TSS_overlapped_rOCRs.bed > TSS_overlapped_rOCRs_GCcontent_ID.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0,"ubi-rOCRs"}else{print $0,"remaining_rOCRs"}}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed TSS_overlapped_rOCRs_GCcontent_ID.txt > TSS_overlapped_rOCRs_GCcontent_ID_annotated.txt
### TSS in basic
cd /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/
intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -b hg38_v28_basic_TSS_filtered_uniq_labeled.bed -wa -u > rOCRs_overlap_hg38_v28_basic_TSS.bed
bedtools getfasta -fi /home/fankaili/genome/hg38.fa -bed rOCRs_overlap_hg38_v28_basic_TSS.bed > rOCRs_overlap_hg38_v28_basic_TSS.fa
python ${dailyCodeDir}count_GC_content.py /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/rOCRs_overlap_hg38_v28_basic_TSS.fa /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/rOCRs_overlap_hg38_v28_basic_TSS_GCcontent.txt
# label
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{id=$1":"$2"-"$3;if(b[id]){print $1,$2,$3,$4,"rOCRs_overlap_TSS",a[id]}}}' rOCRs_overlap_hg38_v28_basic_TSS_GCcontent.txt /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed > tmp.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $1,$2,$3,$4,"ubi-rOCRs_overlap_TSS",$6}else{print $0}}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed tmp.txt > rOCRs_overlap_hg38_v28_basic_TSS_GCcontent_ID.txt
rm tmp.txt


# 2. CpG islands
cd /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/
intersectBed -a rOCRs_overlap_hg38_v28_basic_TSS.bed -b /home/fankaili/genome/hg38_cpgIsland_sorted.txt -wa -u | wc -l
intersectBed -a ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed -b /home/fankaili/genome/hg38_cpgIsland_sorted.txt -wa -u | wc -l




# 3. pre-process
## 1) get master_list of peaks
sort -k1,1 -k2,2n /data/zusers/zhangx/projects/rampage/0_rampage_peak/merge/merged_peak.txt | \
awk '{FS=OFS="\t"}{print $1,$2,$3,"p_"NR,$5,$6,$7}' > all_promoter_peak_list.txt
## 2) get promoter shapes file
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    sample=`awk '{print $4}' <<< ${line}`
    echo ${expID}
    #
    cut -f 1-4,13,6,17 ${promoterShapeResult}${expID}_rampage_entropy_promoter_type.txt | sort -k1,1 -k2,2n > tmp.bed
    intersectBed -a tmp.bed -b all_promoter_peak_list.txt -f 0.5 -F 0.5 -e -wa -wb | \
    awk '{FS=OFS="\t"}{if($5==$13){print $1,$2,$3,$11,$5,$6,$7}}' \
    | sort -k1,1 -k2,2n > ./promoter_shape_file/${expID}_promoter_type0.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt
rm tmp.bed
## 3) count peaks occurence in all samples
if [ -f tmp.txt ]; then rm tmp.txt; fi
#
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    #
    cut -f 4,7 ./promoter_shape_file/${expID}_promoter_type0.txt | sort -u >> tmp.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt
#
cut -f 1 tmp.txt | sort | uniq -c | awk '{OFS="\t"}{print $2,$1}' > tmp2.txt
grep "np" tmp.txt | sort | uniq -c | awk '{OFS="\t"}{print $2,$3,$1}' > tmp3.txt
grep "bp" tmp.txt | sort | uniq -c | awk '{OFS="\t"}{print $2,$3,$1}' > tmp4.txt
grep "wp" tmp.txt | sort | uniq -c | awk '{OFS="\t"}{print $2,$3,$1}' > tmp5.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$3}else{if(a[$1]){print $0,a[$1]}else{print $0,0}}}' tmp3.txt tmp2.txt \
> tmp6.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$3}else{if(a[$1]){print $0,a[$1]}else{print $0,0}}}' tmp4.txt tmp6.txt \
> tmp7.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$3}else{if(a[$1]){print $0,a[$1]}else{print $0,0}}}' tmp5.txt tmp7.txt \
> tmp8.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$0;b[$1]=1}else{if(b[$4]){print $0,a[$4]}}}' tmp8.txt \
all_promoter_peak_list.txt | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{print $1,$2,$3,"peak_"NR,$5,$6,$7,$8,$9,$10,$11,$12}' \
> promoter_peak_master_list.txt
sed -i '1s/^/chr\tstart\tend\tpeak\tn\tstrand\tsignal\tname\ttotal\tnp\tbp\twp\n/' promoter_peak_master_list.txt
## 4) change the peak name of promoter shape file
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    #
    awk '{FS=OFS="\t"}{}' promoter_peak_master_list.txt \
    ./promoter_shape_file/${expID}_promoter_type0.txt > ./promoter_shape_file/${expID}_promoter_type.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt



# 3. peaks with rOCRs and ubi-rOCRs
intersectBed -a promoter_peak_master_list.txt -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -f 0.5 -F 0.5 \
-e -wa -u > tmp.promoter_peaks_overlappng_rOCRs.txt
intersectBed -a promoter_peak_master_list.txt -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed \
-f 0.5 -F 0.5 -e -wa -u > tmp.promoter_peaks_overlappng_ubi-rOCRs.txt









# promoter shapes











cut -f 1-3,17 ${promoterShapeResult}ENCSR855LXJ_rampage_entropy_promoter_type.txt > ss.txt
intersectBed -a TSS_overlapped_rOCRs_GCcontent_ID_annotated.txt -b ss.txt -wa -wb > ss2.txt
cut -f 4,6,10 ss2.txt | sort -u > ss3.txt
