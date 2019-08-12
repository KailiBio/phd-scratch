#!/bin/bash

# -- Kaili
# This script is for analyzing ubi-rOCRs & promoter shapes.
# 1. GC content: ubi-rOCRs vs. rOCRs with TSSs overlapping
# 2. CpG islands
# 3. get master peaks for all samples, get ubiquitous promoter peaks
# 4. enrichment of ubi-rOCRs in ubi-peaks
# 5. statistical analysis for each sample

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/promoter_shape/"
promoterShapeResult="/data/zusers/zhangx/projects/rampage/0_rampage_peak/peak_uniq/"

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


# 3. get master peaks for all samples, get ubiquitous promoter peaks
## 1) get merged peak file
# pool all reads in all samples together to call peaks
# resized to contain 95% of the reads
sort -k1,1 -k2,2n /data/zusers/zhangx/projects/rampage/0_rampage_peak/merge/merged_peak_resized.bed | \
awk '{FS=OFS="\t"}{print $1,$2,$3,"p_"NR,$5,$6}' > all_promoter_peak_list.txt
## 2) match promoter shape in each sample to merged peaks
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    echo ${expID}
    #
    awk '{FS=OFS="\t"}{if($9=="255,0,0"){print $1,$2,$3,"sp",$5,$6}else{print $1,$2,$3,"bp",$5,$6}}' ${promoterShapeResult}${expID}_rampage_peak_high.bed | sort -k1,1 -k2,2n > tmp.bed
    intersectBed -a all_promoter_peak_list.txt -b tmp.bed -wa -wb -s | awk 'BEGIN{FS=OFS="\t";name="";sp=0;bp=0}{if(NR==1){name=$4;if($10=="bp"){np+=1}else{sp+=1}}else{if($4==name){if($10=="bp"){np+=1}else{sp+=1}}else{print name,sp,np;name=$4;sp=0;bp=0;if($10=="bp"){bp+=1}else{sp+=1}}}}END{print name,sp,bp}' | awk 'BEGIN{FS=OFS="\t"}{if($2==0){print $0,"bp"}else if($3==0){print $0,"sp"}else{print $0,"mixed"}}' > ./promoter_shape_file/${expID}_promoter_type.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt
rm tmp.bed
## 3) count peaks occurence in all samples
awk '{FS=OFS="\t"}{print $4,0,0,0}' all_promoter_peak_list.txt > peak_type_count.txt
while read line
do
    expID=`awk '{print $1}' <<< ${line}`
    echo ${expID}
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){if($4=="sp"){a[$1]=1}else if($4=="np"){b[$1]=1}else{c[$1]=1}}else{if(a[$1]){print $1,$2+1,$3,$4}else if(b[$1]){print $1,$2,$3+1,$4}else if(c[$1]){print $1,$2,$3,$4+1}else{print $0}}}' ./promoter_shape_file/${expID}_promoter_type.txt peak_type_count.txt > tmp.txt
    mv tmp.txt peak_type_count.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_list.txt
awk '{FS=OFS="\t"}{print $0,$2+$3+$4}' peak_type_count.txt > tmp.txt
mv tmp.txt peak_type_count.txt
sed -i '1s/^/id\tsp\tnp\tmix\ttotal\n/' peak_type_count.txt
## 4) generate master peak lists, relate to rOCRs
# get master rOCRs file with label.
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0,"ubi-rOCR"}else{print $0,"non-ubi-rOCR"}}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed | sort -k1,1 -k2,2n > /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs_withLable.bed
# get peaks that overlap ubi-rOCRs
sort -k1,1 -k2,2n all_promoter_peak_list.txt > all_promoter_peak_list.bed
intersectBed -a all_promoter_peak_list.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs_withLable.bed -wa -wb -F 0.1 | sort -u > peak_overlap_rOCRs.txt
# get master list
awk '{FS=OFS="\t"}{if(NR>1){if($5>0){print $0}}}' peak_type_count.txt | sort -u | wc -l
awk '{FS=OFS="\t"}{if(NR==FNR){if($5>0 && NR>1){a[$1]=1;b[$1]=$0}}else{if(a[$4]){print $0,b[$4]}}}' peak_type_count.txt peak_overlap_rOCRs.txt | cut -f 1-6,10-11,13-16 | sort -u | sort -k1,1 -k2,2n > peak_master_list.txt

# Rscript make_peak_shape_figures.R

# peaks that not overlap rOCRs
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]!=1 && $5>0 && FNR>1){print $0}}}' peak_overlap_rOCRs.txt peak_type_count.txt | sort -u > peak_not-in-rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$0}else{if(a[$4]){print $0,b[$4]}}}' peak_not-in-rOCRs.txt all_promoter_peak_list.bed | cut -f 1-6,8-11 | sort -k1,1 -k2,2n > peak_not-in-rOCRs.bed
#
intersectBed -a peak_not-in-rOCRs.bed -b /home/fankaili/genome/hg38.blacklist.bed -u | head

# 4. enrichment of ubi-rOCRs in ubi-peaks
awk '{FS=OFS="\t"}{if($9>$10 && $9>$11){print $0,"sp"}else if($10>$9 && $10>$11){print $0,"bp"}else if($11>$9 && $11>$10){print $0,"mix"}else{print $0,"unclear"}}' peak_master_list.txt > peak_master_list_shape.txt
## 1) divide peaks into ubi-peaks & non-ubi-peaks
awk '{if($12>=145){print $0}}' peak_master_list_shape.txt > ubi_peaks.txt
awk '{if($12<145){print $0}}' peak_master_list_shape.txt > non-ubi_peaks.txt
## 2) ubi-rOCRs/non-ubi-rOCRs enrichment in ubi-peaks & cts-peaks
awk '{FS=OFS="\t"}{if($8=="ubi-rOCR"){print $4}}' ubi_peaks.txt | sort -u | wc -l
awk '{FS=OFS="\t"}{if($8=="ubi-rOCR"){print $4}}' non-ubi_peaks.txt | sort -u | wc -l
## 3) promoter shape of ubi-peaks
cut -f 4,13 ubi_peaks.txt | sort -u | cut -f 2 | sort | uniq -c
cut -f 4,13 non-ubi_peaks.txt | sort -u | cut -f 2 | sort | uniq -c
## 4) broad peaks & ubi-rOCRs_EDGEid
cut -f 4,13 peak_master_list_shape.txt | sort -u | cut -f 2 | sort | uniq -c
awk '{FS=OFS="\t"}{if($13=="bp"){print $0}}' ubi_peaks.txt > tmp.bp.txt
awk '{FS=OFS="\t"}{if($8=="ubi-rOCR"){print $4}}' tmp.bp.txt | sort -u | wc -l
awk '{FS=OFS="\t"}{if($13=="sp"){print $0}}' ubi_peaks.txt > tmp.sp.txt
awk '{FS=OFS="\t"}{if($8=="ubi-rOCR"){print $4}}' tmp.sp.txt | sort -u | wc -l
## 5) 3-way Venn: ubi-peaks, ubi-rOCRs, bp
awk '{FS=OFS="\t"}{if($13=="bp" && $8=="ubi-rOCR"){print $4}}' ubi_peaks.txt | sort -u | wc -l
awk '{FS=OFS="\t"}{if($13!="bp" && $8=="ubi-rOCR"){print $4}}' ubi_peaks.txt | sort -u | wc -l
awk '{FS=OFS="\t"}{if($13!="bp" && $8!="ubi-rOCR"){print $4}}' ubi_peaks.txt | sort -u | wc -l
awk '{FS=OFS="\t"}{if($13=="bp" && $8=="ubi-rOCR"){print $4}}' peak_master_list_shape.txt | sort -u | wc -l


# 5. statistical analysis for each sample
# Jun08
if [ -f fisher_matrix_each_sample.txt ];then rm fisher_matrix_each_sample.txt; fi
while read line
do
    dnase_expID=`awk '{print $1}' <<< ${line}`
    rampage_expID=`awk '{print $7}' <<< ${line}`
    sample=`awk '{print $3}' <<< ${line}`
    echo ${sample}
    #
    awk '{FS=OFS="\t"}{if($9=="255,0,0"){print $1,$2,$3,"np",$5,$6}else{print $1,$2,$3,"bp",$5,$6}}' ${promoterShapeResult}${rampage_expID}_rampage_peak_high.bed | awk '{FS=OFS="\t"}{print $0,"peak_"NR}' > tmp.bed
    # ubi-rOCRs
    intersectBed -a tmp.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb | cut -f 4,7 | sort -u > tmp.txt
    ubi_bp=`grep "bp" tmp.txt | wc -l`
    ubi_np=`grep "np" tmp.txt | wc -l`
    # non-ubi rOCRs
    intersectBed -a tmp.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/${dnase_expID}_OCR.bed -wa -wb | cut -f 4,7 | sort -u > tmp2.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=1}else{if(a[$2]!=1){print $0}}}' tmp.txt tmp2.txt > tmp3.txt
    non_ubi_bp=`grep "bp" tmp3.txt | wc -l`
    non_ubi_np=`grep "np" tmp3.txt | wc -l`
    #
    echo -e ${sample}"\t"${ubi_bp}"\t"${ubi_np}"\t"${non_ubi_bp}"\t"${non_ubi_np} >> fisher_matrix_each_sample.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_DNase_matched_list.txt

# Rscript make_peak_shape_figures.R

### rOCRs overlapping TSSs only
if [ -f fisher_matrix_each_sample_2.txt ];then rm fisher_matrix_each_sample_2.txt; fi
while read line
do
    dnase_expID=`awk '{print $1}' <<< ${line}`
    rampage_expID=`awk '{print $7}' <<< ${line}`
    sample=`awk '{print $3}' <<< ${line}`
    echo ${sample}
    #
    awk '{FS=OFS="\t"}{if($9=="255,0,0"){print $1,$2,$3,"np",$5,$6}else{print $1,$2,$3,"bp",$5,$6}}' ${promoterShapeResult}${rampage_expID}_rampage_peak_high.bed | awk '{FS=OFS="\t"}{print $0,"peak_"NR}' > tmp.bed
    # ubi-rOCRs
    intersectBed -a tmp.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed -wa -wb | cut -f 4,7 | sort -u > tmp.txt
    ubi_bp=`grep "bp" tmp.txt | wc -l`
    ubi_np=`grep "np" tmp.txt | wc -l`
    # non-ubi rOCRs
    intersectBed -a /data/zusers/fankaili/ccre/hg38_ubi-rDHS/non_ubi_active_OCR/${dnase_expID}_OCR.bed -b /home/fankaili/genome/mm10_vM18_basic_TSS_filtered.bed -wa -u | sort -k1,1 -k2,2n > tmp2.bed
    intersectBed -a tmp.bed -b tmp2.bed -wa -wb | cut -f 4,7 | sort -u > tmp2.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=1}else{if(a[$2]!=1){print $0}}}' tmp.txt tmp2.txt > tmp3.txt
    non_ubi_bp=`grep "bp" tmp3.txt | wc -l`
    non_ubi_np=`grep "np" tmp3.txt | wc -l`
    #
    echo -e ${sample}"\t"${ubi_bp}"\t"${ubi_np}"\t"${non_ubi_bp}"\t"${non_ubi_np} >> fisher_matrix_each_sample_2.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_DNase_matched_list.txt

### all rOCRs
if [ -f fisher_matrix_each_sample_3.txt ];then rm fisher_matrix_each_sample_3.txt; fi
while read line
do
    dnase_expID=`awk '{print $1}' <<< ${line}`
    rampage_expID=`awk '{print $7}' <<< ${line}`
    sample=`awk '{print $3}' <<< ${line}`
    echo ${sample}
    #
    awk '{FS=OFS="\t"}{if($9=="255,0,0"){print $1,$2,$3,"np",$5,$6}else{print $1,$2,$3,"bp",$5,$6}}' ${promoterShapeResult}${rampage_expID}_rampage_peak_high.bed | awk '{FS=OFS="\t"}{print $0,"peak_"NR}' > tmp.bed
    # ubi-rOCRs
    intersectBed -a tmp.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/ubi-rOCRs_overlap_hg38_v28_basic_TSS.bed -wa -wb | cut -f 4,7 | sort -u > tmp.txt
    ubi_bp=`grep "bp" tmp.txt | wc -l`
    ubi_np=`grep "np" tmp.txt | wc -l`
    # non-ubi rOCRs
    intersectBed -a tmp.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -wa -wb | cut -f 4,7 | sort -u > tmp2.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$2]=1}else{if(a[$2]!=1){print $0}}}' tmp.txt tmp2.txt > tmp3.txt
    non_ubi_bp=`grep "bp" tmp3.txt | wc -l`
    non_ubi_np=`grep "np" tmp3.txt | wc -l`
    #
    echo -e ${sample}"\t"${ubi_bp}"\t"${ubi_np}"\t"${non_ubi_bp}"\t"${non_ubi_np} >> fisher_matrix_each_sample_3.txt
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_DNase_matched_list.txt








while read line
do
    dnase_expID=`awk '{print $1}' <<< ${line}`
    rampage_expID=`awk '{print $7}' <<< ${line}`
    sample=`awk '{print $3}' <<< ${line}`
    echo ${sample}
    #
    awk '{FS=OFS="\t"}{if($9=="255,0,0"){print $1,$2,$3,"np",$5,$6}else{print $1,$2,$3,"bp",$5,$6}}' ${promoterShapeResult}${rampage_expID}_rampage_peak_high.bed | awk '{FS=OFS="\t"}{print $0,"peak_"NR}' > tmp.bed
    wc -l ${promoterShapeResult}${rampage_expID}_rampage_peak_high.bed
    grep "bp" tmp.bed | wc -l
    grep "np" tmp.bed | wc -l
done < /data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_RAMPAGE_DNase_matched_list.txt


## Jul21
new_peak_shape_enrichment.sh



define_RAMPAGE_peak_cutoff.sh



########
# Aug04
# make track hub
cd /data/public_html_users/fankaili/GRCh38_rOCRs/
awk '{FS=OFS="\t"}{print $0,"1",".",$2,$3,"140,140,140"}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed | sort -k1,1 -k2,2n > GRCh38-rOCRs_bed9.bed0
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $1,$2,$3,$4,$5,$6,$7,$8,"6,218,14"}else{print $0}}}' GRCh38_ubi-rOCRs.bed GRCh38-rOCRs_bed9.bed0 | sort -k1,1 -k2,2n > GRCh38-rOCRs_bed9.bed
bedToBigBed GRCh38-rOCRs_bed9.bed /home/fankaili/genome/hg38.chrom.sizes.clean GRCh38-rOCRs_bed9.bb
vim trackDb_ubi-OCR_bed.txt
