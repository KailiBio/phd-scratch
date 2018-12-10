#!/bin/bash

# -- Kaili
# This script is for analyzing 9,448 Hi-C loops.
# 0. pre-processing
# 1. loops overlapped with ubi-rOCRs
# 2. loops overlapped with rOCRs
# 3. loops overlapped with CTCF peak
# 4. loop anchors overlapped with CTCF motif
# 5. loops overlapped with GM12878 specific active-rOCRs
# 6. CTCF singal between ubi-rOCRs overlapped peaks loci and remaining active rOCRs

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/loop/"

# 0. pre-processing
cd /data/zusers/fankaili/ccre/tf/loop/
## 1) get loops
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary%2Breplicate_HiCCUPS_looplist.txt.gz
#
awk '{FS=OFS="\t"}{if(NR>1){print "chr"$1,$2,$3,"loop_"(NR-1)"_a","loop_"(NR-1)}}' \
GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt > GM12878_loop_a_hg19.txt
awk '{FS=OFS="\t"}{if(NR>1){print "chr"$4,$5,$6,"loop_"(NR-1)"_b","loop_"(NR-1)}}' \
GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt > GM12878_loop_b_hg19.txt

## 2) liftover to hg38
liftOver GM12878_loop_a_hg19.txt /home/fankaili/genome/hg19ToHg38.over.chain \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/loop/GM12878_loop_a_hg38.txt \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/loop/GM12878_loop_a_hg38_unMapped.txt
#
liftOver GM12878_loop_b_hg19.txt /home/fankaili/genome/hg19ToHg38.over.chain \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/loop/GM12878_loop_b_hg38.txt \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/loop/GM12878_loop_b_hg38_unMapped.txt
# merge
cd ${workDir}
awk '{FS=OFS="\t"}{if(NR==FNR){a[$5]=$0;b[$5]=1}else{if(b[$5]){print a[$5],$0}}}' GM12878_loop_a_hg38.txt GM12878_loop_b_hg38.txt > \
hg38_GM12878_loop.txt
### get 9441 hic loops
cut -f 1-5 hg38_GM12878_loop.txt > hg38_GM12878_loop_a.bed
cut -f 6-10 hg38_GM12878_loop.txt > hg38_GM12878_loop_b.bed
### get peak loci
cat hg38_GM12878_loop_a.bed hg38_GM12878_loop_b.bed | cut -f 1-3 | sort -u | sort -k1,1 -k2,2n > GRCh38_loop_all_peak_loci.txt
bedtools merge -i GRCh38_loop_all_peak_loci.txt | awk '{FS=OFS="\t"}{print $0,"peak_loci_"NR}' > GRCh38_loop_peak_loci.txt
#
intersectBed -a hg38_GM12878_loop_a.bed -b GRCh38_loop_peak_loci.txt -wa -wb > GRCh38_loop_peak_loci_ID.txt
intersectBed -a hg38_GM12878_loop_b.bed -b GRCh38_loop_peak_loci.txt -wa -wb >> GRCh38_loop_peak_loci_ID.txt



# 1. loops overlapped with ubi-rOCRs
intersectBed -a hg38_GM12878_loop_a.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb > GRCh38_loop_a_overlapped_ubi-rOCRs.txt
intersectBed -a hg38_GM12878_loop_b.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa -wb > GRCh38_loop_b_overlapped_ubi-rOCRs.txt
#
cat GRCh38_loop_a_overlapped_ubi-rOCRs.txt GRCh38_loop_b_overlapped_ubi-rOCRs.txt > GRCh38_loop_overlapped_ubi-rOCRs.txt
cut -f 5 GRCh38_loop_overlapped_ubi-rOCRs.txt | sort -u | wc -l
# 1857 loops with at least one side overlapped with ubi-rOCRs
cut -f 10,11 GRCh38_loop_overlapped_ubi-rOCRs.txt | sort -u | wc -l
# 2165 ubi-rOCRs overlapped with loop anchors
cut -f 10,11 GRCh38_loop_overlapped_ubi-rOCRs.txt | sort -u | cut -f 2 | sort | uniq -c
# 1781 PLS, 376 ELS, 8 CTCF-only
cut -f 4,5 GRCh38_loop_overlapped_ubi-rOCRs.txt | sort -u | cut -f 2 | sort | uniq -d | wc -l
#

# 2. loops overlapped with rOCRs
intersectBed -a hg38_GM12878_loop_a.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -wa -wb > GRCh38_loop_a_overlapped_rOCRs.txt
intersectBed -a hg38_GM12878_loop_b.bed -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed -wa -wb > GRCh38_loop_b_overlapped_rOCRs.txt
#
cat GRCh38_loop_a_overlapped_rOCRs.txt GRCh38_loop_b_overlapped_rOCRs.txt > GRCh38_loop_overlapped_rOCRs.txt
cut -f 5 GRCh38_loop_overlapped_rOCRs.txt | sort -u | wc -l
# 9341 loops with at least one side overlapped with rOCRs_EDGEid
cut -f 9 GRCh38_loop_overlapped_rOCRs.txt | sort -u | wc -l
# 99,167 rOCRs overlapped with loop anchors
cut -f 4,5 GRCh38_loop_overlapped_rOCRs.txt | sort -u | cut -f 2 | sort | uniq -d | wc -l
# 9262 loops both end overlapped with rOCRs


# 3. loops overlapped with CTCF peak
## get CTCF peak file
cp /data/projects/encode/data/ENCSR000DRZ/ENCFF239OQV.bed.gz ./
gzip -d ENCFF239OQV.bed.gz
mv ENCFF239OQV.bed GM12878_CTCF_peaks_ENCSR000DRZ.bed
## intersect with peak
intersectBed -a hg38_GM12878_loop_a.bed -b GM12878_CTCF_peaks_ENCSR000DRZ.bed -wa -wb > \
GRCh38_loop_a_overlapped_CTCF_peaks.txt
intersectBed -a hg38_GM12878_loop_b.bed -b GM12878_CTCF_peaks_ENCSR000DRZ.bed -wa -wb > \
GRCh38_loop_b_overlapped_CTCF_peaks.txt
#
cat GRCh38_loop_a_overlapped_CTCF_peaks.txt GRCh38_loop_b_overlapped_CTCF_peaks.txt > \
GRCh38_loop_overlapped_CTCF_peaks.txt
cut -f 5 GRCh38_loop_overlapped_CTCF_peaks.txt | sort -u | wc -l
#  9188 loops have at least one anchor overlapped with CTCF peaks
cut -f 4,5 GRCh38_loop_overlapped_CTCF_peaks.txt | sort -u | cut -f 2 | sort | uniq -d | wc -l
# 7364 with both ends
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{print a[$4]}}' GRCh38_loop_peak_loci_ID.txt \
GRCh38_loop_overlapped_CTCF_peaks.txt | sort -u | wc -l
# 11327 peak loci


# 4. loop anchors overlapped with CTCF motif
## get CTCF motif file
cp /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/mm10_ctcf_motif_file.motif ./CTCF_motif_file.motif
# get CTCF motif loci in GRCh38
scanMotifGenomeWide.pl CTCF_motif_file.motif /home/fankaili/genome/hg38.fa -bed > GRCh38_CTCF_motif_region.bed
# intersect
intersectBed -a hg38_GM12878_loop_a.bed -b GRCh38_CTCF_motif_region.bed -wa -wb \
> GRCh38_loop_a_overlapped_CTCF_motif.txt
intersectBed -a hg38_GM12878_loop_b.bed -b GRCh38_CTCF_motif_region.bed -wa -wb \
> GRCh38_loop_b_overlapped_CTCF_motif.txt
#
cat GRCh38_loop_a_overlapped_CTCF_motif.txt GRCh38_loop_b_overlapped_CTCF_motif.txt > \
GRCh38_loop_overlapped_CTCF_motif.txt
cut -f 5 GRCh38_loop_overlapped_CTCF_motif.txt | sort -u | wc -l
#  7636 loops have at least one anchor overlapped with CTCF motif
cut -f 4,5 GRCh38_loop_overlapped_CTCF_motif.txt | sort -u | cut -f 2 | sort | uniq -d | wc -l
# 3555 with both ends
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{print a[$4]}}' GRCh38_loop_peak_loci_ID.txt \
GRCh38_loop_overlapped_CTCF_motif.txt | sort -u | wc -l
# 7742 peak loci

scanMotifGenomeWide.pl CTCF_motif_file_jasper.motif /home/fankaili/genome/hg38.fa -bed > GRCh38_CTCF_motif_region_2.bed
#
intersectBed -a hg38_GM12878_loop_a.bed -b GRCh38_CTCF_motif_region_2.bed \
> GRCh38_loop_a_overlapped_CTCF_motif_2.txt
intersectBed -a hg38_GM12878_loop_b.bed -b GRCh38_CTCF_motif_region_2.bed \
> GRCh38_loop_b_overlapped_CTCF_motif_2.txt
#
cat GRCh38_loop_a_overlapped_CTCF_motif_2.txt GRCh38_loop_b_overlapped_CTCF_motif_2.txt > \
GRCh38_loop_overlapped_CTCF_motif_2.txt
cut -f 5 GRCh38_loop_overlapped_CTCF_motif_2.txt | sort -u | wc -l
#  8316 loops have at least one anchor overlapped with CTCF motif
cut -f 4,5 GRCh38_loop_overlapped_CTCF_motif_2.txt | sort -u | cut -f 2 | sort | uniq -d | wc -l
# 4576 with both ends
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{print a[$4]}}' GRCh38_loop_peak_loci_ID.txt \
GRCh38_loop_overlapped_CTCF_motif_2.txt | sort -u | wc -l
# 8947 peak loci


# 5. loops overlapped with GM12878 specific active-rOCRs
## 1) get GM12878 specific active rOCRs
bigWigAverageOverBed /data/projects/encode/data/ENCSR000EMT/ENCFF915DFR.bigWig \
/data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed GRCh38_GM12878_rOCRs_DNase_signal.txt
## calculate z-score
python ${scriptDir}zscore-normalization.py GRCh38_GM12878_rOCRs_DNase_signal.txt 5 > \
GRCh38_GM12878_rOCRs_DNase_z-score.txt
sed -i 's/ \t/\t/g' GRCh38_GM12878_rOCRs_DNase_z-score.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){if($2>1.64){a[$1]=1}}else{if(a[$4]){print $0}}}' \
GRCh38_GM12878_rOCRs_DNase_z-score.txt /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38-rOCRs.bed > \
GRCh38_GM12878_active_rOCRs.bed

## 2)
intersectBed -a hg38_GM12878_loop_a.bed -b GRCh38_GM12878_active_rOCRs.bed -wa -wb > GRCh38_loop_a_overlapped_GM12878_active_rOCRs.txt
intersectBed -a hg38_GM12878_loop_b.bed -b GRCh38_GM12878_active_rOCRs.bed -wa -wb > GRCh38_loop_b_overlapped_GM12878_active_rOCRs.txt
#
cat GRCh38_loop_a_overlapped_GM12878_active_rOCRs.txt GRCh38_loop_b_overlapped_GM12878_active_rOCRs.txt \
> GRCh38_loop_overlapped_GM12878_active_rOCRs.txt
#
cut -f 5 GRCh38_loop_overlapped_GM12878_active_rOCRs.txt | sort -u | wc -l
# 7175 loops with at least one side overlapped with GM12878 active-rOCRs
cut -f 9 GRCh38_loop_overlapped_GM12878_active_rOCRs.txt | sort -u | wc -l
# 18,010 active-rOCRs overlapped with loop anchors
cut -f 4,5 GRCh38_loop_overlapped_GM12878_active_rOCRs.txt | sort -u | cut -f 2 | sort | uniq -d | wc -l
# 3560 loops both end overlapped with active-rOCRs


# 6. CTCF singal between ubi-rOCRs overlapped peaks loci and remaining active rOCRs
## 1) get peak loci CTCF z-score
bigWigAverageOverBed /data/projects/encode/data/ENCSR000DRZ/ENCFF852CRG.bigWig GRCh38_loop_peak_loci.txt \
GRCh38_peakLoci_GM12878_CTCF_signal.txt
## calculate z-score
python ${scriptDir}zscore-normalization.py GRCh38_peakLoci_GM12878_CTCF_signal.txt 5 > \
GRCh38_peakLoci_GM12878_CTCF_z-score.txt
sed -i 's/ \t/\t/g' GRCh38_peakLoci_GM12878_CTCF_z-score.txt

## 2) mark peak loci
awk '{FS=OFS="\t"}{print $0,"remaining"}' GRCh38_loop_peak_loci.txt > GRCh38_loop_peak_loci_marked.txt
## active rOCRs
intersectBed -a GRCh38_loop_peak_loci.txt -b GRCh38_GM12878_active_rOCRs.bed -wa > \
GRCh38_peakLoci_overlapped_active_rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $1,$2,$3,$4,"active-rOCR"}else{print $0}}}' \
GRCh38_peakLoci_overlapped_active_rOCRs.txt GRCh38_loop_peak_loci_marked.txt > temp.txt
## ubi-rOCRs
intersectBed -a GRCh38_loop_peak_loci.txt -b /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed -wa > \
GRCh38_peakLoci_overlapped_ubi-rOCRs.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $1,$2,$3,$4,"ubi-rOCR"}else{print $0}}}' \
GRCh38_peakLoci_overlapped_ubi-rOCRs.txt temp.txt > GRCh38_loop_peak_loci_marked.txt

## 3) get peak loci mark with z-zscore, make figures
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$4]}}' GRCh38_peakLoci_GM12878_CTCF_z-score.txt \
GRCh38_loop_peak_loci_marked.txt >  GRCh38_peakLoci_GM12878_CTCF_z-score_marked.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $0,a[$4]}}' GRCh38_peakLoci_GM12878_CTCF_signal.txt \
GRCh38_loop_peak_loci_marked.txt >  GRCh38_peakLoci_GM12878_CTCF_signal_marked.txt
# Rscript make_histgram_CTCF_hic_peakLoci.R
