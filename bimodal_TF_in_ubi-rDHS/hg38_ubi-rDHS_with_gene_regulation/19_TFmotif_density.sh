#!/bin/bash

# -- Kaili
# This script is for comparing TF motif density between ubi-rOCRs & rOCRs.
# 1. get sequence
# 2. get motif list
# 3. run FIMO to find motifs in rOCRs
# 4. get corresponding ChIP-seq data
# 5. filter motif by overlapping peaks
# 6. calculate TF density in each rOCRs
# 7. make figures

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/hg38_ubi-rDHS_with_gene_regulation/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/"

cd ${workDir}

# 1. get sequence
rOCRs_sequence_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/rOCRs_overlap_hg38_v28_basic_TSS.fa"

# 2. get motif list
## JASPAR core 2018 non-redundence
motif_file="/home/fankaili/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme"
grep "MOTIF" ${motif_file} | awk '{split($3,a,"(");print a[1]}' | sort -u > JASPAR_motif_list.txt
awk '{if($1=="MOTIF"){split($3,a,"(");$3=a[1];print $0}else{print $0}}' ${motif_file} > JASPAR2018_CORE_vertebrates_non-redundant2.meme

# 3. run FIMO to find motifs in rOCRs
# http://meme-suite.org/info/status?service=FIMO&id=appFIMO_5.0.51554319418644-1030637590

#--------------------------
# install MEME
cd ~/bin/
wget http://meme-suite.org/meme-software/5.0.5/meme-5.0.5.tar.gz
tar zxf meme-5.0.5.tar.gz
cd meme-5.0.5
./configure --prefix=$HOME/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
make
make test
make install
export PATH=$HOME/meme/bin:$PATH
## error
cd scripts
perl ./scripts/dependencies.pl
#--------------------------
motif="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/JASPAR2018_CORE_vertebrates_non-redundant2.meme"
rOCRs_sequence_file="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/rOCRs_overlap_hg38_v28_basic_TSS.fa"
fimo ${motif_file} ${rOCRs_sequence_file}
#
awk '{FS=OFS="\t"}{split($3,a,":");split(a[2],b,"-");print a[1],b[1],b[2],"ID",1,$6,$4,$5,$1,$2,$7,$8,$9}' ./fimo_out/fimo.tsv > rOCRs_motif_bed13.bed
sed -i '$ d' rOCRs_motif_bed13.bed
sed -i '$ d' rOCRs_motif_bed13.bed
sed -i '$ d' rOCRs_motif_bed13.bed
sed -i '$ d' rOCRs_motif_bed13.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1":"$2"-"$3]=$4}else{print $1,$2,$3,a[$1":"$2"-"$3],$7,$8,$9,$10,$11}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/rOCRs_overlap_hg38_v28_basic_TSS.bed rOCRs_motif_bed13.bed > rOCRs_motif_bed13_withID.bed
# get motif region
awk '{FS=OFS="\t"}{if(NR>1){print $1,$2+$5,$2+$6,$4,$7,$8}}' rOCRs_motif_bed13_withID.bed | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{split($6,a,"(");print $1,$2,$3,$4,a[1]}' | sort -k4,4 > motif_in_rOCRs.bed

# 4. get corresponding ChIP-seq data
python ${scriptDir}get_TF_peak_file_ENCODE.py

# 5. filter motif by overlapping peaks
# ENCSR402IDP ENCSR720HUL
sort -k3 TF_motif_peak_filelist.txt > TF_motif_peak_filelist_sort.txt
python ${scriptDir}get_merged_TF_peaks.py "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/TF_motif_peak_filelist_sort.txt"  "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/encode_TF_peak_merged.bed"
sort -k1,1 -k2,2n encode_TF_peak_merged.bed > encode_TF_peak_merged_sorted.bed
intersectBed -a motif_in_rOCRs.bed -b encode_TF_peak_merged_sorted.bed -wa -wb | awk '{if($5==$9){print $0}}' | sort -u | sort -k4,4 > motif_in_rOCRs_with_peak.txt

# 6. calculate TF density in each rOCRs
python ${scriptDir}calculate_TF_density.py "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/motif_in_rOCRs_with_peak.txt" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/rOCRs_TF_coverage.txt" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/rOCRs_TF_count.txt"
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$4]){print $0,$3-$2,a[$4],a[$4]/($3-$2)}else{print $0,$3-$2,0,0}}}' rOCRs_TF_coverage.txt /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/rOCRs_overlap_hg38_v28_basic_TSS.bed > tmp.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$4]){print $0,a[$4]}else{print $0,0}}}' rOCRs_TF_count.txt tmp.bed > rOCRs_TF_density.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0,"ubi-rOCRs"}else{print $0,"rOCRs"}}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed rOCRs_TF_density.txt > rOCRs_TF_density_withlable.txt

#### Apr18
# calculate motif desity in each rOCRs (only motif, without peak limitation)
python ${scriptDir}calculate_TF_density.py "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/motif_in_rOCRs.bed" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/rOCRs_TF_coverage_553.txt" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/TFmotif/rOCRs_TF_count_553.txt"
#
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$4]){print $0,$3-$2,a[$4],a[$4]/($3-$2)}else{print $0,$3-$2,0,0}}}' rOCRs_TF_coverage_553.txt /data/zusers/fankaili/ccre/hg38_ubi-rDHS/basic_annotation/rOCRs_overlap_hg38_v28_basic_TSS.bed > tmp.bed
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$4]){print $0,a[$4]}else{print $0,0}}}' rOCRs_TF_count_553.txt tmp.bed > rOCRs_TF_density_553.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $0,"ubi-rOCRs"}else{print $0,"rOCRs"}}}' /data/zusers/fankaili/ccre/hg38_ubi-rDHS/GRCh38_ubi-rOCRs_EDGEid.bed rOCRs_TF_density_553.txt > rOCRs_TF_density_withlable_553.txt

# 7. make figures
# Rscript make_TF_density_figures.R



awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{print $0,a[$4]}}' rOCRs_TF_density_withlable_553.txt motif_in_rOCRs.bed | sort -k5,5 -k4,4 > motif_in_rOCRs_withlabel.bed
cut -f 5 motif_in_rOCRs_withlabel.bed | sort -u > motif_list.txt
#
echo -e "motif\tubi_TF\tnon_ubi_TF\tubi\tnon_ubi" > motif_contigency_table.txt
while read motif
do
    echo ${motif}
    ubi_TF=`awk -v motif="$motif" '{if($5==motif && $6=="ubi-rOCRs"){print $0}}' motif_in_rOCRs_withlabel.bed | cut -f 4 | sort -u | wc -l`
    non_ubi_TF=`awk -v motif="$motif" '{if($5==motif && $6=="rOCRs"){print $0}}' motif_in_rOCRs_withlabel.bed | cut -f 4 | sort -u | wc -l`
    echo -e ${motif}"\t"${ubi_TF}"\t"${non_ubi_TF}"\t7543\t26635" >> motif_contigency_table.txt
done < motif_list.txt
