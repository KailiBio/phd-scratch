#!/bin/bash

# This script is for given GENCODE gtf file, get seperate annotation file.
# 1. gene info
# 2. transcript info
# 3. TSS bed file
# 4. filtered TSS bed file

# INPUT: gtf file
#              assemby info. (hg38/GRCh38, mm10 et al.)
#               gencode version, start with 'v'. (for example: v28)
#               type of gtf file, comprehensive or basic.
#               path for output file.
# OUTPUT: four files: gene, transcript, all TSS and filtered TSS.
# EXP: bash process_GENCODE_gtf.sh /home/fankaili/genome/gencode.v28.basic.annotation.gtf hg38 v28 basic
#           /home/fankaili/genome/

gtf_file=$1
assembly=$2
version=$3
gtf_type=$4
outPath=$5

#######################

#### gene
# gene gene info
awk '{FS=OFS="\t"}
{if($3=="gene"){split($9,a,";");
for(i=1;i<=length(a);i++){
    if(a[i]~/gene_id/){split(a[i],b,"\"");id=b[2]};
    if(a[i]~/gene_type/){split(a[i],c,"\"");type=c[2]};
    if(a[i]~/gene_name/){split(a[i],d,"\"");name=d[2]}
};
print $1,$4,$5,id,$7,type,name}}' ${gtf_file} > ${outPath}${assembly}_${version}_${gtf_type}_gene.txt
sed -i '1ichr\tstart\tend\tgene_id\tstrand\tgene_type\tgene_name' \
${outPath}${assembly}_${version}_${gtf_type}_gene.txt
# filtering gene
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(FNR>1 && a[$6]){print $0}}}' \
/home/fankaili/genome/transcript_type_for_filter_TSS.txt ${outPath}${assembly}_${version}_${gtf_type}_gene.txt > ${outPath}${assembly}_${version}_${gtf_type}_gene_filtered.txt
sed -i '1ichr\tstart\tend\tgene_id\tstrand\tgene_type\tgene_name' \
${outPath}${assembly}_${version}_${gtf_type}_gene_filtered.txt
# get protein-coding genes
awk '{FS=OFS="\t"}{if(NR==1 || $6=="protein_coding"){print $0}}' ${outPath}${assembly}_${version}_${gtf_type}_gene.txt > ${outPath}${assembly}_${version}_${gtf_type}_gene_protein_coding.txt

#### transcript
# get transcript info
awk '{FS=OFS="\t"}
{if($3=="transcript"){split($9,a,";");
for(i=1;i<=length(a);i++){
    if(a[i]~/gene_id/){split(a[i],b,"\"");id=b[2]};
    if(a[i]~/gene_type/){split(a[i],c,"\"");type=c[2]};
    if(a[i]~/gene_name/){split(a[i],d,"\"");name=d[2]};
    if(a[i]~/transcript_id/){split(a[i],e,"\"");id2=e[2]};
    if(a[i]~/transcript_type/){split(a[i],f,"\"");type2=f[2]};
    if(a[i]~/transcript_name/){split(a[i],g,"\"");name2=g[2]};
};
print $1,$4,$5,id2,$7,id,type,name,type2,name2}}' ${gtf_file}  > ${outPath}${assembly}_${version}_${gtf_type}_transcript.txt
sed -i '1ichr\tstart\tend\ttranscript_id\tstrand\tgene_id\tgene_type\tgene_name\ttranscript_type\ttranscript_name' \
${outPath}${assembly}_${version}_${gtf_type}_transcript.txt
# filtering transcript
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(FNR==1){print $0}else{if(a[$9]){print $0}}}}' /home/fankaili/genome/transcript_type_for_filter_TSS.txt ${outPath}${assembly}_${version}_${gtf_type}_transcript.txt > ${outPath}${assembly}_${version}_${gtf_type}_transcript_filtered.txt
# transcripts of protein-coding gene
awk '{FS=OFS="\t"}{if(NR==1 || $9=="protein_coding"){print $0}}' ${outPath}${assembly}_${version}_${gtf_type}_transcript.txt > ${outPath}${assembly}_${version}_${gtf_type}_transcript_protein_coding.txt

#### TSS
# get TSS
awk '{FS=OFS="\t"}{if(NR>1){if($5=="+"){print $1,$2,$2,$4,".",$5,$6}else{print $1,$3,$3,$4,".",$5,$6}}}' ${outPath}${assembly}_${version}_${gtf_type}_transcript.txt | \
sort -k1,1 -k2,2n > ${outPath}${assembly}_${version}_${gtf_type}_TSS.bed
# filtering TSS
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(FNR>1 && a[$9]){if($5=="+"){print $1,$2,$2,$4,".",$5,$6}else{print $1,$3,$3,$4,".",$5,$6}}}}' \
/home/fankaili/genome/transcript_type_for_filter_TSS.txt ${outPath}${assembly}_${version}_${gtf_type}_transcript.txt \
| sort -k1,1 -k2,2n > ${outPath}${assembly}_${version}_${gtf_type}_TSS_filtered.bed
# TSSs of protein-coding gene
awk '{FS=OFS="\t"}{if(NR>1){if($5=="+"){print $1,$2,$2,$4,".",$5,$6}else{print $1,$3,$3,$4,".",$5,$6}}}' ${outPath}${assembly}_${version}_${gtf_type}_transcript_protein_coding.txt > ${outPath}${assembly}_${version}_${gtf_type}_TSS_protein_coding.bed
