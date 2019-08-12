#!/bin/bash

# -- Kaili
# This script is for converting signal file into zipped wig for running ChromImpute.

sample=$1

signalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_25bp_bins/"
ordered_bed_file="/data/zusers/fankaili/ideas/imputation_comparison/mm10_25bp_inOrder.bed"
outDir="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute_wholeGenome/input_data/"
outDir_unzip="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute_wholeGenome/input_data2/"
chrom_size_file="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute_wholeGenome/mm10.chrom.sizes.clean2"

for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC CTCF DNAme
do
    name=${sample}"_"${mark}"_25bp"
    #
    paste ${ordered_bed_file} ${signalDir}${name}.txt > tmp.convert_zippedWig.txt
    #
    echo "track type=wiggle_0 name="${name} > ${outDir_unzip}${name}_25bp.wig
    for c in {1..19} X Y
    do
        echo "fixedStep  chrom=chr"${c}" start=0 step=25 span=25" >> ${outDir_unzip}${name}_25bp.wig
        max=`awk -v chr="$c" '{if($1=="chr"chr){print $2}}' ${chrom_size_file}`
        awk -v max="${max}" -v chr="$c" '{if($1=="chr"chr){m=$2/25;a[m]=1;b[m]=$5}}END{for(i=1;i<=(max/25);i++){if(a[i]){print b[i]}else{print 0}}}' tmp.convert_zippedWig.txt >> ${outDir_unzip}${name}_25bp.wig
    done
    #
    gzip -c ${outDir_unzip}${name}_25bp.wig > ${outDir}chr19_${name}_25bp.wig.gz
done
