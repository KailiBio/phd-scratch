#!/bin/bash

# -- Kaili
# This script is for generating signal in matched resolution for ChromImpute and Avocado.

mark=$1

workDir="/data/zusers/fankaili/ideas/imputation_comparison/"
oldSignalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
signalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_chr19_25bp_dhsResolution/"

cd ${workDir}

####
if [ $mark == "DNAme" ]; then col=6; else col=5; fi

for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    awk -v col="$col" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$col}else{print $1,$2,$3,a[$4]}}' ${oldSignalDir}${sample}_${mark}_dhs.tab mm10_OCR_chr19_bins.bed | sort -k1,1 -k2,2n > tmp.${mark}.bedGraph
    bedGraphToBigWig tmp.${mark}.bedGraph /home/fankaili/genome/mm10.chrom.sizes.clean ${signalDir}${sample}_${mark}_dhs.bw
    bigWigAverageOverBed ${signalDir}${sample}_${mark}_dhs.bw mm10_chr19_25bp.bed ${signalDir}${sample}_${mark}_dhs.tab
    awk -v col="$col" '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$col}else{print a[$4]}}' ${signalDir}${sample}_${mark}_dhs.tab mm10_chr19_25bp.bed > ${signalDir}${sample}_${mark}_dhs.txt
done

rm tmp.${mark}.bedGraph
