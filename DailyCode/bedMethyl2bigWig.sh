#!/bin/bash

# -- Kaili
# This script is for converting CpG bedMethyl file from Bismark to bigWig.
### CpG sites with read coverage lower than 5 are removed.
### signal value are methylated percentage.
### need bedGraphToBigWig.

# INPUT: CpG bedMethylfile. (or CHH/CHG bedMethyl file, if you want.)
#               chrom.size file.
#               the path&name of output file. (it would be better to end with .bigWig or .bw)
# OUTPUT: bigWig file in 1bp resolution.
# EXP: bash bedMethyl2bigWig.sh ENCFF306DKT.bed /home/fankaili/genome/mm10.chrom.sizes mm10_embryonic-facial-prominence_11.5_DNAme.bigWig

bedFile=$1
chromSize=$2
output=$3

############
awk '{FS=OFS="\t"}{if($10>5){print $1,$2,$3,$11/100}}' ${bedFile} | sort -k1,1 -k2,2n > tmp_bedMethyl2bigWig_${bedFile}.bedGraph
bedGraphToBigWig tmp_bedMethyl2bigWig_${bedFile}.bedGraph ${chromSize} ${output}
rm tmp_bedMethyl2bigWig_${bedFile}.bedGraph
