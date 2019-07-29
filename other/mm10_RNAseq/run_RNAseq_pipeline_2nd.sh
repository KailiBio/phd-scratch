#!/bin/bash

# -- Kaili
# This script is for running STAR & RSEM to calculate gene expression.

read1=$1 #gzipped fastq file for read1
read2=$2 #gzipped fastq file for read1, use "" if single-end
STARgenomeDir=$3
RSEMrefDir=$4
dataType=$5 # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE
nThreadsSTAR=$6 # number of threads for STAR
nThreadsRSEM=$7 # number of threads for RSEM
workDir=$8
prefix=$9

# output: all in the working directory, fixed names
# Aligned.sortedByCoord.out.bam                 # alignments, standard sorted BAM, agreed upon formatting
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Quant.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
# Quant.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
# Quant.pdf                                     # RSEM diagnostic plots
# Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data
# Signal.{Unique,UniqueMultiple}.unstranded.bw  # 2 bigWig files for unstranded data

# read1="/data/zusers/fankaili/ccre/mm10_rnaseq/ENCFF581SPK.fastq"
# STARgenomeDir="/home/fankaili/STARgenome"
# RSEMrefDir="/home/fankaili/RSEMgenome/RSEMref_mm10"
# dataType="unstr_SE" # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE
# nThreadsSTAR=8 # number of threads for STAR
# nThreadsRSEM=8 # number of threads for RSEM
# workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/ss/"
# prefix="ss"

cd ${workDir}

# executables
STAR=STAR
RSEM=rsem-calculate-expression
bedGraphToBigWig=bedGraphToBigWig

SECONDS=0

######### RSEM
#### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
trBAMsortRAM=60G

#mv Aligned.toTranscriptome.out.bam Tr.bam

case "$dataType" in
str_SE|unstr_SE)
      # single-end data
      cat <( samtools view -H Aligned.toTranscriptome.out.bam ) <( samtools view -@ $nThreadsRSEM Aligned.toTranscriptome.out.bam | sort -S $trBAMsortRAM -T ./ ) | samtools view -@ $nThreadsRSEM -bS - > Aligned.toTranscriptome.out.sorted.bam
      ;;
str_PE|unstr_PE)
      # paired-end data, merge mates into one line before sorting, and un-merge after sorting
      cat <( samtools view -H Aligned.toTranscriptome.out.bam ) <( samtools view -@ $nThreadsRSEM Aligned.toTranscriptome.out.bam | awk '{printf "%s", $0 " "; getline; print}' | sort -S $trBAMsortRAM -T ./ | tr ' ' '\n' ) | samtools view -@ $nThreadsRSEM -bS - > Aligned.toTranscriptome.out.sorted.bam
      ;;
esac

#'rm' Tr.bam

#--------------
echo "4"

# RSEM parameters: common
RSEMparCommon="--bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345"

# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun=" -p $nThreadsRSEM --ci-memory 30000 "

# RSEM parameters: data type dependent

case "$dataType" in
str_SE)
      #OPTION: stranded single end
      RSEMparType="--forward-prob 0"
      ;;
str_PE)
      #OPTION: stranded paired end
      RSEMparType="--paired-end --forward-prob 0"
      ;;
unstr_SE)
      #OPTION: unstranded single end
      RSEMparType=""
      ;;
unstr_PE)
      #OPTION: unstranded paired end
      RSEMparType="--paired-end"
      ;;
esac


###### RSEM command
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.sorted.bam $RSEMrefDir ${prefix} >& Log.rsem
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.sorted.bam $RSEMrefDir ${prefix} >& Log.rsem

#--------------
echo "5"

###### RSEM diagnostic plot creation
# Notes:
# 1. rsem-plot-model requires R (and the Rscript executable)
# 2. This command produces the file Quant.pdf, which contains multiple plots
echo rsem-plot-model ${prefix} ${prefix}.pdf
rsem-plot-model ${prefix} ${prefix}.pdf


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

# ----
echo "Done. Cheers!"
