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


SECONDS=0

set -e

if [ -d /tmp/${prefix} ];then rm -rf /tmp/${prefix}/; fi
mkdir /tmp/${prefix}/
cd /tmp/${prefix}/

MAXWAIT=600 #MAXWAIT is the maximum wait second
sleep $((RANDOM % MAXWAIT))

cp ${read1} /tmp/${prefix}/${prefix}.fastq

# executables
STAR=STAR
RSEM=rsem-calculate-expression
bedGraphToBigWig=bedGraphToBigWig

# STAR parameters: common
STARparCommon=" --genomeDir $STARgenomeDir  --readFilesIn /tmp/${prefix}/${prefix}.fastq   --outSAMunmapped Within --outFilterType BySJout \
 --outSAMattributes NH HI AS NM MD    --outFilterMultimapNmax 20   --outFilterMismatchNmax 999   \
 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   \
 --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1"

# STAR parameters: run-time, controlled by DCC
STARparRun=" --runThreadN $nThreadsSTAR --genomeLoad LoadAndKeep  --limitBAMsortRAM 10000000000"

# STAR parameters: type of BAM output: quantification or sorted BAM or both
#     OPTION: sorted BAM output && transcritomic BAM for quantification
STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"


# STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM

case "$dataType" in
str_SE|str_PE)
      #OPTION: stranded data
      STARparStrand=""
      STARparWig="--outWigStrand Stranded"
      ;;
      #OPTION: unstranded data
unstr_SE|unstr_PE)
      STARparStrand="--outSAMstrandField intronMotif"
      STARparWig="--outWigStrand Unstranded"
      ;;
esac

# STAR parameters: metadata
STARparsMeta="--outSAMheaderCommentFile commentsENCODElong.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate"



###### STAR command
echo $STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta
$STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta

#--------------
echo "step 1, all set!"

###### bedGraph generation, now decoupled from STAR alignment step
# working subdirectory for this STAR run
mkdir Signal

echo $STAR --runMode inputAlignmentsFromBAM   --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
$STAR --runMode inputAlignmentsFromBAM   --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr

# move the signal files from the subdirectory
mv Signal/Signal*bg .

#--------------
echo "step 2, all set!"

###### bigWig conversion commands
# exclude spikeins
grep ^chr $STARgenomeDir/chrNameLength.txt > chrNL.txt

case "$dataType" in
str_SE|str_PE)
      # stranded data
      str[1]=-; str[2]=+;
      for istr in 1 2
      do
      for imult in Unique UniqueMultiple
      do
          grep ^chr Signal.$imult.str$istr.out.bg | sort -k1,1 -k2,2n > sig.tmp
          $bedGraphToBigWig sig.tmp  chrNL.txt Signal.$imult.strand${str[istr]}.bw
      done
      done
      ;;
unstr_SE|unstr_PE)
      # unstranded data
      for imult in Unique UniqueMultiple
      do
          grep ^chr Signal.$imult.str1.out.bg | sort -k1,1 -k2,2n > sig.tmp
          $bedGraphToBigWig sig.tmp chrNL.txt  Signal.$imult.unstranded.bw
      done
      ;;
esac

#--------------
echo "step 3, all set!"

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
echo "step 4, all set!"

#################################
# Jul16
# add size sanity check
###
RED='\033[0;31m'
NC='\033[0m' # No Color
echo ""
echo "### sanity check!"
#
# check if the number of reads in STAR output are the same as STAR log
x=`grep "Number of input reads" Log.final.out | cut -f 2`
echo "number of input reads in STAR log: "$x
y=`samtools view Aligned.sortedByCoord.out.bam | cut -f 1 | sort -u | wc -l`
echo "number of reads in bam file: "$y
if [ "$x" != "$y" ];then
	echo -e "${RED}bam file not match/! Error\!\!\! ${NC}";
	#rm -rf /tmp/${prefix};
	#exit 1;
else
	echo "bam file matched. Continue...";
fi
# check if the bam file are truncated
sanity_check=`samtools quickcheck -v Aligned.toTranscriptome.out.bam`
if [ "$sanity_check" == "Aligned.toTranscriptome.out.bam" ];then
	echo -e "${RED}Truncated file/! Error\!\!\! ${NC}";
	#rm -rf /tmp/${prefix};
	#exit 2;
else
	echo "transcript bam file is intact. Continue...";
fi
#
sanity_check2=`samtools quickcheck -v Aligned.toTranscriptome.out.sorted.bam`
if [ "$sanity_check2" == "Aligned.toTranscriptome.out.sorted.bam" ];then
	echo -e "${RED}Truncated sorted ile/! Error\!\!\! ${NC}";
	#rm -rf /tmp/${prefix};
	#exit 3;
else
	echo "sorted transcript bam file is intact. Continue...";
fi
echo ""




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
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir ${prefix} >& Log.rsem
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.sorted.bam $RSEMrefDir ${prefix} >& Log.rsem

#--------------
echo "step 5, all set!"


############################
# Jul16
# check RSEM log
RED='\033[0;31m'
NC='\033[0m' # No Color
#
log_line=`awk 'NR==2' Log.rsem`
if [ "$log_line" != "Parsed 1000000 entries" ];
then 
	echo -e "${RED}Warning: Need double check!${NC}"; 
else
	echo "RSEM seems correct."
fi



###### RSEM diagnostic plot creation
# Notes:
# 1. rsem-plot-model requires R (and the Rscript executable)
# 2. This command produces the file Quant.pdf, which contains multiple plots
echo rsem-plot-model ${prefix} ${prefix}.pdf
rsem-plot-model ${prefix} ${prefix}.pdf


### move and clean folder
rm /tmp/${prefix}/${prefix}.fastq
mv /tmp/${prefix}/* ${workDir}
rm -rf /tmp/${prefix}

### running time 
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

# ----
echo "Done. Cheers!"
