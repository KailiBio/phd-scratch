#!/bin/bash

# -- Kaili
# This script is for running ChromImpute to do CTCF imputation.

workDir=$1
chrom_size_file=$2

cd ${workDir}

master_table=${workDir}"master_table.txt"

dataDir=${workDir}"input_data"
#convertDir="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute/converted_data"
distanceDir=${workDir}"distance/"
trainDataDir=${workDir}"train_data/"
predictorDir=${workDir}"predictor/"
outputDir=${workDir}"output/"

#########

SECONDS=0
echo ""
echo "Start here!"
echo ""

# 1. convert data from 25bp to desired resolution
## Convert
# java -mx4000M -jar ${chromimputPath}ChromImpute.jar Convert ${dataDir} ${master_table} ${chrom_size} ${convertDir}
# echo ""
# echo "finish step1: Convert"
# echo ""

# 2. calculate global distance between datasets
## ComputeGlobalDist
java -mx4000M -jar ${chromimputPath}ChromImpute.jar ComputeGlobalDist ${dataDir} ${master_table} ${chrom_size} ${distanceDir}
echo ""
echo "finish step2: ComputeGlobalDist"
echo ""

# 3. generate the features for the training
## GenerateTrainData
java -mx4000M -jar ${chromimputPath}ChromImpute.jar GenerateTrainData ${dataDir} ${distanceDir} ${master_table} ${chrom_size} ${trainDataDir} CTCF
echo ""
echo "finish step3: GenerateTrainData"
echo ""

# 4. generate the trianed predictors for a specific mark in a specific sample type of interest
## Train
java -mx4000M -jar ${chromimputPath}ChromImpute.jar Train ${trainDataDir} ${master_table} ${predictorDir} liver_14.5 CTCF
java -mx4000M -jar ${chromimputPath}ChromImpute.jar Train ${trainDataDir} ${master_table} ${predictorDir} lung_14.5 CTCF
echo ""
echo "finish step4: Train"
echo ""

# 5. generate the imputed signal track
## Apply
java -mx4000M -jar ${chromimputPath}ChromImpute.jar Apply ${dataDir} ${distanceDir} ${predictorDir} ${master_table} ${chrom_size} ${outputDir} liver_14.5 CTCF
java -mx4000M -jar ${chromimputPath}ChromImpute.jar Apply ${dataDir} ${distanceDir} ${predictorDir} ${master_table} ${chrom_size} ${outputDir} lung_14.5 CTCF
echo ""
echo "finish step5: Apply"
echo ""

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
