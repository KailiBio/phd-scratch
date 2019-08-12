#!/bin/bash

# -- Kaili
# This script is for running ChromImpute to do CTCF imputation.

cd /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/

chromimputPath="/home/fankaili/bin/ChromImpute/"

master_table="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute/master_table4.txt"
chrom_size="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute/mm10_chr19_chrom.size"

dataDir="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute/input_data"
#convertDir="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute/converted_data"
distanceDir="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute/distance4/"
trainDataDir="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute/train_data4/"
predictorDir="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute/predictor4/"
outputDir="/data/zusers/fankaili/ideas/imputation_comparison/chromimpute/output4/"

SECONDS=0
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
java -mx4000M -jar ${chromimputPath}ChromImpute.jar Train ${trainDataDir} ${master_table} ${predictorDir} forebrain_0 CTCF
java -mx4000M -jar ${chromimputPath}ChromImpute.jar Train ${trainDataDir} ${master_table} ${predictorDir} midbrain_0 CTCF
java -mx4000M -jar ${chromimputPath}ChromImpute.jar Train ${trainDataDir} ${master_table} ${predictorDir} hindbrain_0 CTCF
echo ""
echo "finish step4: Train"
echo ""

# 5. generate the imputed signal track
## Apply
java -mx4000M -jar ${chromimputPath}ChromImpute.jar Apply ${dataDir} ${distanceDir} ${predictorDir} ${master_table} ${chrom_size} ${outputDir} forebrain_0 CTCF
java -mx4000M -jar ${chromimputPath}ChromImpute.jar Apply ${dataDir} ${distanceDir} ${predictorDir} ${master_table} ${chrom_size} ${outputDir} midbrain_0 CTCF
java -mx4000M -jar ${chromimputPath}ChromImpute.jar Apply ${dataDir} ${distanceDir} ${predictorDir} ${master_table} ${chrom_size} ${outputDir} hindbrain_0 CTCF
echo ""
echo "finish step5: Apply"
echo ""

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
