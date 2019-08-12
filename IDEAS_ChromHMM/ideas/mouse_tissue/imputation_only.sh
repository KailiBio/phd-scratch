#!/bin/bash

# -- Kaili
# This script is for running IDEAS for imputation only.


##############
# shell

IDEAS_job_name=DHS_reps
script_dir=/home/fankaili/git/IDEAS_2018/
working_dir=/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/
output_dir=/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_result/
binfile=mm10_OCR-center_bins_v3_signal_based_space.bed

time Rscript ./bin/runme.R $IDEAS_job_name'.input' $IDEAS_job_name'.parafile' $output_dir

# step 1. get parafile, put into command
Rscript bin/ideaspipe.R DHS_reps.input mm10_OCR-center_bins_v3_signal_based_space.bed -o /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_result/DHS_reps -sample 20 5 -thread 32 -randstart 50 500000 -split mm10_OCR-center_bins_v3_signal_based_space.bed.inv -impute None -norm -G 42 -C 100 -minerr 0.5 -cap 16

# step 2. run ideaspipe.R
Rscript bin/ideaspipe.R DHS_reps.input  -o /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_result/DHS_reps -sample 20 5 -thread 32 -randstart 50 500000 -split mm10_OCR-center_bins_v3_signal_based_space.bed.inv -impute None -norm -G 42 -C 100 -minerr 0.5 -cap 16

## R
runideas(targs, tout, head, tinv);
## add "-otherstate .state .cluster"
## add "-K 1"

ideas DHS_reps.input -impute none -norm -G 42 -C 100 -minerr 0.5 -cap 16 -sample 20 5 -thread 32 -o  /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_result/DHS_reps.tmp.1

ideas DHS_reps.input -impute none -norm -G 42 -C 100 -minerr 0.5 -cap 16 -sample 20 5 -thread 32 -o  /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_result/DH
S_reps.tmp.1 -inv 0 1320


##############
# R

datafile = "DHS_reps.input"
parafile = "DHS_reps.parafile"
tmpfolder = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_result/"

args=c("DHS_reps.input", "mm10_OCR-center_bins_v3_signal_based_space.bed", "-o", "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_result/DHS_reps" "-sample", "20", "5", "-thread", "32", "-randstart", "50", "500000", "-split", "mm10_OCR-center_bins_v3_signal_based_space.bed.inv", "-impute", "None", "-norm", "-G", "42", "-C", "100", "-minerr", "0.5", "-cap", "16")


###
args=c("DHS_reps.input",  "-o"," /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/DHS_reps_result/DHS_reps", "-sample", "20", "5", "-thread", "32", "-randstart", "50", "500000", "-split", "mm10_OCR-center_bins_v3_signal_based_space.bed.inv", "-impute", "None", "-norm", "-G", "42", "-C", "100", "-minerr", "0.5", "-cap", "16")

runideas(targs, tout, head, tinv)
