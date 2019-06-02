#!/bin/bash

# -- Kaili
# This script is for runing IDEAS for testing imputation between reps using partial model.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

# 1. ctcf_9sample_to_11sample_rep2
cp ctcf_9sample_to_11sample.sh ctcf_9sample_to_11sample_rep2.sh
vim ctcf_9sample_to_11sample_rep2.sh
cp ctcf_9sample_to_11sample.parafile ctcf_9sample_to_11sample_rep2.parafile
vim ctcf_9sample_to_11sample_rep2.parafile
cp ctcf_9sample_to_11sample.input ctcf_9sample_to_11sample_rep2.input
vim ctcf_9sample_to_11sample_rep2.input

nohup bash ctcf_9sample_to_11sample_rep2.sh > ./nohup/nohup.ctcf_9sample_to_11sample_rep2.out 2>&1&
