#!/bin/bash

# -- Kaili
# This script is for impuating CTCF states.
# 1. run 66 samples: 10 marks + CTCF
# 2. run 11 samples all with CTCF
# 3. rerun 11 samples 11 marks only for rep1 data

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
signalDir="/data/zusers/fankaili/ideas/signal/macs2_pvalue/"

mkdir /data/zusers/fankaili/ideas/CTCF_impute/
workDir="/data/zusers/fankaili/ideas/CTCF_impute/"

cd ${workDir}

# 1. run 66 samples: 10 marks + CTCF
mkdir ${workDir}all_samples
cd ${workDir}all_samples

# .input file
cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input all_samples_CTCF.input
grep "CTCF" /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input >> all_samples_CTCF.input

# .sh file
cp /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/e14.5p0_8hm_ATAC_DNAme.sh all_samples_CTCF.sh
vim all_samples_CTCF.sh

# .parafile file
cp /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/e14.5p0_8hm_ATAC_DNAme.parafile all_samples_CTCF.parafile
vim all_samples_CTCF.parafile

# other
cp /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/mm10.bed ./
cp -r /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/bin ./
cp -r /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/data ./

nohup bash all_samples_CTCF.sh > ./nohup.all_samples_CTCF.out 2>&1&



# 2. run 11 samples all with CTCF
mkdir ${workDir}ctcf_samples
cd ${workDir}ctcf_samples

# .input file
grep "CTCF" /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input | \
awk '{FS=OFS=" "}{print $1}' > CTCF_sample_list.txt
#
awk '{FS=OFS=" "}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' CTCF_sample_list.txt \
/data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme_CTCF_1/e14.5p0_8hm_ATAC_DNAme_CTCF_1.input > ctcf_samples.input

# .sh file
cp ${workDir}all_samples/all_samples_CTCF.sh ctcf_samples.sh
vim ctcf_samples.sh

# .parafile file
cp ${workDir}all_samples/all_samples_CTCF.parafile ctcf_samples.parafile
vim ctcf_samples.parafile

# other
cp /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/mm10.bed ./
cp -r /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/bin ./
cp -r /data/zusers/fankaili/ideas/e14.5p0_8hm_ATAC_DNAme_CTCF/e14.5p0_8hm_ATAC_DNAme/data ./

nohup bash ctcf_samples.sh > ./nohup.ctcf_samples.out 2>&1&



# 3. rerun 11 samples 11 marks only for rep1 data
# .input file
cp ctcf_samples.input ctcf_samples_old.input
awk '{FS=OFS=" "}{print $1,$2,"/data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/"$1"_"$2"_normal.txt"}' \
ctcf_samples_old.input > ctcf_samples.input
# bed file
cut -f 1 /data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/embryonic-facial-prominence_11.5_H3K36me3_normal.tab > tmp.txt
awk '{FS=OFS=" "}{if(NR==FNR){a[$4]=$0}else{print a[$1]}}' mm10.bed tmp.txt > tmp2.txt
mv tmp2.txt mm10.bed
#
nohup bash ctcf_samples.sh > ./nohup.ctcf_samples_Jan17.out 2>&1&
