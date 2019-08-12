#!/bin/bash

# -- Kaili
# This script is for doing imputation only.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

# 0. preparation
## 1) get bed file order
merge_state_cluster_file(){
    prefix=$1
    bed_file=$2
    #
    awk '{print $1}' ${bed_file} | uniq > bed.order
    # merge state file
    if [ -f ./${prefix}_result/${prefix}state ];then rm ./${prefix}_result/${prefix}state; fi
    for i in {1..19} X Y
    do
        cat ./${prefix}_result/${prefix}.chr${i}.state | awk '{if(NR>1){print $0}}' >> ./${prefix}_result/${prefix}state
    done
    awk '{FS=OFS}{if(NR==FNR){a[$1]=$0}else{print a[$4]}}' ./${prefix}_result/${prefix}state ${bed_file} > ./${prefix}_result/tmp.state
    head -1 ./${prefix}_result/${prefix}.chr1.state > ./${prefix}_result/${prefix}.state
    cat ./${prefix}_result/tmp.state >> ./${prefix}_result/${prefix}.state
    # merge cluster
    if [ -f ./${prefix}_result/${prefix}.cluster ]; then rm ./${prefix}_result/${prefix}.cluster; fi
    while read chr
    do
        cat ./${prefix}_result/${prefix}.${chr}.cluster >> ./${prefix}_result/${prefix}.cluster
    done < bed.order
}

prefix="dhs_ctcf"
bed_file="mm10_OCR-center_bins_v3_signal_based_space.bed"
merge_state_cluster_file ${prefix} ${bed_file}

## 2) change the code to custom
## /data/zusers/fankaili/ideas/dhs_ctcf/bin/ideaspipe.R
# str=paste("ideas", paste(args, collapse=" "), "-thread", threadn, "-o", out, "-otherstate /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf.state /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf.cluster")
# str=paste("ideas", paste(args, collapse=" "), "-thread", threadn, "-o", out, "-otherstate /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf.state /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf.cluster", "-K 1")


# 1. 9-impute-11, 100
cp ctcf_9sample_impute_11sample.parafile 9-impute-11_only_100.parafile
vim 9-impute-11_only_100.parafile
cp ctcf_9sample_impute_11sample.sh 9-impute-11_only_100.sh
vim 9-impute-11_only_100.sh
cp ctcf_9sample_impute_11sample.input 9-impute-11_only_100.input

nohup bash 9-impute-11_only_100.sh > ./nohup/nohup.9-impute-11_only_100.out 2>&1&

# 2. 9-impute-11, 10K
awk '{FS=OFS}{if(NR==1){print $0}else{printf 10000*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 10000*$i};printf "\n"}}' /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf.para0 > /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf_10k.para0

cp 9-impute-11_only_100.parafile 9-impute-11_only_10k.parafile
vim 9-impute-11_only_10k.parafile
cp 9-impute-11_only_100.sh 9-impute-11_only_10k.sh
vim 9-impute-11_only_10k.sh
cp 9-impute-11_only_100.input 9-impute-11_only_10k.input

nohup bash 9-impute-11_only_10k.sh > ./nohup/nohup.9-impute-11_only_10k.out 2>&1&

# 3. 9-impute-11, 1m
awk '{FS=OFS}{if(NR==1){print $0}else{printf 1000000*$1;for(i=2;i<=NF;i++){printf " ";printf "%.2f", 1000000*$i};printf "\n"}}' /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf.para0 > /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf_1m.para0

cp 9-impute-11_only_100.parafile 9-impute-11_only_1m.parafile
vim 9-impute-11_only_1m.parafile
cp 9-impute-11_only_100.sh 9-impute-11_only_1m.sh
vim 9-impute-11_only_1m.sh
cp 9-impute-11_only_100.input 9-impute-11_only_1m.input

nohup bash 9-impute-11_only_1m.sh > ./nohup/nohup.9-impute-11_only_1m.out 2>&1&






## 100
./bin/ideas 9-impute-11_only_100.input.1 mm10_OCR-center_bins_v3_signal_based_space.bed -impute none -norm -C 100 -minerr 0.5 -cap 16 -otherpara /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf_100fold.para0 -sample 20 5 -thread 32 -o /data/zusers/fankaili/ideas/dhs_ctcf/9-impute-11_only_100_result/9-impute-11_only_100.tmp.1 -otherstate /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf.state -G


## 100 & one-cluster
./bin/ideas 9-impute-11_only_100.input.1 mm10_OCR-center_bins_v3_signal_based_space.bed -impute none -norm -C 100 -minerr 0.5 -cap 16 -otherpara /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf_100fold.para0 -sample 20 5 -thread 32 -o /data/zusers/fankaili/ideas/dhs_ctcf/9-impute-11_only_100_result/9-impute-11_only_100.tmp.2 -otherstate /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf.state -K1


## 10k
./bin/ideas 9-impute-11_only_100.input.1 mm10_OCR-center_bins_v3_signal_based_space.bed -impute none -norm -C 100 -minerr 0.5 -cap 16 -otherpara /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf_10k.para0 -sample 20 5 -thread 32 -o /data/zusers/fankaili/ideas/dhs_ctcf/9-impute-11_only_10k_result/9-impute-11_only_10k.tmp.1 -otherstate /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf.state

## 1m
./bin/ideas 9-impute-11_only_100.input.1 mm10_OCR-center_bins_v3_signal_based_space.bed -impute none -norm -C 100 -minerr 0.5 -cap 16 -otherpara /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf_1m.para0 -sample 20 5 -thread 32 -o /data/zusers/fankaili/ideas/dhs_ctcf/9-impute-11_only_1m_result/9-impute-11_only_1m.tmp.1 -otherstate /data/zusers/fankaili/ideas/dhs_ctcf/dhs_ctcf_result/dhs_ctcf.state


#################
# -K 1
cp ctcf_9sample_impute_11sample.sh ctcf_9sample_impute_11sample_K.sh
vim ctcf_9sample_impute_11sample_K.sh
cp ctcf_9sample_impute_11sample.parafile ctcf_9sample_impute_11sample_K.parafile
vim ctcf_9sample_impute_11sample_K.parafile
cp ctcf_9sample_impute_11sample.input ctcf_9sample_impute_11sample_K.input

nohup bash ctcf_9sample_impute_11sample_K.sh > ./nohup/nohup.ctcf_9sample_impute_11sample_K.out 2>&1&




# bash validate_imputation-only.sh
