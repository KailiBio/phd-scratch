#!/bin/bash

# -- Kaili
# This script is for testing IDEAS using different input.
# Just using chr1.
# INPUT:
# 1. yu's input
# 2. MACS2 FE
# 3. MACS2 p-value
# 4. negBinomial


####################

# 1. yu's input
cd /data/zusers/fankaili/ideas/test_ideas/yu_input/
mkdir signal
#
for file in `ls /data/zusers/fankaili/ideas/yu_input/signal/`
do
    echo ${file};
    #paste /data/zusers/fankaili/ideas/yu_input/mm10.bed /data/zusers/fankaili/ideas/yu_input/signal/${file} | grep 'chr1' | awk '{print $5}' > ./signal/${file} ;
    paste /data/zusers/fankaili/ideas/yu_input/mm10.bed /data/zusers/fankaili/ideas/yu_input/signal/${file} | grep 'chr2' | awk '{print $5}' >> ./signal/${file} ;
done
#
grep 'chr1' /data/zusers/fankaili/ideas/yu_input/mm10.bed > ./mm10_chr1.bed
grep 'chr2' /data/zusers/fankaili/ideas/yu_input/mm10.bed >> ./mm10_chr1.bed
#
cp /data/zusers/fankaili/ideas/yu_input/IDEAS_yu_input_2.input test_ideas_1.input
sed -i 's/\/data\/zusers\/fankaili\/ideas\/yu_input\//\/data\/zusers\/fankaili\/ideas\/test_ideas\/yu_input\//g' test_ideas_1.input
#
cp /data/zusers/fankaili/ideas/yu_input/IDEAS_yu_input_2.parafile test_ideas_1.parafile
vim test_ideas_1.parafile
#
cp /data/zusers/fankaili/ideas/yu_input/IDEAS_yu_input_2.sh test_ideas_1.sh
vim test_ideas_1.sh
#
nohup bash test_ideas_1.sh > nohup.test_ideas_1.out 2>&1&



# 2. MACS2 FE
cd /data/zusers/fankaili/ideas/test_ideas/fe/
mkdir signal
#
for file in `ls /data/zusers/fankaili/ideas/run_ideas_p_value/signal/`
do
    if [[ "${file: -3}" == "txt" ]];then
        echo ${file} ;
        paste /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed /data/zusers/fankaili/ideas/run_ideas_p_value/signal/${file} | grep 'chr1' | awk '{print $5}' > ./signal/${file} ;
        paste /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed /data/zusers/fankaili/ideas/run_ideas_p_value/signal/${file} | grep 'chr2' | awk '{print $5}' >> ./signal/${file} ;
    fi
done

cp /data/zusers/fankaili/ideas/run_ideas/signal/*_200bin_group1.txt ./signal/
#
cp /data/zusers/fankaili/ideas/run_ideas/mm10_200bin_group1.bed ./
#
cp /data/zusers/fankaili/ideas/run_ideas/IDEAS_group1.input ./test_ideas_2.input
sed -i 's/\/data\/zusers\/fankaili\/ideas\/run_ideas\//\/data\/zusers\/fankaili\/ideas\/test_ideas\/fe\//g' test_ideas_2.input
#
cp /data/zusers/fankaili/ideas/run_ideas/IDEAS_group1.parafile ./test_ideas_2.parafile
vim test_ideas_2.parafile
#
cp /data/zusers/fankaili/ideas/run_ideas/IDEAS_group1.sh ./test_ideas_2.sh
vim test_ideas_2.sh
#
nohup bash test_ideas_2.sh > nohup.test_ideas_2.out 2>&1&



# 3. MACS2 p-value
cd /data/zusers/fankaili/ideas/test_ideas/pvalue/
mkdir signal
#
for file in `ls /data/zusers/fankaili/ideas/run_ideas_p_value/signal/`
do
    if [[ "${file: -3}" == "txt" ]];then
        echo ${file} ;
        paste /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed /data/zusers/fankaili/ideas/run_ideas_p_value/signal/${file} | grep 'chr1' | awk '{print $5}' > ./signal/${file} ;
        paste /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed /data/zusers/fankaili/ideas/run_ideas_p_value/signal/${file} | grep 'chr2' | awk '{print $5}' >> ./signal/${file} ;
    fi
done
#
grep "chr1" /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed > mm10_chr.bed
grep "chr2" /data/zusers/fankaili/ideas/run_ideas_p_value/mm10.bed >> mm10_chr.bed
#
cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input test_ideas_3.input
sed -i 's/\/data\/zusers\/fankaili\/ideas\/run_ideas_p_value\//\/data\/zusers\/fankaili\/ideas\/test_ideas\/pvalue\//g' test_ideas_3.input
#
cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.parafile test_ideas_3.parafile
vim test_ideas_3.parafile
#
cp /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.sh test_ideas_3.sh
vim test_ideas_3.sh
#
nohup bash test_ideas_3.sh > nohup.test_ideas_3.out 2>&1&



# 4. negBinomial
cd /data/zusers/fankaili/ideas/test_ideas/negbinomial/
mkdir signal
#
for file in `ls /data/zusers/fankaili/ideas/run_ideas_repeat/signal/`
do
    if [[ "${file: -3}" == "txt" ]];then
        echo ${file} ;
        paste /data/zusers/fankaili/ideas/run_ideas_repeat/mm10.bed /data/zusers/fankaili/ideas/run_ideas_repeat/signal/${file} | grep 'chr1' | awk '{print $5}' > ./signal/${file} ;
        paste /data/zusers/fankaili/ideas/run_ideas_repeat/mm10.bed /data/zusers/fankaili/ideas/run_ideas_repeat/signal/${file} | grep 'chr2' | awk '{print $5}' >> ./signal/${file} ;
    fi
done
#
grep "chr1" /data/zusers/fankaili/ideas/run_ideas_repeat/mm10.bed > mm10_chr.bed
grep "chr2" /data/zusers/fankaili/ideas/run_ideas_repeat/mm10.bed >> mm10_chr.bed
#
cp /data/zusers/fankaili/ideas/run_ideas_repeat/run_IDEAS_8hm_atac_dname_repeat.input test_ideas_4.input
sed -i 's/\/data\/zusers\/fankaili\/ideas\/run_ideas_repeat\//\/data\/zusers\/fankaili\/ideas\/test_ideas\/negbinomial\//g' test_ideas_4.input
sed -i 's/\/data\/zusers\/fankaili\/ideas\/run_ideas_p_value\//\/data\/zusers\/fankaili\/ideas\/test_ideas\/pvalue\//g' test_ideas_4.input
#
cp /data/zusers/fankaili/ideas/run_ideas_repeat/run_IDEAS_8hm_atac_dname_repeat.parafile test_ideas_4.parafile
vim test_ideas_4.parafile
#
cp /data/zusers/fankaili/ideas/run_ideas_repeat/run_IDEAS_8hm_atac_dname_repeat.sh test_ideas_4.sh
vim test_ideas_4.shls
#
nohup bash test_ideas_4.sh > nohup.test_ideas_4.out 2>&1&
