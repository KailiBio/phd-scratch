#!/bin/bash

# -- Kaili
# This script is for cutting code.sh into array jobs.

# INPUT: your command file (normal .sh file)
#              path for separate code file
#              name prefix of separate code file
#              num of separate files you want. (i)
#              the maximum command line in each file. (n; total command line should >= i*n)
# OUTPUT: separate command into i file, then run it.
# EXP: bash cut_code_into_array.sh /data/zusers/fankaili/ideas/get_DNAme_bw_code.sh /data/zusers/fankaili/ideas/
#          get_DNAme_bw_code 11 24

codeFile=$1
codePath=$2
prefix=$3
num_of_files=$4
num_of_command=$5

#################

for ((i=1; i<=num_of_files; i++))
do
    awk -v i="$i" -v n="$num_of_command" '{if(NR>((i-1)*n) && NR<=(n*i)){print $0}}' ${codeFile} \
    > ${codePath}${prefix}_${i}.sh
    nohup bash ${codePath}${prefix}_${i}.sh > ${codePath}nohup.${prefix}_${i}.out 2>&1&
done
