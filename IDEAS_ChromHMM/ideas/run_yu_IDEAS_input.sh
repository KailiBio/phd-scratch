#!/bin/bash

# -- Kaili
# This script is for geting Yu's IDEAS input and run IDEAS.
# INPUT
# OUTPUT

# 1. get bed file & make input file
for line in `awk '{print $2}' /data/projects/segmentation/IDEAS_Input/md5sum.txt`
do
    path="/data/projects/segmentation/IDEAS_Input/"${line#./};
    filename=${line#./bx.psu.edu/~yuzhang/me66/Input/};
    echo $filename;
    gunzip -c ${path} > /data/zusers/fankaili/ideas/yu_input/signal/${filename%.gz} ;
    tissue=`cut -d "_" -f1,2 <<< ${filename}` ;
    mark=`echo ${filename} | cut -d "_" -f3 | cut -d "." -f2` ;
    if [ ${mark} = 'cpm' ]; then
        mark="ATAC" ;
    elif [[ $mark == ENCSR* ]]; then
        mark="DNAme" ;
    fi
    echo $mark ;
    echo ${tissue}" "${mark}" /data/zusers/fankaili/ideas/yu_input/signal/"${filename%.gz} >> /data/zusers/fankaili/ideas/yu_input/IDEAS_yu_input.input
done

# 2. make .sh & .parafile
vim IDEAS_yu_input.sh
vim IDEAS_yu_input.parafile

cp -r /home/fankaili/git/IDEAS_2018/bin /data/zusers/fankaili/ideas/yu_input/
cp -r /home/fankaili/git/IDEAS_2018/data /data/zusers/fankaili/ideas/yu_input/

# 3. make
### @ z018
vim run_IDEAS_yu_input.sh
sbatch run_IDEAS_yu_input.sh
