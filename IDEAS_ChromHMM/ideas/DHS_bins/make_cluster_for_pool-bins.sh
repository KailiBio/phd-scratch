#!/bin/bash

# -- Kaili
# This script is for making clusering figures and merging.

mark=$1

mkdir /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/pool-bins/hclust/
cd /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/pool-bins/hclust/
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"

#########
echo ${mark}
#
state=0
Rscript ${scriptDir}do_hclust_for_pool-bins.R ${state} ${mark}
cp state${state}_${mark}_hlcust.pdf ${mark}_hlcust.pdf
#
for state in {1..43}
do
    echo ${state}
    Rscript ${scriptDir}do_hclust_for_pool-bins.R ${state} ${mark}
    pdfjam --outfile tmp.pdf ${mark}_hlcust.pdf state${state}_${mark}_hlcust.pdf
    mv tmp.pdf ${mark}_hlcust.pdf
done
