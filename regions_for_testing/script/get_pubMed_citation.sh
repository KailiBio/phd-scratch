#!/bin/bash

# -- Kaili
# This script is for getting gene pubMed citation number as long with link.

if  [ $# -lt 1 ];then
	echo -e "$0 your.gene.table(ensemble_id in the 2nd column and gene_name in the 3rd column) out.file"
	exit 1
fi
#
cat /home/fankaili/genome/gene2ensembl | awk 'BEGIN{FS=OFS="\t"} {if(NR==FNR){a[$3]=$2}else{split($1,b,".");print a[b[1]]?a[b[1]]:"NA",$2,$1}}' - $1 > /data/zusers/fankaili/regions/$3/tmp.gene.list
python /data/zusers/fankaili/github/weng-lab/Kaili/regions_for_testing/script/getGeneRefNumber.py /data/zusers/fankaili/regions/$3/tmp.gene.list > $2
rm /data/zusers/fankaili/regions/$3/tmp.gene.list
