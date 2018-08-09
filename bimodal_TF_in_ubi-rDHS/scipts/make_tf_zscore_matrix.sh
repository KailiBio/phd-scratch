#!/bin/bash

# -- Kaili
# This script is for making TF zscore matrix.
# INPUT: 1) rDHS list:
#        2) ubi-rDHS list:
#        3) zscore file folder:
# OUTPUT: matrix of all zscore of TF in all biosamples.

#########################
help_info(){
	echo "usage:"
	echo "sh make_tf_zscore_matrix.sh <option>* [-l list] [-d ourDir] [-m outMatrix] [-t tf] [-f zscore_folder]"
	echo ""
	echo "This file is for makding TF zscore matrix of all biosamples."
	echo ""
	echo "Arguments:"
	echo "-l bed file contails the list of ID."
	echo "-d Dir of output file."
	echo "-m The name of the output matrix."
	echo "-t which tf: CTCF, SMC3, RAD21 et al."
	echo "-f Path of folder contains all the zscore file from zscore-normalization.py."
	echo ""
	echo "Any questions, please contact me."
	echo "		--Kaili Fan (fankaili.bio@gmail.com)"
	echo ""
}

if [ $# -lt 5 ];then
	help_info
	exit 1
fi

while getopts "l:d:m:t:f:" Arg
do
	case $Arg in
		l)	list=$OPTARG;;
		d)	outDir=$OPTARG;;
		m)	outMatrix=$OPTARG;;
		t)	tf=$OPTARG;;
		f)	zscore_folder=$OPTARG;;
		?)	echo "Wrong parameter!!!"
			exit 1;;
	esac
done

#########################
cd ${outDir}

echo "id" > ${outMatrix}
awk '{print $4}' ${list} >> ${outMatrix}
#
for file in `ls ${zscore_folder}*_${tf}_signal_zscore.txt`
do
    filename0=${file%_${tf}_signal_zscore.txt} ;
    filename=${filename0#${zscore_folder}hg19_rDHS_} ;
    echo $filename ;
    echo -e "id\t"$filename > temp.txt ;
    awk '{FS=" \t";OFS="\t"}{print $1,$2}' $file >> temp.txt ;
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print $0,b[$1]}}}' temp.txt ${outMatrix} > temp_matrix.txt ;
    mv temp_matrix.txt ${outMatrix} ;
    rm temp*.txt;
done
