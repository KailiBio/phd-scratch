#!/bin/bash

# -- Kaili
# This script is for calculating zscore for all biosamples in given list.
# INPUT:
# OUROUT:

#########################
help_info(){
	echo "usage:"
	echo "sh calculate_zscore.sh <option>* [-l list] [-d ourDir] [-t tf]"
	echo ""
	echo "This file is for makding TF zscore matrix of all biosamples."
	echo ""
	echo "Arguments:"
	echo "-l list of celllines for zscore calculating."
	echo "-d Dir of output file."
	echo "-t which tf: CTCF, SMC3, RAD21 et al."
	echo ""
	echo "Any questions, please contact me."
	echo "		--Kaili Fan (fankaili.bio@gmail.com)"
	echo ""
}

if [ $# -lt 5 ];then
	help_info
	exit 1
fi

while getopts "l:d:t:" Arg
do
	case $Arg in
		l)	list=$OPTARG;;
		d)	outDir=$OPTARG;;
		t)	tf=$OPTARG;;
		?)	echo "Wrong parameter!!!"
			exit 1;;
	esac
done

#########################

cd ${outDir}
#
for cellline in `cat ${list}`
do
    echo ${cellline} ;
    line=`grep ${tf} /data/zusers/fankaili/ccre/tf/encode_tf_file_list/encode_hg19_${cellline}_tf_id_list.txt` ;
    if [ "$line" != "" ]; then
        n=`echo $line | awk '{print NF}'` ;
        nn=`expr $(($n/3))` ;
        for ((i=1;i<=$nn;i++))
        do
            a=`expr $((1+($i-1)*3))`;
            b=`expr $((2+($i-1)*3))`;
            id=`awk -v i="$a" '{print $i}' <<< $line`;
            file_id=`awk -v i="$b" '{print $i}' <<< $line`;
            bigWigAverageOverBed /data/projects/encode/data/${id}/${file_id}.bigWig \
            /data/zusers/moorej3/ENCODE-Registry/hg19/V4/hg19-rDHSs.bed ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_${tf}_signal.tab ;
            python /data/zusers/fankaili/ccre/tf/zscore-normalization.py ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_${tf}_signal.tab > ${outDir}hg19_rDHS_${cellline}_${id}_${file_id}_${tf}_signal_zscore.txt ;
        done
    fi
done
