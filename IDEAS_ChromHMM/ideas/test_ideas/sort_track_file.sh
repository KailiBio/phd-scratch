#!/bin/bash

# -- Kaili
# This script is for sorting the track file by short label.
# EXP: bash sort_track_file.sh "/data/zusers/fankaili/ideas/yu_input/IDEAS_yu_input_2_result/Tracks/" "trackDb_IDEAS_yu_input_2"

# workDir="/data/zusers/fankaili/ideas/yu_input/IDEAS_yu_input_2_result/Tracks/"
# file="trackDb_IDEAS_yu_input_2"

#########################
help_info(){
	echo "usage:"
	echo "bash sort_track_file.sh -d workDir -f file"
	echo ""
	echo "Arguments:"
	echo "-d working directory."
	echo "-f file prefix."
	echo ""
	echo "Any questions, please contact me."
	echo "		--Kaili Fan (fankaili.bio@gmail.com)"
	echo ""
}

if [ $# -lt 2 ];then
	help_info
	exit 1
fi

while getopts "d:f:" Arg
do
	case $Arg in
		d)	workDir=$OPTARG;;
		f)	file=$OPTARG;;
		?)	echo "Wrong parameter!!!"
			exit 1;;
	esac
done


########
cd ${workDir}

cp ${file}.txt ${file}_old.txt

awk '{FS=" ";OFS="\t"}{if($1=="priority"){printf $2}else{if($1=="shortLabel"){printf "\t"$2"\n"}}}' ${file}_old.txt | sort -k2 > track_order.txt

awk -v RS= '{print > ("whatever-" NR ".txt")}' ${file}_old.txt

n=0
for i in `cut -f 1 track_order.txt`
do
    echo ${i};
    n=$(($n+1));
    echo ${n};
    sed -i "s/priority\ $i/priority\ $n/g" whatever-${i}.txt;
    cat whatever-${i}.txt >> ${file}.txt;
    echo "" >> ${file}.txt;
done

rm whatever-*.txt
