#!/bin/bash/

# -- Kaili
# This script is for running IDEAS for MACS2 FE file.

cd /data/zusers/fankaili/ideas/run_ideas_fe/

cp /data/zusers/fankaili/ideas/run_ideas/IDEAS_8hm_atac_dname_2* ./
cp /data/zusers/fankaili/ideas/run_ideas/mm10_200bin.bed ./

##################
mv IDEAS_8hm_atac_dname_2.input run_ideas_fe_1.input
sed -i 's/\/data\/zusers\/fankaili\/ideas\/run_ideas\//\/data\/zusers\/fankaili\/ideas\/run_ideas_fe\//g' run_ideas_fe_1.input
mv IDEAS_8hm_atac_dname_2.sh run_ideas_fe_1.sh
vim run_ideas_fe_1.sh
mv IDEAS_8hm_atac_dname_2.parafile run_ideas_fe_1.parafile
vim run_ideas_fe_1.parafile

# cp /data/zusers/fankaili/ideas/run_ideas/signal/*_200bin.tab ./signal/
for file in `ls /data/zusers/fankaili/ideas/run_ideas/signal/`
do
  if [[ "${file: -11}" == "_200bin.tab" ]];then
    echo ${file};
    filename=${file%.tab}
    echo ${filename}
    if [[ "${file: -17}" == "_DNAme_200bin.tab" ]];then
        echo "DNAme";
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print a[$4]}}' /data/zusers/fankaili/ideas/run_ideas/signal/${file} mm10_200bin.bed > ./signal/${filename}.txt;
    else
        echo "ATAC or ChIP";
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print a[$4]}}' /data/zusers/fankaili/ideas/run_ideas/signal/${file} mm10_200bin.bed > ./signal/${filename}.txt;
    fi
  fi
done

nohup bash run_ideas_fe_1.sh > nohup.run_ideas_fe_1.out 2>&1&

##################
cp run_ideas_fe_1.input run_ideas_fe_2.input
cp run_ideas_fe_1.parafile run_ideas_fe_2.parafile
vim run_ideas_fe_2.parafile
cp run_ideas_fe_1.sh run_ideas_fe_2.sh
vim run_ideas_fe_2.sh

nohup bash run_ideas_fe_2.sh > nohup.run_ideas_fe_2.out 2>&1&
