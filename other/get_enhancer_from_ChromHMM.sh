
cd /data/zusers/fankaili/shared/for_qin/
for file in `ls /data/public_html_users/vanderva/trackhub/chromhmmpaper/18states_10marks/*.bigBed`
do
    filename0=${file%.bigBed} ;
    filename=${filename0#/data/public_html_users/vanderva/trackhub/chromhmmpaper/18states_10marks/};
    echo ${filename};
    #bigBedToBed ${file} ./ChromHMM_18status_10mark_bed/${filename}.bed;
    awk '{FS=OFS="\t"}{if($4=="Enh"){print $0}}' ./ChromHMM_18status_10mark_bed/${filename}.bed > ./enhancer/${filename}_enhancer.bed ;
    awk '{FS=OFS="\t"}{if($4=="Tss" || $4=="TssFlnk"){print $0}}' ./ChromHMM_18status_10mark_bed/${filename}.bed > ./promoter/${filename}_promoter.bed ;
done
