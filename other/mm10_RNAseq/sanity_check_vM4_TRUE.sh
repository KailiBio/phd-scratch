
# -- Kaili
# This script is for sanity check of re-runing DEcall in vM4.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/"
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"

cd ${workDir}

#############

sanity_check(){
    outfile=$1
    DEG_folder=$2
    #
    if [ -f ${outfile} ];then rm ${outfile}; fi
    #
    for file in `ls ./${DEG_folder}/`
    do
        sample1=`awk '{split($0,a,"_VS_");print a[1]}' <<< ${file}`
        sample2=`awk '{split($1,a,"_VS_");split(a[2],b,".txt");print b[1]}' <<< ${file}`
        echo ${sample1}" VS "${sample2}
        # get match sample names
        if [[ $sample1 == "facial"* ]]; then
            sample11=${sample1/facial/embryonic_facial_prominence}
        elif [[ $sample1 == "neural.tube"* ]];then
            sample11=${sample1/neural.tube/neural_tube}
        else
            sample11=$sample1
        fi
        #
        if [[ $sample2 == *"facial"* ]]; then
            sample22=${sample2/facial/embryonic_facial_prominence}
        elif [[ $sample2 == "neural.tube"* ]];then
            sample22=${sample2/neural.tube/neural_tube}
        else
            sample22=$sample2
        fi
        # for Junko's
        if [ -f /data/projects/screen/Version-4/mouse_epigenome/de_all_pairs/data/${sample11}_VS_${sample22}.txt.gz ];then
            cp /data/projects/screen/Version-4/mouse_epigenome/de_all_pairs/data/${sample11}_VS_${sample22}.txt.gz ./
            gzip -d ${sample11}_VS_${sample22}.txt.gz
            #
            echo -e "id\tlog2FoldChange\tpadj" > tmp.matrix.txt
            awk '{FS=OFS="\t"}{if(NR>1){print $1,$3,$7}}' ${sample11}_VS_${sample22}.txt >> tmp.matrix.txt
            awk '{if(NR>1 && $7<0.05 && ($3>2 || $3<-2)){print $1}}' ${sample11}_VS_${sample22}.txt > tmp1.significant.txt
            #
            rm ${sample11}_VS_${sample22}.txt
        elif [ -f /data/projects/screen/Version-4/mouse_epigenome/de_all_pairs/data/${sample22}_VS_${sample11}.txt.gz ];then
            cp /data/projects/screen/Version-4/mouse_epigenome/de_all_pairs/data/${sample22}_VS_${sample11}.txt.gz ./
            gzip -d ${sample22}_VS_${sample11}.txt.gz
            #
            echo -e "id\tlog2FoldChange\tpadj" > tmp.matrix.txt
            awk '{FS=OFS="\t"}{if(NR>1){print $1,$3,$7}}' ${sample22}_VS_${sample11}.txt >> tmp.matrix.txt
            awk '{if(NR>1 && $7<0.05 && ($3>2 || $3<-2)){print $1}}' ${sample22}_VS_${sample11}.txt > tmp1.significant.txt
            #
            rm ${sample22}_VS_${sample11}.txt
        fi
        #
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$3;b[$1]=$7}else{if(FNR==1){print $0,"my_log2FoldChange","my_padj"}else{print $0,a[$1],b[$1]}}}' ./${DEG_folder}/${file} tmp.matrix.txt > tmp.matrix2.txt
        awk '{if(NR>1 && $7<0.05 && ($3>2 || $3<-2)){print $1}}' ./${DEG_folder}/${file} > tmp2.significant.txt
        #
        cat tmp1.significant.txt tmp2.significant.txt | sort | uniq -d > tmp.co_sig.txt
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(FNR==1){print $0,"sig"}else{if(a[$1]){print $0,"co_sig"}else{print $0,"non_sig"}}}}' tmp.co_sig.txt tmp.matrix2.txt > tmp.matrix3.txt
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(FNR==1){print $0}else{if(a[$1] && $6=="co_sig"){print $0}else if(a[$1] && $6=="non_sig"){print $1,$2,$3,$4,$5,"sig_Junko_only"}else{print $0}}}}' tmp1.significant.txt tmp.matrix3.txt > tmp.matrix4.txt
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(FNR==1){print $0}else{if(a[$1] && $6=="co_sig"){print $0}else if(a[$1] && $6=="non_sig"){print $1,$2,$3,$4,$5,"sig_mine_only"}else{print $0}}}}' tmp2.significant.txt tmp.matrix4.txt > tmp.matrix5.txt
        mv tmp.matrix5.txt ./sanity_check_TRUE/${sample1}_VS_${sample2}_matrix.txt
        #
        Rscript ${scriptDir}make_DE_scatter_volcano.R ${sample1}"_VS_"${sample2}
    done
}

outfile="DE_call_comparison_M4_read-counts.txt"
DEG_folder="de_M4_read_TRUE"
#
sanity_check ${outfile} ${DEG_folder}
