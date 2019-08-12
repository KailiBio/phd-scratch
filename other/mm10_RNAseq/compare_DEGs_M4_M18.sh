
# -- Kaili
# This script is for sanity check of re-runing DEcall in vM4.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/"
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"

cd ${workDir}

#############

sanity_check(){
    DEG_folder=$1
    outfile=$2
    #
    echo -e "sample\tM4\tM18\toverlap" > ${outfile}
    #
    for file in `ls ./${DEG_folder}/`
    do
        sample1=`awk '{split($0,a,"_VS_");print a[1]}' <<< ${file}`
        sample2=`awk '{split($1,a,"_VS_");split(a[2],b,".txt");print b[1]}' <<< ${file}`
        echo ${sample1}" VS "${sample2}
        # M4
        echo -e "id\tM4\tlog2FoldChange\tpadj\tM18" > tmp.matrix.txt
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$3;b[$1]=$7}else{if(FNR>1){print $1,$2,a[$2],b[$2],$3}}}' ./de_M4_read_TRUE/${file} mm10_M4_M18_overlap_genelist2.txt >> tmp.matrix.txt
        # M18
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$3;b[$1]=$7}else{if(FNR==1){print $0,"my_log2FoldChange","my_padj"}else{print $0,a[$5],b[$5]}}}' ./${DEG_folder}/${file} tmp.matrix.txt > tmp.matrix2.txt
        #### significant
        cut -f 1,3,4,6,7 tmp.matrix2.txt > tmp.matrix3.txt
        awk '{if(NR>1 && $3<0.05 && ($2>2 || $2<-2)){print $1}}'  tmp.matrix3.txt > tmp.M4_significant.txt
        awk '{if(NR>1 && $5<0.05 && ($4>2 || $4<-2)){print $1}}'  tmp.matrix3.txt > tmp.M18_significant.txt
        #
        cat tmp.M4_significant.txt tmp.M18_significant.txt | sort | uniq -d > tmp.co_sig.txt
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(FNR==1){print $0,"sig"}else{if(a[$1]){print $0,"co_sig"}else{print $0,"non_sig"}}}}' tmp.co_sig.txt tmp.matrix3.txt > tmp.matrix4.txt
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(FNR==1){print $0}else{if(a[$1] && $6=="co_sig"){print $0}else if(a[$1] && $6=="non_sig"){print $1,$2,$3,$4,$5,"sig_M4_only"}else{print $0}}}}' tmp.M4_significant.txt tmp.matrix4.txt > tmp.matrix5.txt
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=1}else{if(FNR==1){print $0}else{if(a[$1] && $6=="co_sig"){print $0}else if(a[$1] && $6=="non_sig"){print $1,$2,$3,$4,$5,"sig_M18_only"}else{print $0}}}}' tmp.M18_significant.txt tmp.matrix5.txt > tmp.matrix6.txt
        mv tmp.matrix6.txt ./compare_DEGs_M4_M18/${sample1}_VS_${sample2}_matrix.txt
        #
        Rscript ${scriptDir}make_DE_scatter_volcano.R ${sample1}"_VS_"${sample2}
        #
        a=`wc -l tmp.M4_significant.txt | awk '{print $1}'`
        b=`wc -l tmp.M18_significant.txt | awk '{print $1}'`
        c=`wc -l tmp.co_sig.txt | awk '{print $1}'`
        echo -e ${sample1}"_VS_"${sample2}"\t"${a}"\t"${b}"\t"${c} >> ${outfile}
    done
}

DEG_folder="de_M18_readCounts"
outfile="/data/zusers/fankaili/ccre/mm10_rnaseq/compare_DEGs_M4_M18/compare_DEGs_M4_M18_matrix.txt"
#
sanity_check ${DEG_folder} ${outfile}


rm tmp*.txt
