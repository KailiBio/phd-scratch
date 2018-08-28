#!/bin/bash

# -- Kaili
# This script is for running t-SNE for all the 21 biosamples.
# !!! all the ubi-rDHS files in the same orders!!!

masterFile="/data/zusers/moorej3/ENCODE-Registry/hg19/V4/Cell-Type-Specific/Master-Cell-List.txt"
zscoreDir="/data/zusers/moorej3/ENCODE-Registry/hg19/V4/signal-output/"
cell_type_specificDir="/data/zusers/moorej3/ENCODE-Registry/hg19/V4/Cell-Type-Specific/Five-Group/"
scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/bimodal_TF_in_ubi-rDHS/scipts/"

workDir="/data/zusers/fankaili/ccre/tf/tsne/"

#####################
# for 10,921 ubi-rDHS
#####################


# 1. get the matrix
cd ${workDir}

while read line
do
    echo ${line} ;
    dnaseExp=`awk '{FS=OFS="\t"}{print $1}' <<< ${line}` ;
    dnaseID=`awk '{FS=OFS="\t"}{print $2}' <<< ${line}` ;
    h3k4me3Exp=`awk '{FS=OFS="\t"}{print $3}' <<< ${line}` ;
    h3k4me3ID=`awk '{FS=OFS="\t"}{print $4}' <<< ${line}` ;
    h3k27acExp=`awk '{FS=OFS="\t"}{print $5}' <<< ${line}` ;
    h3k27acID=`awk '{FS=OFS="\t"}{print $6}' <<< ${line}` ;
    ctcfExp=`awk '{FS=OFS="\t"}{print $7}' <<< ${line}` ;
    ctcfID=`awk '{FS=OFS="\t"}{print $8}' <<< ${line}` ;
    cellline=`awk '{FS=OFS="\t"}{print $9}' <<< ${line}` ;
    # log10(distance)
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$4]=$5}else{if(a[$1]){print $1,b[$1]}}}' hg19_ubi-rDHS_distance2TSS_log10.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt > temp.distance.txt ;
    # DNase
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${dnaseExp}-${dnaseID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt > temp.dnase.txt ;
    # H3K4me3
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${h3k4me3Exp}-${h3k4me3ID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt > temp.h3k4me3.txt ;
    # H3K27ac
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${h3k27acExp}-${h3k27acID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt > temp.h3k27ac.txt ;
    # CTCF
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${ctcfExp}-${ctcfID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list.txt > temp.ctcf.txt ;
    #
    paste temp.distance.txt temp.dnase.txt temp.h3k4me3.txt temp.h3k27ac.txt temp.ctcf.txt > temp.matrix.txt ;
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$10}else{if(a[$1]){print $0,a[$1]}}}' \
    ${cell_type_specificDir}${dnaseID}_${h3k4me3ID}_${h3k27acID}_${ctcfID}.cREs.bed temp.matrix.txt > ./matrix/hg19_ubi-rDHS_${cellline}_five_dimension_matrix.txt
    sed -i '1iid\tdistance\tdnase\th3k4me3\th3k27ac\tctcf\ttype' ./matrix/hg19_ubi-rDHS_${cellline}_five_dimension_matrix.txt
    rm temp.*.txt ;
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_master.txt

# 2. get agonotic definition
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{if(a[$1]){print $1,a[$1]}}}' /data/zusers/fankaili/ccre/hg19-cREs.bed \
/data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt > /data/zusers/fankaili/ccre/temp.hg19_ubi_ccRE_ccreid_agnostic.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$3}else{print $0,a[$2]}}' /data/zusers/fankaili/ccre/color_list.txt \
/data/zusers/fankaili/ccre/temp.hg19_ubi_ccRE_ccreid_agnostic.txt > /data/zusers/fankaili/ccre/hg19_ubi_ccRE_ccreid_agnostic.txt
#
rm /data/zusers/fankaili/ccre/temp.hg19_ubi_ccRE_ccreid_agnostic.txt

# 3. get DNAme
## make cell-line DNAme list: /data/zusers/fankaili/ccre/tf/hg19_four_dimension_cellline_DNAme.txt
cd ${workDir}

echo "id" > hg19_ubi-rDHS_DNAme_matrix.txt
cat /data/zusers/fankaili/ccre/ubi_ccREs_hg19_list_ccreid.txt >> hg19_ubi-rDHS_DNAme_matrix.txt
while read line
do
    cellline=`awk '{FS=OFS="\t"}{print $1}' <<< ${line}` ;
    id=`awk '{FS=OFS="\t"}{print $2}' <<< ${line}` ;
    if [[ ${id} != *"NA"* ]];then
        echo ${cellline} ;
        echo ${cellline} > temp.txt ;
        n=`awk -v id="$id" '{FS=OFS="\t"}{if(NR==1){for(i=1;i<=NF;i++){if($i==id){print i}}}}' /data/zusers/fankaili/ccre_old/hg19_ccREs_DNAme.txt` ;
        echo $n ;
        awk -v n="$n" '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$4]=$n}else{if(a[$1]){print b[$1]}}}' \
        /data/zusers/fankaili/ccre_old/hg19_ccREs_DNAme.txt hg19_ubi-rDHS_DNAme_matrix.txt >> temp.txt ;
        paste hg19_ubi-rDHS_DNAme_matrix.txt temp.txt > temp2.txt ;
        mv temp2.txt hg19_ubi-rDHS_DNAme_matrix.txt ;
        rm temp.txt ;
    fi
done < /data/zusers/fankaili/ccre/tf/hg19_four_dimension_cellline_DNAme.txt

# 4. get all 21 cell-type list
awk '{FS=OFS="\t"}{print $9}' /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_master.txt > /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt

# 5. run t-SNR and make figures
while read line
do
    echo ${line};
    mkdir /data/zusers/fankaili/ccre/tf/tsne/matrix/${line}
    Rscript ${scriptDir}run_tSNE.R /data/public_html_users/fankaili/tSNE_10921/ /data/zusers/fankaili/ccre/tf/tsne/matrix/ \
    /data/zusers/fankaili/ccre/hg19_ubi_ccRE_ccreid_agnostic.txt /data/zusers/fankaili/ccre/tf/tsne/hg19_ubi-rDHS_DNAme_matrix.txt ${line}
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt

# outDir="/data/public_html_users/fankaili/tSNE_10921/"
# matrixDir="/data/zusers/fankaili/ccre/tf/tsne/matrix/"
# agnosticFile="/data/zusers/fankaili/ccre/hg19_ubi_ccRE_ccreid_agnostic.txt"
# dnameFile="/data/zusers/fankaili/ccre/tf/tsne/hg19_ubi-rDHS_DNAme_matrix.txt"
# celltype="2_GM12878"

# 6. show figures using html
cd /data/public_html_users/fankaili/tSNE_10921/

# link for one biosample
while read line
do
    echo ${line} ;
    head -7 /data/public_html_users/fankaili/ccREs_TF/index.html > ./${line}/index.html
    # echo '<div style="text-align:center;font-size:20px;color:red;">tSNE for '${line}'</div>' >> ./${line}/index.html ;
    echo '<div style="text-align:center"><embed width="1600" height="1000" src="hg19_ubi-rDHS_'${line}'_tSNE_in_'${line}'.pdf"> </embed></div>' >> ./${line}/index.html ;
    while read ll
    do
        if [[ "$ll" != "$line" ]]; then
            echo '<div style="text-align:center"><embed width="1600" height="1000" src="hg19_ubi-rDHS_'${line}'_tSNE_in_'${ll}'.pdf"> </embed></div>' >> ./${line}/index.html ;
        fi
    done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt
    tail -3 /data/public_html_users/fankaili/ccREs_TF/index.html >> ./${line}/index.html
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt

# link for all
rm index.html
while read line
do
    echo ${line} ;
    echo '<div style="text-align:left;font-size:20px;color:red;">'${line}'</div>' >> index.html ;
    echo '<p><a href="./'${line}'/index.html">tSNE in '${line}'</a ></p >' >> index.html ;
    i=0
    while read ll
    do
        a=$((i+1));
        i=${a} ;
        echo '<p><a href="./'${line}'/hg19_ubi-rDHS_'${line}'_tSNE_in_'${ll}'.pdf">'${i}'. '${ll}'</a ></p >' >> index.html ;
    done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt





#####################
# for 29,773 ubi-rDHS
#####################


# 1. get the matrix
cd ${workDir}

while read line
do
    echo ${line} ;
    dnaseExp=`awk '{FS=OFS="\t"}{print $1}' <<< ${line}` ;
    dnaseID=`awk '{FS=OFS="\t"}{print $2}' <<< ${line}` ;
    h3k4me3Exp=`awk '{FS=OFS="\t"}{print $3}' <<< ${line}` ;
    h3k4me3ID=`awk '{FS=OFS="\t"}{print $4}' <<< ${line}` ;
    h3k27acExp=`awk '{FS=OFS="\t"}{print $5}' <<< ${line}` ;
    h3k27acID=`awk '{FS=OFS="\t"}{print $6}' <<< ${line}` ;
    ctcfExp=`awk '{FS=OFS="\t"}{print $7}' <<< ${line}` ;
    ctcfID=`awk '{FS=OFS="\t"}{print $8}' <<< ${line}` ;
    cellline=`awk '{FS=OFS="\t"}{print $9}' <<< ${line}` ;
    # log10(distance)
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$4]=$5}else{if(a[$1]){print $1,b[$1]}}}' hg19_ubi-rDHS_2_distance2TSS_log10.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list_ccreid.txt > temp.distance.txt ;
    # DNase
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${dnaseExp}-${dnaseID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list.txt > temp.dnase.txt ;
    # H3K4me3
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${h3k4me3Exp}-${h3k4me3ID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list.txt > temp.h3k4me3.txt ;
    # H3K27ac
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${h3k27acExp}-${h3k27acID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list.txt > temp.h3k27ac.txt ;
    # CTCF
    awk '{FS=OF"\t"}{if(NR==FNR){a[$1]=1;b[$1]=$2}else{if(a[$1]){print b[$1]}}}' ${zscoreDir}${ctcfExp}-${ctcfID}.txt \
    /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list.txt > temp.ctcf.txt ;
    #
    paste temp.distance.txt temp.dnase.txt temp.h3k4me3.txt temp.h3k27ac.txt temp.ctcf.txt > temp.matrix.txt ;
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$10}else{if(a[$1]){print $0,a[$1]}}}' \
    ${cell_type_specificDir}${dnaseID}_${h3k4me3ID}_${h3k27acID}_${ctcfID}.cREs.bed temp.matrix.txt > ./matrix_2/hg19_ubi-rDHS_2_${cellline}_five_dimension_matrix.txt
    sed -i '1iid\tdistance\tdnase\th3k4me3\th3k27ac\tctcf\ttype' ./matrix_2/hg19_ubi-rDHS_2_${cellline}_five_dimension_matrix.txt
    rm temp.*.txt ;
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_master.txt

# 2. get agonotic definition
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$9}else{if(a[$1]){print $1,a[$1]}}}' /data/zusers/fankaili/ccre/hg19-cREs.bed \
/data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list_ccreid.txt > /data/zusers/fankaili/ccre/temp.hg19_ubi_ccRE_2_ccreid_agnostic.txt
#
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$3}else{print $0,a[$2]}}' /data/zusers/fankaili/ccre/color_list.txt \
/data/zusers/fankaili/ccre/temp.hg19_ubi_ccRE_2_ccreid_agnostic.txt > /data/zusers/fankaili/ccre/hg19_ubi_ccRE_2_ccreid_agnostic.txt
#
rm /data/zusers/fankaili/ccre/temp.hg19_ubi_ccRE_2_ccreid_agnostic.txt

# 3. get DNAme
## make cell-line DNAme list: /data/zusers/fankaili/ccre/tf/hg19_four_dimension_cellline_DNAme.txt
cd ${workDir}

echo "id" > hg19_ubi-rDHS_2_DNAme_matrix.txt
cat /data/zusers/fankaili/ccre/ubi_ccREs_2_hg19_list_ccreid.txt >> hg19_ubi-rDHS_2_DNAme_matrix.txt
while read line
do
    cellline=`awk '{FS=OFS="\t"}{print $1}' <<< ${line}` ;
    id=`awk '{FS=OFS="\t"}{print $2}' <<< ${line}` ;
    if [[ ${id} != *"NA"* ]];then
        echo ${cellline} ;
        echo ${cellline} > temp.txt ;
        n=`awk -v id="$id" '{FS=OFS="\t"}{if(NR==1){for(i=1;i<=NF;i++){if($i==id){print i}}}}' /data/zusers/fankaili/ccre_old/hg19_ccREs_DNAme.txt` ;
        echo $n ;
        awk -v n="$n" '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1;b[$4]=$n}else{if(a[$1]){print b[$1]}}}' \
        /data/zusers/fankaili/ccre_old/hg19_ccREs_DNAme.txt hg19_ubi-rDHS_2_DNAme_matrix.txt >> temp.txt ;
        paste hg19_ubi-rDHS_2_DNAme_matrix.txt temp.txt > temp2.txt ;
        mv temp2.txt hg19_ubi-rDHS_2_DNAme_matrix.txt ;
        rm temp.txt ;
    fi
done < /data/zusers/fankaili/ccre/tf/hg19_four_dimension_cellline_DNAme.txt


# 5. run t-SNR and make figures
while read line
do
    echo ${line};
    mkdir /data/public_html_users/fankaili/tSNE_29733/${line}
    Rscript ${scriptDir}run_tSNE.R /data/public_html_users/fankaili/tSNE_29733/ /data/zusers/fankaili/ccre/tf/tsne/matrix_2/ \
    /data/zusers/fankaili/ccre/hg19_ubi_ccRE_2_ccreid_agnostic.txt /data/zusers/fankaili/ccre/tf/tsne/hg19_ubi-rDHS_2_DNAme_matrix.txt 2_${line}
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt

# outDir="/data/public_html_users/fankaili/tSNE_29733/"
# matrixDir="/data/zusers/fankaili/ccre/tf/tsne/matrix_2/"
# agnosticFile="/data/zusers/fankaili/ccre/hg19_ubi_ccRE_2_ccreid_agnostic.txt"
# dnameFile="/data/zusers/fankaili/ccre/tf/tsne/hg19_ubi-rDHS_2_DNAme_matrix.txt"
# celltype0="2_GM12878"

# 6. show figures using html

cd /data/public_html_users/fankaili/tSNE_29733/

# link for one biosample
while read line
do
    echo ${line} ;
    head -7 /data/public_html_users/fankaili/ccREs_TF/index.html > ./${line}/index.html
    # echo '<div style="text-align:center;font-size:20px;color:red;">tSNE for '${line}'</div>' >> ./${line}/index.html ;
    echo '<div style="text-align:center"><embed width="1600" height="1000" src="hg19_ubi-rDHS_'${line}'_tSNE_in_'${line}'.pdf"> </embed></div>' >> ./${line}/index.html ;
    while read ll
    do
        if [[ "$ll" != "$line" ]]; then
            echo '<div style="text-align:center"><embed width="1600" height="1000" src="hg19_ubi-rDHS_'${line}'_tSNE_in_'${ll}'.pdf"> </embed></div>' >> ./${line}/index.html ;
        fi
    done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt
    tail -3 /data/public_html_users/fankaili/ccREs_TF/index.html >> ./${line}/index.html
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt

# link for all
rm index.html
while read line
do
    echo ${line} ;
    echo '<div style="text-align:left;font-size:20px;color:red;">'${line}'</div>' >> index.html ;
    echo '<p><a href="./'${line}'/index.html">tSNE in '${line}'</a ></p >' >> index.html ;
    i=0
    while read ll
    do
        a=$((i+1));
        i=${a} ;
        echo '<p><a href="./'${line}'/hg19_ubi-rDHS_'${line}'_tSNE_in_'${ll}'.pdf">'${i}'. '${ll}'</a ></p >' >> index.html ;
    done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt
done < /data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt
