#!/bin/bash

# -- Kaili
# This script is for figuring out why samples with CTCF data are different.
# using dhs_ctcf _result
# state12 as an example for CTCF states

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

# Rscript make_figs_validate_CTCFstates_with_CTCFdata.R

# 1. signal of lung_0 in state 12
awk '{FS=OFS="\t"}{if($5=="12"){print $0}}' ./state_bed/lung_0_state.bed > ./ctcfstate_with_ctcfdata/lung0_state12.bed
if [ -f ./ctcfstate_with_ctcfdata/lung0_state12_signal.txt ];then rm ./ctcfstate_with_ctcfdata/lung0_state12_signal.txt; fi
#
while read line
do
    sample=`awk '{print $1}' <<< $line`
    if [ $sample == "lung_0" ];then
        mark=`awk '{print $2}' <<< $line`
        path0=`awk '{print $3}' <<< $line`
        path=${path0%.txt}
        echo $path
        #
        if [ $mark == "DNAme" ];then
            awk -v mark="$mark" '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$6,mark}}}' ./ctcfstate_with_ctcfdata/lung0_state12.bed ${path}.tab >> ./ctcfstate_with_ctcfdata/lung0_state12_signal.txt
        else
            awk -v mark="$mark" '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5,mark}}}' ./ctcfstate_with_ctcfdata/lung0_state12.bed ${path}.tab >> ./ctcfstate_with_ctcfdata/lung0_state12_signal.txt
        fi
    fi
done < dhs_ctcf.input
####
if [ -f ./ctcfstate_with_ctcfdata/lung0_allstates_signal.txt ];then rm ./ctcfstate_with_ctcfdata/lung0_allstates_signal.txt; fi
#
while read line
do
    sample=`awk '{print $1}' <<< $line`
    if [ $sample == "lung_0" ];then
        mark=`awk '{print $2}' <<< $line`
        path0=`awk '{print $3}' <<< $line`
        path=${path0%.txt}
        echo $path
        #
        if [ $mark == "DNAme" ];then
            awk -v mark="$mark" '{FS=OFS="\t"}{print $1,$6,mark}' ${path}.tab >> ./ctcfstate_with_ctcfdata/lung0_allstates_signal.txt
        else
            awk -v mark="$mark" '{FS=OFS="\t"}{print $1,$5,mark}' ${path}.tab >> ./ctcfstate_with_ctcfdata/lung0_allstates_signal.txt
        fi
    fi
done < dhs_ctcf.input


# 2. number of CTCF states in different samples
if [ -f ./ctcfstate_with_ctcfdata/sample_state12_count.txt ];then rm ./ctcfstate_with_ctcfdata/sample_state12_count.txt; fi
#
while read sample
do
    echo $sample
    #
    num=`awk '{FS=OFS="\t"}{if($5=="12"){print $0}}' ./state_bed/${sample}_state.bed | wc -l | awk '{print $1}'`
    echo -e $sample"\t"$num >> ./ctcfstate_with_ctcfdata/sample_state12_count.txt
done < all_66_sample.txt

# 3. regression
## randomly pick 100,000 bins, using the signal of other 10 marks to interpret CTCF signals
random_bins="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_signal_based_random.bed"
## lung_0 as example
## 1) get signal
echo "id" > ./ctcfstate_with_ctcfdata/lung0_signal_random100k.txt
cut -f 4 ${random_bins} >> ./ctcfstate_with_ctcfdata/lung0_signal_random100k.txt
#
while read line
do
    sample=`awk '{print $1}' <<< $line`
    if [ $sample == "lung_0" ];then
        mark=`awk '{print $2}' <<< $line`
        path0=`awk '{print $3}' <<< $line`
        path=${path0%.txt}
        echo $path
        #
        if [ $mark == "DNAme" ];then
            echo -e "id\t"$mark > tmp.txt
            awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$6}else{print $4,a[$4]}}' ${path}.tab ${random_bins} >> tmp.txt
        else
            echo -e "id\t"$mark > tmp.txt
            awk  '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $4,a[$4]}}' ${path}.tab ${random_bins} >> tmp.txt
        fi
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$1]}}' tmp.txt ./ctcfstate_with_ctcfdata/lung0_signal_random100k.txt > tmp2.txt
        mv tmp2.txt ./ctcfstate_with_ctcfdata/lung0_signal_random100k.txt
    fi
done < dhs_ctcf.input
rm tmp.txt

# 4. ARI in HOX region
# chr2:74,644,287-74,860,586
echo -e "chr2\t74644287\t74860586" > ./ctcfstate_with_ctcfdata/text_hox_region.bed
intersectBed -a /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_signal_based.bed -b ./ctcfstate_with_ctcfdata/text_hox_region.bed -wa -u | cut -f 4 > ./ctcfstate_with_ctcfdata/hox_region.bed
#
head -1 ctcf_9sample_impute.state > ./ctcfstate_with_ctcfdata/hox_region.state
awk '{FS=OFS}{if(NR==FNR){a[$1]=1}else{if(a[$1]){print $0}}}' ./ctcfstate_with_ctcfdata/hox_region.bed ctcf_9sample_impute.state >> ./ctcfstate_with_ctcfdata/hox_region.state
#
# sbatch state_conservation_hox_region_ward.sh
## z018 8432624
pdfjam --outfile state_conservation_hox_region_state.pdf state_conservation_hox_region_state_0.pdf state_conservation_hox_region_state_1.pdf state_conservation_hox_region_state_2.pdf state_conservation_hox_region_state_3.pdf state_conservation_hox_region_state_4.pdf state_conservation_hox_region_state_5.pdf state_conservation_hox_region_state_6.pdf state_conservation_hox_region_state_7.pdf state_conservation_hox_region_state_8.pdf state_conservation_hox_region_state_9.pdf state_conservation_hox_region_state_10.pdf state_conservation_hox_region_state_11.pdf state_conservation_hox_region_state_12.pdf state_conservation_hox_region_state_13.pdf state_conservation_hox_region_state_14.pdf state_conservation_hox_region_state_15.pdf state_conservation_hox_region_state_16.pdf state_conservation_hox_region_state_17.pdf state_conservation_hox_region_state_18.pdf state_conservation_hox_region_state_19.pdf state_conservation_hox_region_state_20.pdf state_conservation_hox_region_state_21.pdf state_conservation_hox_region_state_22.pdf state_conservation_hox_region_state_23.pdf state_conservation_hox_region_state_24.pdf state_conservation_hox_region_state_25.pdf state_conservation_hox_region_state_26.pdf state_conservation_hox_region_state_27.pdf state_conservation_hox_region_state_28.pdf state_conservation_hox_region_state_29.pdf state_conservation_hox_region_state_30.pdf state_conservation_hox_region_state_31.pdf state_conservation_hox_region_state_32.pdf state_conservation_hox_region_state_33.pdf state_conservation_hox_region_state_34.pdf state_conservation_hox_region_state_35.pdf state_conservation_hox_region_state_36.pdf state_conservation_hox_region_state_37.pdf state_conservation_hox_region_state_38.pdf state_conservation_hox_region_state_39.pdf state_conservation_hox_region_state_40.pdf state_conservation_hox_region_state_41.pdf state_conservation_hox_region_state_42.pdf state_conservation_hox_region_state_43.pdf state_conservation_hox_region_state_44.pdf state_conservation_hox_region_state_45.pdf state_conservation_hox_region_state_46.pdf

Rscript /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/make_state_conservation_heatmap_ARI_ward_2.R /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/ari_ward/ /data/zusers/fankaili/ideas/dhs_ctcf/ctcfstate_with_ctcfdata/hox_region.state 100 66 hox_region
