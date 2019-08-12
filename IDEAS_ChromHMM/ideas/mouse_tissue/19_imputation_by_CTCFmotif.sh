#!/bin/bash

# -- Kaili
# This is script is for analyzing CTCF imputation by CTCF motif.
# 1. randomly pick 100k bins, do regression for predicted power
# 2. Will parafile fold influence prediction?
# 3. how about imputating new samples? 8impute3
# 4. compare with other imputation methods
# 5. add DNAme signal track

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/dhs_ctcf/"

cd ${workDir}

# 1. randomly pick 100k bins, do regression for predicted power
random_bins="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/mm10_OCR-center_bins_v3_signal_based_random.bed"
## 1) get signal
## lung_0
echo "id" > ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt
cut -f 4 ${random_bins} >> ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt
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
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$1]}}' tmp.txt ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt > tmp2.txt
        mv tmp2.txt ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt
    fi
done < dhs_ctcf.input
## CTCF motif
echo -e "id\tCTCFmotif" > tmp.txt
awk  '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$4]){print $4,a[$4]}else{print $4,0}}}' CTCFmotif_bin_withScore.txt ${random_bins} >> tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$1]}}' tmp.txt ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt > tmp2.txt
mv tmp2.txt ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt
## CTCF signal average
echo -e "id\tCTCFaverage" > tmp.txt
paste mm10_OCR-center_bins_v3_signal_based_space_bedformat.bed 9sample_averageCTCFsignal.txt > 9sample_averageCTCFsignal.tab
awk  '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$5}else{print $4,a[$4]}}' 9sample_averageCTCFsignal.tab ${random_bins} >> tmp.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$1]}}' tmp.txt ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt > tmp2.txt
mv tmp2.txt ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt
## CTCF peak summit
intersectBed -a ${random_bins} -b ./peaks_validation_9To11/lung_0_ctcf_peak_center_withSignal.bed -wa -u > tmp.txt
echo -e "id\tCTCFpeakSummit" > tmp2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(FNR>1){if(a[$1]){print $1,1}else{print $1,0}}}}' tmp.txt ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt >> tmp2.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{print $0,a[$1]}}' tmp2.txt ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt > tmp3.txt
mv tmp3.txt ./ctcfstate_with_ctcfdata/lung0_13signal_random100k.txt
#
rm tmp*.txt
## 2) do regression
# Rscript do_regression_CTCFprediction.R


# 2. Will parafile fold influence prediction? run 9impute66
bash ${scriptDir}choose_parafile_fold.sh
bash ${scriptDir}impute_all66samples.sh


# 3. how about imputating new samples? 8impute3
bash ${scriptDir}impute_CTCFstates_for_new_tissues.sh


# 4. compare with other imputation methods
bash ${scriptDir}compare_imputation_tools.sh
bash ${scriptDir}compare_imputation_tools2.sh
bash ${scriptDir}compare_imputation_tools_wholeGenome.sh


# 5. add DNAme signal track
mkdir /data/public_html_users/fankaili/IDEAS/DNAme_signal_hub/
cd /data/public_html_users/fankaili/IDEAS/DNAme_signal_hub/
#
cp /data/projects/encode/data/ENCSR286OOJ/ENCFF757KNC.bigWig ./stomach_0.bigWig
cp /data/projects/encode/data/ENCSR409HKJ/ENCFF046KKY.bigWig ./lung_0.bigWig
cp /data/projects/encode/data/ENCSR128HOP/ENCFF527LDJ.bigWig ./kidney_0.bigWig
cp /data/projects/encode/data/ENCSR353IFP/ENCFF017WVA.bigWig ./intestine_0.bigWig
cp /data/projects/encode/data/ENCSR874LLF/ENCFF648PDH.bigWig ./midbrain_0.bigWig
cp /data/projects/encode/data/ENCSR550CYA/ENCFF987BYS.bigWig ./liver_0.bigWig
cp /data/projects/encode/data/ENCSR397YEG/ENCFF667GAV.bigWig ./heart_0.bigWig
cp /data/projects/encode/data/ENCSR168RTO/ENCFF927KWX.bigWig ./hindbrain_0.bigWig
cp /data/projects/encode/data/ENCSR020GGM/ENCFF486IAI.bigWig ./forebrain_0.bigWig
cp /data/projects/encode/data/ENCSR191UKH/ENCFF087EHH.bigWig ./lung_14.5.bigWig
cp /data/projects/encode/data/ENCSR788XSZ/ENCFF269OVC.bigWig ./liver_14.5.bigWig
#
cp ../CTCF_signal_hub/genomes_e14.5p0_CTCF_signal.txt genomes_DNAme_signal.txt
vim genomes_DNAme_signal.txt
cp ../CTCF_signal_hub/hub_e14.5p0_CTCF_signal.txt hub_DNAme_signal.txt
vim hub_DNAme_signal.txt
cp ../CTCF_signal_hub/trackDb_e14.5p0_CTCF_signal.txt trackDb_DNAme_signal.txt
vim trackDb_DNAme_signal.txt
### https://users.wenglab.org/fankaili/IDEAS/DNAme_signal_hub/hub_DNAme_signal.txt
