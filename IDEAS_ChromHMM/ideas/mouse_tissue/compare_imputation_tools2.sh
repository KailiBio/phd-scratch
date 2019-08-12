#!bin/bash

# -- Kaili
# This script is for comparing imputation result with other tools. (given signal input in the same resolution)

# 0. make bigWig and signal file to match the resolution.
# 1.run ChromImpute
# 2. run Avocado
# 3. run IDEAS in 25bp bins
# 4. make track hub
# 5. calculate AUPR

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/imputation_comparison/"
oldSignalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_dhs_bins/"
signalDir="/data/zusers/fankaili/ideas/signal/rep1_signal_chr19_25bp_dhsResolution/"

cd ${workDir}

# 0. make bigWig and signal file to match the resolution.
for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC CTCF DNAme
do
    nohup bash ${scriptDir}get_signal_in_DHSbins_resolution.sh ${mark} > ./nohup/nohup.get_signal_in_DHSbins_resolution_${mark}.out 2>&1&
done
# z001: 14629-14639


# 1.run ChromImpute
## 1) make wig.gz file
for file in `ls ${signalDir}*.txt`
do
    name0=${file%_dhs.txt}
    name=${name0#/data/zusers/fankaili/ideas/signal/rep1_signal_chr19_25bp_dhsResolution/}
    echo $name
    #
    paste mm10_chr19_25bp.bed ${file} > tmp.txt
    #
    echo "track type=wiggle_0 name="${name} > /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/input_data0/${name}_dhs_25bp.wig
    echo "fixedStep  chrom=chr19 start=0 step=25 span=25" >> /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/input_data0/${name}_dhs_25bp.wig
    awk '{m=$2/25;a[m]=1;b[m]=$5}END{for(i=1;i<=2457263;i++){if(a[i]){print b[i]}else{print 0}}}' tmp.txt >> /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/input_data0/${name}_dhs_25bp.wig
    #
    gzip -c /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/input_data0/${name}_dhs_25bp.wig > /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/input_data/chr19_${name}_dhs_25bp.wig.gz
done
rm tmp.txt
## 2) get master table
if [ -f ./chromimpute_dhsResolution/master_table.txt ];then rm ./chromimpute_dhsResolution/master_table.txt; fi
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC CTCF DNAme
    do
        echo -e ${sample}"\t"${mark}"\t"${sample}"_"${mark}"_dhs_25bp" >> ./chromimpute_dhsResolution/master_table.txt
    done
done
vim ./chromimpute_dhsResolution/master_table.txt
## 3) run ChromImpute
cp ./chromimpute/mm10_chr19_chrom.size ./chromimpute_dhsResolution/
#
nohup bash ${scriptDir}run_chromimpute_dhsResolution.sh > ./nohup/nohup.run_chromimpute_dhsResolution.out 2>&1&
# z001: 26003


# 2. run Avocado
## convert signal file into compressed numpy
python ${scriptDir}convert_signal_compressed_numpy_dhsResolution.py
## run avocado using default
nohup python ${scriptDir}run_avocado_dhsResolution.sh > ./nohup/nohup.run_avocado_dhsResolution.out 2>&1&
# z001: 33878


# 3. run IDEAS in 25bp bins
# get bins
bedtools merge -i mm10_OCR_chr18-19_bins.bed > mm10_chr18-19.bed
bedtools makewindows -b mm10_chr18-19.bed -w 25 | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{print $0,"bin_"NR}' > mm10_chr18-19_25bp.bed
# get signal
cd /data/zusers/fankaili/ideas/
encode_data_path="/data/projects/encode/data/"
signal_path="/data/zusers/fankaili/ideas/signal/rep1_signal_chr18-19_25bp/"
mkdir ${signal_path}
#
if [ -f ENCODE_rep1_signal_chr18-19_25bp-bins_code.sh ]; then rm ENCODE_rep1_signal_chr18-19_25bp-bins_code.sh; fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "bigWigAverageOverBed "${encode_data_path}${expID}"/"${fileID}".bigWig /data/zusers/fankaili/ideas/imputation_comparison/IDEAS/mm10_chr18-19_25bp.bed "${signal_path}${sample}"_"${assay}"_chr19_25bp.tab" >> ENCODE_rep1_signal_chr18-19_25bp-bins_code.sh
    echo "awk '{print ""\$""5}' "${signal_path}${sample}"_"${assay}"_chr19_25bp.tab > "${signal_path}${sample}"_"${assay}"_chr19_25bp.txt" >> ENCODE_rep1_signal_chr18-19_25bp-bins_code.sh
done < ENCODE_mouse_rep1_CTCFsample_filelist1.txt
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    assay=`awk '{print $2}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    #
    echo "bigWigAverageOverBed /data/zusers/fankaili/ideas/signal/mouse_66samples_DNAme_bigWig/"${fileID}".bigWig /data/zusers/fankaili/ideas/imputation_comparison/IDEAS/mm10_chr18-19_25bp.bed "${signal_path}${sample}"_"${assay}"_chr19_25bp.tab" >> ENCODE_rep1_signal_chr18-19_25bp-bins_code.sh
    echo "awk '{print ""\$""6}' "${signal_path}${sample}"_"${assay}"_chr19_25bp.tab > "${signal_path}${sample}"_"${assay}"_chr19_25bp.txt" >> ENCODE_rep1_signal_chr18-19_25bp-bins_code.sh
done < ENCODE_mouse_rep1_CTCFsample_filelist2.txt
###
for i in {1..11}
do
    awk -v i="$i" '{if(NR>((i-1)*22) && NR<=(22*i)){print $0}}' ENCODE_rep1_signal_chr18-19_25bp-bins_code.sh > \
    /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_chr18-19_25bp-bins_code_${i}.sh
    nohup bash /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_chr18-19_25bp-bins_code_${i}.sh > \
    /data/zusers/fankaili/ideas/code/get_rep1_signal/nohup.ENCODE_rep1_signal_chr18-19_25bp-bins_code_${i}.out 2>&1&
done
# CTCFmotif
intersectBed -a mm10_chr18-19_25bp.bed -b /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/mm10_CTCF_motif_region.bed -wa -wb | cut -f 4,9 | sort -k1,1 -k2,2nr | awk 'BEGIN{FS=OFS="\t";bin=""}{if(NR==1){bin=$1;print $0}else{if(bin!=$1){bin=$1;print $0}}}' > mm10_chr18-19_25bp_bin_withCTCFmotifScore.txt
awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=1}else{if(b[$4]){print a[$4]}else{print 0}}}' mm10_chr18-19_25bp_bin_withCTCFmotifScore.txt mm10_chr18-19_25bp.bed > mm10_chr18-19_25bp_bin_CTCFmotifScore.txt
# CTCFaverage
awk '{FS=OFS="\t"}{print $4,0}' mm10_chr18-19_25bp.bed > chr18-19_25bp_9sample_sumCTCFsignal.txt
#
while read sample
do
    echo $sample
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{print $1,$2+a[$1]}}' /data/zusers/fankaili/ideas/signal/rep1_signal_chr18-19_25bp/${sample}_CTCF_chr19_25bp.tab chr18-19_25bp_9sample_sumCTCFsignal.txt > tmp.txt
    mv tmp.txt chr18-19_25bp_9sample_sumCTCFsignal.txt
done < /data/zusers/fankaili/ideas/dhs_ctcf/9_ctcf_sample.txt
#
awk '{FS=OFS="\t"}{print $2/9}' chr18-19_25bp_9sample_sumCTCFsignal.txt > chr18-19_25bp_9sample_averageCTCFsignal.txt
#####
cp ./IDEAS/test_IDEAS_chr18-19.sh ./IDEAS/test_IDEAS_chr18-19_25bp.sh
vim ./IDEAS/test_IDEAS_chr18-19_25bp.sh
cp ./IDEAS/test_IDEAS_chr18-19.parafile ./IDEAS/test_IDEAS_chr18-19_25bp.parafile
vim ./IDEAS/test_IDEAS_chr18-19_25bp.parafile
####
if [ -f ./IDEAS/test_IDEAS_chr18-19_25bp.input ];then rm ./IDEAS/test_IDEAS_chr18-19_25bp.input;fi
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 liver_14.5 lung_0 lung_14.5 stomach_0
do
    echo ${sample}
    for mark in H3K4me1 H3K4me2 H3K4me3 H3K9me3 H3K9ac H3K27me3 H3K27ac H3K36me3 ATAC DNAme
    do
        echo ${sample}" "${mark}" /data/zusers/fankaili/ideas/imputation_comparison/IDEAS/signal/"${sample}_${mark}_dhs_signal_chr18-19.txt >> ./IDEAS/test_IDEAS_chr18-19_25bp.input
    done
    #
    echo ${sample}" CTCFmotif /data/zusers/fankaili/ideas/imputation_comparison/IDEAS/mm10_chr18-19_25bp_bin_CTCFmotifScore.txt" >> ./IDEAS/test_IDEAS_chr18-19_25bp.input
    echo ${sample}" CTCFaverage /data/zusers/fankaili/ideas/imputation_comparison/IDEAS/chr18-19_25bp_9sample_averageCTCFsignal.txt" >> ./IDEAS/test_IDEAS_chr18-19_25bp.input
done
#
for sample in forebrain_0 midbrain_0 hindbrain_0 heart_0 intestine_0 kidney_0 liver_0 lung_0 stomach_0
do
    echo ${sample}" CTCF /data/zusers/fankaili/ideas/imputation_comparison/IDEAS/signal/"${sample}_CTCF_dhs_signal_chr18-19.txt >> ./IDEAS/test_IDEAS_chr18-19_25bp.input
done
########
awk '{OFS=" "}{print $1,$2,$3,$4}' mm10_chr18-19_25bp.bed > mm10_chr18-19_25bp_space.txt

nohup bash ./IDEAS/test_IDEAS_chr18-19_25bp.sh > ./nohup/nohup.test_IDEAS_chr18-19_25bp.out 2>&1&
# z008 5658



# 4. make track hub
## 1) get bigWig file
gzip -d /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/output/chr19_impute_liver_14.5_CTCF.wig.gz
wigToBigWig /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/output/chr19_impute_liver_14.5_CTCF.wig /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/mm10_chr19_chrom.size2 /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/output/chr19_impute_liver_14.5_CTCF.bw
cp /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/output/chr19_impute_liver_14.5_CTCF.bw /data/public_html_users/fankaili/IDEAS/other_imputation_DHSbins_resolution/ChromImpute_liver_14.5_sameResolution.bigWig
#
gzip -d /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/output/chr19_impute_lung_14.5_CTCF.wig.gz
wigToBigWig /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/output/chr19_impute_lung_14.5_CTCF.wig /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/mm10_chr19_chrom.size2 /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/output/chr19_impute_lung_14.5_CTCF.bw
cp /data/zusers/fankaili/ideas/imputation_comparison/chromimpute_dhsResolution/output/chr19_impute_lung_14.5_CTCF.bw /data/public_html_users/fankaili/IDEAS/other_imputation_DHSbins_resolution/ChromImpute_lung_14.5_sameResolution.bigWig
####
paste /data/zusers/fankaili/ideas/imputation_comparison/mm10_chr19_25bp.bed /data/zusers/fankaili/ideas/imputation_comparison/avocado_dhsResolution/avocado_predict_liver_14.5_CTCF.txt | cut -f 1-3,5 | sort -k1,1 -k2,2n > /data/zusers/fankaili/ideas/imputation_comparison/avocado_dhsResolution/avocado_predict_liver_14.5_CTCF.bedGraph
bedGraphToBigWig /data/zusers/fankaili/ideas/imputation_comparison/avocado_dhsResolution/avocado_predict_liver_14.5_CTCF.bedGraph /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/mm10_chr19_chrom.size2 /data/zusers/fankaili/ideas/imputation_comparison/avocado_dhsResolution/avocado_predict_liver_14.5_CTCF.bw
cp /data/zusers/fankaili/ideas/imputation_comparison/avocado_dhsResolution/avocado_predict_liver_14.5_CTCF.bw /data/public_html_users/fankaili/IDEAS/other_imputation_DHSbins_resolution/Avocado_liver_14.5_sameResolution.bigWig
#
paste /data/zusers/fankaili/ideas/imputation_comparison/mm10_chr19_25bp.bed /data/zusers/fankaili/ideas/imputation_comparison/avocado_dhsResolution/avocado_predict_lung_14.5_CTCF.txt | cut -f 1-3,5 | sort -k1,1 -k2,2n > /data/zusers/fankaili/ideas/imputation_comparison/avocado_dhsResolution/avocado_predict_lung_14.5_CTCF.bedGraph
bedGraphToBigWig /data/zusers/fankaili/ideas/imputation_comparison/avocado_dhsResolution/avocado_predict_lung_14.5_CTCF.bedGraph /data/zusers/fankaili/ideas/imputation_comparison/chromimpute/mm10_chr19_chrom.size2 /data/zusers/fankaili/ideas/imputation_comparison/avocado_dhsResolution/avocado_predict_lung_14.5_CTCF.bw
cp /data/zusers/fankaili/ideas/imputation_comparison/avocado_dhsResolution/avocado_predict_lung_14.5_CTCF.bw /data/public_html_users/fankaili/IDEAS/other_imputation_DHSbins_resolution/Avocado_lung_14.5_sameResolution.bigWig
## 2) get files
cp /data/public_html_users/fankaili/IDEAS/other_imputation/genomes.txt /data/public_html_users/fankaili/IDEAS/other_imputation_DHSbins_resolution/genomes.txt
vim /data/public_html_users/fankaili/IDEAS/other_imputation_DHSbins_resolution/genomes.txt
cp /data/public_html_users/fankaili/IDEAS/other_imputation/hub.txt /data/public_html_users/fankaili/IDEAS/other_imputation_DHSbins_resolution/hub.txt
vim /data/public_html_users/fankaili/IDEAS/other_imputation_DHSbins_resolution/hub.txt
vim /data/public_html_users/fankaili/IDEAS/other_imputation_DHSbins_resolution/trackDb.txt


# 5. calculate AUPR
bash ${scriptDir}calculate_AUPR_sameResolution.sh
