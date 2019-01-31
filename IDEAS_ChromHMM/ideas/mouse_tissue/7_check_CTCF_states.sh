#!bin/bash

# -- Kaili
# This script is for checking and analyzing real CTCF states.
# 1. decide which state is real CTCF states.
# 2. get 11 candidate states in liver14.5
# 3. 11 candidate states with CTCF peak (liver14.5)
# 4. 11 candidate states with CTCF motif (liver14.5)
# 5. remake contigency table and Venn
# 6. composition of the 11states_noPeak/motif
# 7. for all states
# 8. state conservation
# 9. for 6 states
# 10. percentage of states with CTCF peaks, running cut-off

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
dailyCodeDir="/data/zusers/fankaili/github/weng-lab/Kaili/DailyCode/"
workDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/"

cd ${workDir}

# 1. decide which state is real CTCF states.
## get state bed for each sample
mkdir state_bed
#
stateDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/"
prefix="ctcf_samples."
stateBedDir="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_bed/"
nohup bash ${dailyCodeDir}make_each_sample_state_bed_IDEAS.sh ${stateDir} ${prefix} ${stateBedDir} 11 \
> /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/nohup.make_each_sample_state_bed.out 2>&1&
# get CTCF signal for each state (43 states in this run)
mkdir state_CTCF_signal

for ((state=0; state<=42; state++))
do
    echo ${state}
    bash ${scriptDir}get_CTCF_signal_for_given_state.sh /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ ${state}
done

# Rscript make_boxplot_CTCF_signal_for_all_states.R

# average signal of each states
if [ -f ./state_CTCF_signal/state_CTCF_signal_ave.txt ]; then rm ./state_CTCF_signal/state_CTCF_signal_ave.txt; fi
for ((state=0; state<=42; state++))
do
    echo ${state}
    awk -v state="$state" 'BEGIN{FS=OFS="\t";sum=0}{sum+=$4}END{print state, sum/NR}' \
    ./state_CTCF_signal/state_${state}_CTCF_signal.txt >> ./state_CTCF_signal/state_CTCF_signal_ave.txt
done
sort -k2,2nr ./state_CTCF_signal/state_CTCF_signal_ave.txt > ./state_CTCF_signal/state_CTCF_signal_ave_sorted.txt

# 2. get 11 candidate states in liver14.5
mkdir ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region
for i in 21 29 30 40 38 42 27 11 28 26 31
do
    echo ${i}
    awk '{FS=OFS="\t"}{if(NR==FNR){if($2=="liver_14.5"){a[$3]=1}}else{if(a[$4]){print $0}}}' \
    ./state_CTCF_signal/state_${i}_CTCF_signal.txt mm10_tab.bed | sort -k1,1 -k2,2n > \
    ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_${i}_region.bed
done

## average CTCF signal of 11 states in liver14.5
if [ -f ./state_CTCF_signal/liver14.5_state_CTCF_signal_ave.txt ];
then
    rm ./state_CTCF_signal/liver14.5_state_CTCF_signal_ave.txt
fi
for ((state=0; state<=42; state++))
do
    echo ${state}
    awk -v state="$state" 'BEGIN{FS=OFS="\t";sum=0;n=0}{if($2=="liver_14.5"){sum+=$4;n+=1}}END{print state, sum/n}' \
    ./state_CTCF_signal/state_${state}_CTCF_signal.txt >> ./state_CTCF_signal/liver14.5_state_CTCF_signal_ave.txt
done
sort -k2,2nr ./state_CTCF_signal/liver14.5_state_CTCF_signal_ave.txt> ./state_CTCF_signal/liver14.5_state_CTCF_signal_ave_sorted.txt


# 3. 11 candidate states with CTCF peak (liver14.5)
if [ -f ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11state_peak.txt ]
then
    rm ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11state_peak.txt
fi
#
for i in 21 29 30 40 38 42 27 11 28 26 31
do
    echo ${i}
    total=`wc -l ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_${i}_region.bed | awk '{print $1}'`
    num=`intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_${i}_region.bed \
    -b liver_14.5_day_CTCF_peak_sorted.bed -f 1 -F 0.5 -e -wa | sort -u | wc -l | awk '{print $1}'`
    echo -e ${i}"\t"${total}"\t"${num} >> ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11state_peak.txt
done
##### only with strong peak
if [ -f ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11state_peak32.txt ]
then
    rm ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11state_peak32.txt
fi
#
for i in 21 29 30 40 38 42 27 11 28 26 31
do
    echo ${i}
    total=`wc -l ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_${i}_region.bed | awk '{print $1}'`
    num=`intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_${i}_region.bed \
    -b liver_14.5_day_CTCF_peak_over32.bed -f 1 -F 0.5 -e -wa | sort -u | wc -l | awk '{print $1}'`
    echo -e ${i}"\t"${total}"\t"${num} >> ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11state_peak32.txt
done




# 4. 11 candidate states with CTCF motif (liver14.5)
if [ -f ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11state_motif.txt ]
then
    rm ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11state_motif.txt
fi
#
for i in 21 29 30 40 38 42 27 11 28 26 31
do
    echo ${i}
    total=`wc -l ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_${i}_region.bed | awk '{print $1}'`
    num=`intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_${i}_region.bed \
    -b mm10_CTCF_motif_region.bed -wa | sort -u | wc -l | awk '{print $1}'`
    echo -e ${i}"\t"${total}"\t"${num} >> ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11state_motif.txt
done

# Rscript overlapped_11states_CTCF_peak_motif.R

# 5. remake contigency table and Venn
cat ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_*_region.bed > \
./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_region.bed
#
count_contigency(){
    ctcf_bin_file=$1
    all_bin_file=$2
    interested_file=$3
    outpath=$4
    prefix=$5
    signal_file=$6
    #
    a=`wc -l ${ctcf_bin_file} | awk '{print $1}'`
    b=`wc -l ${all_bin_file} | awk '{print $1}'`
    #  get bed
    intersectBed -a ${ctcf_bin_file} -b ${interested_file} -f 1 -F 0.5 -e -wa -u | sort -k1,1 -k2,2n \
    > ${outpath}${prefix}with.bed
    intersectBed -a ${ctcf_bin_file} -b ${interested_file} -f 1 -F 0.5 -e -wa -v | sort -u | sort -k1,1 -k2,2n \
    > ${outpath}${prefix}without.bed
    intersectBed -a ${all_bin_file} -b ${ctcf_bin_file} -wa -v | sort -u > tmp_count_contigency.bed
    intersectBed -a tmp_count_contigency.bed -b ${interested_file} -f 1 -F 0.5 -e -wa -u | sort -k1,1 -k2,2n \
    > ${outpath}${prefix}non_with.bed
    # echo number
    n_11_withpeak=`wc -l ${outpath}${prefix}with.bed | awk '{print $1}'`
    n_11_withoutpeak=`wc -l ${outpath}${prefix}without.bed | awk '{print $1}'`
    n_no11_withpeak=`wc -l ${outpath}${prefix}non_with.bed | awk '{print $1}'`
    echo -e ${n_11_withpeak}"\t"${n_11_withoutpeak}"\t"${n_no11_withpeak}"\t"${a}"\t"${b}
    # get signal
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' \
    ${outpath}${prefix}with.bed ${signal_file} > ${outpath}${prefix}with_signal.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' \
    ${outpath}${prefix}without.bed ${signal_file} > ${outpath}${prefix}without_signal.txt
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$1]){print $1,$5}}}' \
    ${outpath}${prefix}non_with.bed ${signal_file} > ${outpath}${prefix}non_with_signal.txt
}

ctcf_bin_file="./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_region.bed"
all_bin_file="mm10_tab.bed"
peak_file="liver_14.5_day_CTCF_peak_sorted.bed"
motif_file="mm10_CTCF_motif_region.bed"
outpath="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_CTCF_signal/liver14.5_CTCF_state_candidate_region/"
peak_prefix="liver14.5_11states_peak_"
motif_prefix="liver14.5_11states_motif_"
signal_file="/data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/liver_14.5_CTCF_normal.tab"

count_contigency ${ctcf_bin_file} ${all_bin_file} ${peak_file} ${outpath} ${peak_prefix} ${signal_file}
count_contigency ${ctcf_bin_file} ${all_bin_file} ${motif_file} ${outpath} ${motif_prefix} ${signal_file}

# Rscript make_CTCF_signal_3group_boxplot.R

#### Venn
intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_peak_with.bed \
-b mm10_CTCF_motif_region.bed -f 1 -F 0.5 -e -wa -u | wc -l
intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_peak_with.bed \
-b mm10_CTCF_motif_region.bed -f 1 -F 0.5 -e -wa -v | sort -u | wc -l
intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_motif_with.bed \
-b liver_14.5_day_CTCF_peak_sorted.bed -f 1 -F 0.5 -e -wa -v | sort -u | wc -l
#
intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_peak_non_with.bed \
-b mm10_CTCF_motif_region.bed -f 1 -F 0.5 -e -wa -u | wc -l
intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_peak_non_with.bed \
-b mm10_CTCF_motif_region.bed -f 1 -F 0.5 -e -wa -v | sort -u | wc -l

# 6. composition of the 11states_noPeak/motif
if [ -f ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_noPeak_composition.txt ]
then
    rm ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_noPeak_composition.txt
fi
#
if [ -f ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_noMotif_composition.txt ]
then
    rm ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_noMotif_composition.txt
fi
for i in 21 29 30 40 38 42 27 11 28 26 31
do
    echo ${i}
    # peak
    num_peak=`cat ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_peak_without.bed \
    ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_${i}_region.bed | cut -f 4 | sort | uniq -d \
    | wc -l | awk '{print $1}'`
    echo -e ${i}"\t"${num_peak} >> ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_noPeak_composition.txt
    # motif
    num_motif=`cat ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_motif_without.bed \
    ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_${i}_region.bed | cut -f 4 | sort | uniq -d \
    | wc -l | awk '{print $1}'`
    echo -e ${i}"\t"${num_motif} >> ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_noMotif_composition.txt
done

# Rscript make_composition_11states_noPeakMotif_pie.R

# 7. for all states
## 1) get region
mkdir ./state_CTCF_signal/liver14.5_all_states_region
for i in {0..42}
do
    echo ${i}
    awk '{FS=OFS="\t"}{if(NR==FNR){if($2=="liver_14.5"){a[$3]=1}}else{if(a[$4]){print $0}}}' \
    ./state_CTCF_signal/state_${i}_CTCF_signal.txt mm10_tab.bed | sort -k1,1 -k2,2n > \
    ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_state_${i}_region.bed
done
## 2)  peak
if [ -f ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_all_states_peak.txt ]
then
    rm ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_all_states_peak.txt
fi
#
for i in {0..42}
do
    echo ${i}
    total=`wc -l ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_state_${i}_region.bed | awk '{print $1}'`
    num=`intersectBed -a ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_state_${i}_region.bed \
    -b liver_14.5_day_CTCF_peak_sorted.bed -f 1 -F 0.5 -e -wa | sort -u | wc -l | awk '{print $1}'`
    echo -e ${i}"\t"${total}"\t"${num} >> ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_all_states_peak.txt
done
## 3) motif
if [ -f ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_all_states_motif.txt ]
then
    rm ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_all_states_motif.txt
fi
#
for i in {0..42}
do
    echo ${i}
    total=`wc -l ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_state_${i}_region.bed | awk '{print $1}'`
    num=`intersectBed -a ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_state_${i}_region.bed \
    -b mm10_CTCF_motif_region.bed -wa | sort -u | wc -l | awk '{print $1}'`
    echo -e ${i}"\t"${total}"\t"${num} >> ./state_CTCF_signal/liver14.5_all_states_region/liver14.5_all_states_motif.txt
done

# Rscript make_percentage_all_states_peak_motif_barplot.R


# 8. state conservation
## 1) get all states file for all celltypes
for j in {5..15}
do
    echo ${j}
    celltype=`head -1 /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.chr1.state | awk -v i="${j}" '{FS=" "}{print $i}'`
    #
    for s in {0..42}
    do
        if [ -f ./ctcf_states_rep1/mm10_state_${s}_${celltype}.bed ];then
            rm ./ctcf_states_rep1/mm10_state_${s}_${celltype}.bed
        fi
        #
        for i in {1..19} X Y
        do
            awk -v s="$s" -v j="${j}" '{FS=" ";OFS="\t"}{if($j==s){print $2,$3,$4,$1}}' \
            /data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ctcf_samples_result/ctcf_samples.chr${i}.state \
            >> ./ctcf_states_rep1/mm10_state_${s}_${celltype}.bed
        done
        cut -f 1-4 ./ctcf_states_rep1/mm10_state_${s}_${celltype}.bed | sort -u | sort -k1,1 -k2,2n > ./ctcf_states_rep1/mm10_state_${s}_${celltype}_sorted.bed
    done
done

## 2) get CTCF states overlapped matrix.
### need a long time for this step
for s in {0..42}
do
    echo ${s}
    echo -e "lung_14.5\tliver_14.5\tstomach_0\tmidbrain_0\tkidney_0\tliver_0\tintestine_0\tlung_0\theart_0\thindbrain_0\tforebrain_0" \
    > ./ctcf_states_rep1/state_${s}_overlapped.txt
    echo -e "lung_14.5\tliver_14.5\tstomach_0\tmidbrain_0\tkidney_0\tliver_0\tintestine_0\tlung_0\theart_0\thindbrain_0\tforebrain_0" \
    > ./ctcf_states_rep1/state_${s}_union.txt
    #
    for type1 in lung_14.5 liver_14.5 stomach_0 midbrain_0 kidney_0 liver_0 intestine_0 lung_0 heart_0 hindbrain_0 forebrain_0
    do
        echo ${type1}
        for type2 in lung_14.5 liver_14.5 stomach_0 midbrain_0 kidney_0 liver_0 intestine_0 lung_0 heart_0 hindbrain_0 forebrain_0
        do
            num=`cat ./ctcf_states_rep1/mm10_state_${s}_${type1}_sorted.bed ./ctcf_states_rep1/mm10_state_${s}_${type2}_sorted.bed | cut -f 4 | sort | uniq -d | wc -l`
            echo -n -e $num"\t" >> ./ctcf_states_rep1/state_${s}_overlapped.txt
            union=`cat ./ctcf_states_rep1/mm10_state_${s}_${type1}_sorted.bed ./ctcf_states_rep1/mm10_state_${s}_${type2}_sorted.bed | cut -f 4 | sort -u | wc -l`
            echo -n -e $union"\t" >> ./ctcf_states_rep1/state_${s}_union.txt
        done
        echo -n -e "\n" >> ./ctcf_states_rep1/state_${s}_overlapped.txt
        echo -n -e "\n" >> ./ctcf_states_rep1/state_${s}_union.txt
    done
done

# Rscript make_jaccard_heatmap_each_state.R


# 9. for 6 states
cat ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_21_region.bed \
./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_29_region.bed \
./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_30_region.bed \
./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_40_region.bed \
./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_42_region.bed \
./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_state_11_region.bed> \
./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_6states_region.bed
# contigency
ctcf_bin_file="./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_6states_region.bed"
all_bin_file="mm10_tab.bed"
peak_file="liver_14.5_day_CTCF_peak_sorted.bed"
motif_file="mm10_CTCF_motif_region.bed"
outpath="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/state_CTCF_signal/liver14.5_CTCF_state_candidate_region/"
peak_prefix="liver14.5_6states_peak_"
motif_prefix="liver14.5_6states_motif_"
signal_file="/data/zusers/fankaili/ideas/signal/rep1_signal_normal_bins/liver_14.5_CTCF_normal.tab"

count_contigency ${ctcf_bin_file} ${all_bin_file} ${peak_file} ${outpath} ${peak_prefix} ${signal_file}
count_contigency ${ctcf_bin_file} ${all_bin_file} ${motif_file} ${outpath} ${motif_prefix} ${signal_file}
#### Venn
intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_peak_with.bed \
-b mm10_CTCF_motif_region.bed -f 0.5 -F 0.5 -e -wa -u | wc -l
intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_peak_with.bed \
-b mm10_CTCF_motif_region.bed -f 0.5 -F 0.5 -e -wa -v | sort -u | wc -l
intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_motif_with.bed \
-b liver_14.5_day_CTCF_peak_sorted.bed -f 0.5 -F 0.5 -e -wa -v | sort -u | wc -l
#
intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_peak_non_with.bed \
-b mm10_CTCF_motif_region.bed -f 0.5 -F 0.5 -e -wa -u | wc -l
intersectBed -a ./state_CTCF_signal/liver14.5_CTCF_state_candidate_region/liver14.5_11states_peak_non_with.bed \
-b mm10_CTCF_motif_region.bed -f 0.5 -F 0.5 -e -wa -v | sort -u | wc -l

# Rscript make_CTCF_signal_3group_boxplot.R


# 10. percentage of states with CTCF peaks, running cut-off
bash ${scriptDir}get_state_percentage_of_CTCF_peaks_with_running_cutoff.sh
