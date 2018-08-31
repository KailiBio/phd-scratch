#!/bin/bash

## Jul18.testIDEAS.sh
# 1. test IDEAS
# 2. compare IDEAS with ChromHMM

######################
## IDEAS

# 1. install IDEAS
cd /home/fankaili/git/
git clone https://github.com/guanjue/IDEAS_2018.git

# 2. run test sample
cd /data/zusers/fankaili/ideas/ss/

# 1)
cp /home/fankaili/git/IDEAS_2018/test_data/run_IDEAS.sh ./
cp /home/fankaili/git/IDEAS_2018/test_data/run_IDEAS.parafile ./
cp /home/fankaili/git/IDEAS_2018/test_data/run_IDEAS.input ./


# 2)
time bash run_IDEAS.sh

# 3. run our test data
cd /data/zusers/fankaili/ideas/test/

# 1) get .bed file
awk '{FS=OFS="\t"}{print $1}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed | sort -u
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr1"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr1_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr2"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr2_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr3"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr3_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr4"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr4_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr5"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr5_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr6"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr6_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr7"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr7_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr8"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr8_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr9"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr9_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr10"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr10_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr11"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr11_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr12"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr12_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr13"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr13_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr14"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr14_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr15"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr15_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr16"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr16_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr17"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr17_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr18"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr18_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chr19"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chr19_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chrX"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chrX_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chrY"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chrY_200bin.bed
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{if($1=="chrM"){i=i+1;print $1,$2,$3,R""i}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/mm10_chrM_200bin.bed

# 2) get signal file
cd /data/zusers/fankaili/ideas/test/
for i in {1..19} X Y M
do
    echo $i;
    awk '{FS=" ";OFS="\t"}{print $0,1,"."}' ./space_bed/mm10_chr${i}_200bin.bed > ./tab_bed6/mm10_chr${i}_200bin_tab.bed
done

cd /data/zusers/fankaili/ideas/test/
for file in `ls /data/public_html_users/vanderva/trackhub/chromhmmpaper/binned/dnase/`
do
    filename=${file%.bigWig}
    for i in {1..19} X Y M
    do
        bigWigAverageOverBed /data/public_html_users/vanderva/trackhub/chromhmmpaper/binned/dnase/${file} ./tab_bed6/mm10_chr${i}_200bin_tab.bed ./signal/${filename}_chr${i}.tab ;
        awk '{print $5}' ./signal/${filename}_chr${i}.tab > ./signal/${filename}_chr${i}.txt ;
    done
done

# 3) run on Slurm

sed 's/chr1/chr2/g' run_IDEAS_chr1.sh > run_IDEAS_chr2.sh
sed 's/chr1/chr3/g' run_IDEAS_chr1.sh > run_IDEAS_chr3.sh
sed 's/chr1/chr4/g' run_IDEAS_chr1.sh > run_IDEAS_chr4.sh
sed 's/chr1/chr5/g' run_IDEAS_chr1.sh > run_IDEAS_chr5.sh
sed 's/chr1/chr6/g' run_IDEAS_chr1.sh > run_IDEAS_chr6.sh
sed 's/chr1/chr7/g' run_IDEAS_chr1.sh > run_IDEAS_chr7.sh
sed 's/chr1/chr8/g' run_IDEAS_chr1.sh > run_IDEAS_chr8.sh
sed 's/chr1/chr9/g' run_IDEAS_chr1.sh > run_IDEAS_chr9.sh
sed 's/chr1/chr10/g' run_IDEAS_chr1.sh > run_IDEAS_chr10.sh
sed 's/chr1/chr11/g' run_IDEAS_chr1.sh > run_IDEAS_chr11.sh
sed 's/chr1/chr12/g' run_IDEAS_chr1.sh > run_IDEAS_chr12.sh
sed 's/chr1/chr13/g' run_IDEAS_chr1.sh > run_IDEAS_chr13.sh
sed 's/chr1/chr14/g' run_IDEAS_chr1.sh > run_IDEAS_chr14.sh
sed 's/chr1/chr15/g' run_IDEAS_chr1.sh > run_IDEAS_chr15.sh
sed 's/chr1/chr16/g' run_IDEAS_chr1.sh > run_IDEAS_chr16.sh
sed 's/chr1/chr17/g' run_IDEAS_chr1.sh > run_IDEAS_chr17.sh
sed 's/chr1/chr18/g' run_IDEAS_chr1.sh > run_IDEAS_chr18.sh
sed 's/chr1/chr19/g' run_IDEAS_chr1.sh > run_IDEAS_chr19.sh
sed 's/chr1/chrX/g' run_IDEAS_chr1.sh > run_IDEAS_chr20.sh
sed 's/chr1/chrY/g' run_IDEAS_chr1.sh > run_IDEAS_chr21.sh
sed 's/chr1/chrM/g' run_IDEAS_chr1.sh > run_IDEAS_chr22.sh

sed 's/chr1/chr2/g' run_IDEAS_chr1.input > run_IDEAS_chr2.input
sed 's/chr1/chr3/g' run_IDEAS_chr1.input > run_IDEAS_chr3.input
sed 's/chr1/chr4/g' run_IDEAS_chr1.input > run_IDEAS_chr4.input
sed 's/chr1/chr5/g' run_IDEAS_chr1.input > run_IDEAS_chr5.input
sed 's/chr1/chr6/g' run_IDEAS_chr1.input > run_IDEAS_chr6.input
sed 's/chr1/chr7/g' run_IDEAS_chr1.input > run_IDEAS_chr7.input
sed 's/chr1/chr8/g' run_IDEAS_chr1.input > run_IDEAS_chr8.input
sed 's/chr1/chr9/g' run_IDEAS_chr1.input > run_IDEAS_chr9.input
sed 's/chr1/chr10/g' run_IDEAS_chr1.input > run_IDEAS_chr10.input
sed 's/chr1/chr11/g' run_IDEAS_chr1.input > run_IDEAS_chr11.input
sed 's/chr1/chr12/g' run_IDEAS_chr1.input > run_IDEAS_chr12.input
sed 's/chr1/chr13/g' run_IDEAS_chr1.input > run_IDEAS_chr13.input
sed 's/chr1/chr14/g' run_IDEAS_chr1.input > run_IDEAS_chr14.input
sed 's/chr1/chr15/g' run_IDEAS_chr1.input > run_IDEAS_chr15.input
sed 's/chr1/chr16/g' run_IDEAS_chr1.input > run_IDEAS_chr16.input
sed 's/chr1/chr17/g' run_IDEAS_chr1.input > run_IDEAS_chr17.input
sed 's/chr1/chr18/g' run_IDEAS_chr1.input > run_IDEAS_chr18.input
sed 's/chr1/chr19/g' run_IDEAS_chr1.input > run_IDEAS_chr19.input
sed 's/chr1/chrX/g' run_IDEAS_chr1.input > run_IDEAS_chrX.input
sed 's/chr1/chrY/g' run_IDEAS_chr1.input > run_IDEAS_chrY.input
sed 's/chr1/chrM/g' run_IDEAS_chr1.input > run_IDEAS_chrM.input

sed 's/chr1/chr2/g' run_IDEAS_chr1.parafile > run_IDEAS_chr2.parafile
sed 's/chr1/chr3/g' run_IDEAS_chr1.parafile > run_IDEAS_chr3.parafile
sed 's/chr1/chr4/g' run_IDEAS_chr1.parafile > run_IDEAS_chr4.parafile
sed 's/chr1/chr5/g' run_IDEAS_chr1.parafile > run_IDEAS_chr5.parafile
sed 's/chr1/chr6/g' run_IDEAS_chr1.parafile > run_IDEAS_chr6.parafile
sed 's/chr1/chr7/g' run_IDEAS_chr1.parafile > run_IDEAS_chr7.parafile
sed 's/chr1/chr8/g' run_IDEAS_chr1.parafile > run_IDEAS_chr8.parafile
sed 's/chr1/chr9/g' run_IDEAS_chr1.parafile > run_IDEAS_chr9.parafile
sed 's/chr1/chr10/g' run_IDEAS_chr1.parafile > run_IDEAS_chr10.parafile
sed 's/chr1/chr11/g' run_IDEAS_chr1.parafile > run_IDEAS_chr11.parafile
sed 's/chr1/chr12/g' run_IDEAS_chr1.parafile > run_IDEAS_chr12.parafile
sed 's/chr1/chr13/g' run_IDEAS_chr1.parafile > run_IDEAS_chr13.parafile
sed 's/chr1/chr14/g' run_IDEAS_chr1.parafile > run_IDEAS_chr14.parafile
sed 's/chr1/chr15/g' run_IDEAS_chr1.parafile > run_IDEAS_chr15.parafile
sed 's/chr1/chr16/g' run_IDEAS_chr1.parafile > run_IDEAS_chr16.parafile
sed 's/chr1/chr17/g' run_IDEAS_chr1.parafile > run_IDEAS_chr17.parafile
sed 's/chr1/chr18/g' run_IDEAS_chr1.parafile > run_IDEAS_chr18.parafile
sed 's/chr1/chr19/g' run_IDEAS_chr1.parafile > run_IDEAS_chr19.parafile
sed 's/chr1/chrX/g' run_IDEAS_chr1.parafile > run_IDEAS_chrX.parafile
sed 's/chr1/chrY/g' run_IDEAS_chr1.parafile > run_IDEAS_chrY.parafile
sed 's/chr1/chrM/g' run_IDEAS_chr1.parafile > run_IDEAS_chrM.parafile

rm test_IDEAS.err test_IDEAS.out
sbatch test_IDEAS.sh

# 4. run on every two chr
####### Jul20
cd /data/zusers/fankaili/ideas/test/

# 1) get .bed file
awk '{FS=OFS="\t"}{print $1}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed | sort -u
#
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr1"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr1_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr2"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr2_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr3"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr3_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr4"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr4_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr5"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr5_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr6"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr6_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr7"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr7_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr8"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr8_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr9"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr9_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr10"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr10_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr11"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr11_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr12"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr12_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr13"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr13_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr14"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr14_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr15"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr15_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr16"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr16_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr17"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr17_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr18"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr18_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr19"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr19_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chrX"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chrX_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chrY"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chrY_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chrM"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chrM_200bin.bed
#
cat ./space_bed/tmp_mm10_chr1_200bin.bed ./space_bed/tmp_mm10_chr2_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test1.bed
cat ./space_bed/tmp_mm10_chr3_200bin.bed ./space_bed/tmp_mm10_chr4_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test2.bed
cat ./space_bed/tmp_mm10_chr5_200bin.bed ./space_bed/tmp_mm10_chr6_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test3.bed
cat ./space_bed/tmp_mm10_chr7_200bin.bed ./space_bed/tmp_mm10_chr8_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test4.bed
cat ./space_bed/tmp_mm10_chr9_200bin.bed ./space_bed/tmp_mm10_chr10_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test5.bed
cat ./space_bed/tmp_mm10_chr11_200bin.bed ./space_bed/tmp_mm10_chr12_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test6.bed
cat ./space_bed/tmp_mm10_chr13_200bin.bed ./space_bed/tmp_mm10_chr14_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test7.bed
cat ./space_bed/tmp_mm10_chr15_200bin.bed ./space_bed/tmp_mm10_chr16_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test8.bed
cat ./space_bed/tmp_mm10_chr17_200bin.bed ./space_bed/tmp_mm10_chr18_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test9.bed
cat ./space_bed/tmp_mm10_chr19_200bin.bed ./space_bed/tmp_mm10_chrM_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test10.bed
cat ./space_bed/tmp_mm10_chrX_200bin.bed ./space_bed/tmp_mm10_chrY_200bin.bed | awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' > ./space_bed/mm10_200bin_test11.bed
#
rm ./space_bed/tmp_*.bed

# 2) get signal file
cd /data/zusers/fankaili/ideas/test/
for i in {1..11}
do
    echo $i;
    awk '{FS=" ";OFS="\t"}{print $0,1,"."}' ./space_bed/mm10_200bin_test${i}.bed > ./tab_bed6/mm10_200bin_test${i}_tab.bed
done

cd /data/zusers/fankaili/ideas/test/
for file in `ls /data/public_html_users/vanderva/trackhub/chromhmmpaper/binned/dnase/`
do
    filename=${file%.bigWig}
    for i in {1..11}
    do
        bigWigAverageOverBed /data/public_html_users/vanderva/trackhub/chromhmmpaper/binned/dnase/${file} ./tab_bed6/mm10_200bin_test${i}_tab.bed ./signal/${filename}_test${i}.tab ;
        awk '{print $5}' ./signal/${filename}_test${i}.tab > ./signal/${filename}_test${i}.txt ;
    done
done

# 3) file for IDEAS
cd /data/zusers/fankaili/ideas/test/

sed 's/test1/test2/g' run_IDEAS_test1.sh > run_IDEAS_test2.sh
sed 's/test1/test3/g' run_IDEAS_test1.sh > run_IDEAS_test3.sh
sed 's/test1/test4/g' run_IDEAS_test1.sh > run_IDEAS_test4.sh
sed 's/test1/test5/g' run_IDEAS_test1.sh > run_IDEAS_test5.sh
sed 's/test1/test6/g' run_IDEAS_test1.sh > run_IDEAS_test6.sh
sed 's/test1/test7/g' run_IDEAS_test1.sh > run_IDEAS_test7.sh
sed 's/test1/test8/g' run_IDEAS_test1.sh > run_IDEAS_test8.sh
sed 's/test1/test9/g' run_IDEAS_test1.sh > run_IDEAS_test9.sh
sed 's/test1/test10/g' run_IDEAS_test1.sh > run_IDEAS_test10.sh
sed 's/test1/test11/g' run_IDEAS_test1.sh > run_IDEAS_test11.sh

sed 's/test1/test2/g' run_IDEAS_test1.input > run_IDEAS_test2.input
sed 's/test1/test3/g' run_IDEAS_test1.input > run_IDEAS_test3.input
sed 's/test1/test4/g' run_IDEAS_test1.input > run_IDEAS_test4.input
sed 's/test1/test5/g' run_IDEAS_test1.input > run_IDEAS_test5.input
sed 's/test1/test6/g' run_IDEAS_test1.input > run_IDEAS_test6.input
sed 's/test1/test7/g' run_IDEAS_test1.input > run_IDEAS_test7.input
sed 's/test1/test8/g' run_IDEAS_test1.input > run_IDEAS_test8.input
sed 's/test1/test9/g' run_IDEAS_test1.input > run_IDEAS_test9.input
sed 's/test1/test10/g' run_IDEAS_test1.input > run_IDEAS_test10.input
sed 's/test1/test11/g' run_IDEAS_test1.input > run_IDEAS_test11.input

sed 's/test1/test2/g' run_IDEAS_test1.parafile > run_IDEAS_test2.parafile
sed 's/test1/test3/g' run_IDEAS_test1.parafile > run_IDEAS_test3.parafile
sed 's/test1/test4/g' run_IDEAS_test1.parafile > run_IDEAS_test4.parafile
sed 's/test1/test5/g' run_IDEAS_test1.parafile > run_IDEAS_test5.parafile
sed 's/test1/test6/g' run_IDEAS_test1.parafile > run_IDEAS_test6.parafile
sed 's/test1/test7/g' run_IDEAS_test1.parafile > run_IDEAS_test7.parafile
sed 's/test1/test8/g' run_IDEAS_test1.parafile > run_IDEAS_test8.parafile
sed 's/test1/test9/g' run_IDEAS_test1.parafile > run_IDEAS_test9.parafile
sed 's/test1/test10/g' run_IDEAS_test1.parafile > run_IDEAS_test10.parafile
sed 's/test1/test11/g' run_IDEAS_test1.parafile > run_IDEAS_test11.parafile

# 4) run on Slurm
cd /home/fankaili/test_IDEAS/
rm test_IDEAS.err test_IDEAS.out
sbatch test_IDEAS.sh


# 5.
cd /data/zusers/fankaili/ideas/test/
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr1"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr1_200bin.bed
awk 'BEGIN{FS==OFS="\t"}{if($1=="chr2"){print $1,$2,$3}}' /data/zusers/vanderva/vanderva.zlab3/code/git/segmentation/output/chromhmm_18_state_8marks_dname_atac_dnase_0.3/processed/embryonic_facial_prominence_11.5_mm10_18_posterior.bed \
> ./space_bed/tmp_mm10_chr2_200bin.bed

head -972496 ./space_bed/tmp_mm10_chr1_200bin.bed > mm10_200bin_ss.txt
head -27504 ./space_bed/tmp_mm10_chr2_200bin.bed >> mm10_200bin_ss.txt
awk -v R="R" 'BEGIN{FS="\t";OFS=" ";i=0}{i=i+1;print $1,$2,$3,R""i}' mm10_200bin_ss.txt > mm10_200bin_ss.bed
awk '{FS=" ";OFS="\t"}{print $0,1,"."}' mm10_200bin_ss.bed > ./tab_bed6/mm10_200bin_ss_tab.bed

cp run_IDEAS_test1.input run_IDEAS_ss.input
cp run_IDEAS_test1.parafile run_IDEAS_ss.parafile
cp run_IDEAS_test1.sh run_IDEAS_ss.sh
# modify by vim

cd /data/zusers/fankaili/ideas/test/
for file in `ls /data/public_html_users/vanderva/trackhub/chromhmmpaper/binned/dnase/`
do
    filename=${file%.bigWig}
    bigWigAverageOverBed /data/public_html_users/vanderva/trackhub/chromhmmpaper/binned/dnase/${file} ./tab_bed6/mm10_200bin_ss_tab.bed ./signal/${filename}_ss.tab ;
    awk '{print $5}' ./signal/${filename}_ss.tab > ./signal/${filename}_ss.txt ;
done


rm test_IDEAS_ss.err test_IDEAS_ss.out
sbatch test_IDEAS_ss.sh
