#!/bin/bash

# -- Kaili
# This script is for running Juicer Loop on GHPCC.
# from install HICCUPS


# 1. install HICCUPS
# https://docs.google.com/document/d/10599KoPp_C8lEcAdNoGbTpXm24WmalxXdE9Vpn-htIA/edit
cd ~
mkdir github

## 1) get Juicerbox
git clone https://github.com/theaidenlab/juicebox

## 2) get apache-ant
# http://dita-ot.sourceforge.net/doc/ot-userguide13/xhtml/installing/linux_installingant.html
wget http://mirror.metrocast.net/apache//ant/binaries/apache-ant-1.9.13-bin.tar.gz
tar -xzf apache-ant-1.9.13-bin.tar.gz

## 3) env variable setting for build juicebox jars

cd /home/kf27w/github/Juicebox/native_launcher/natives/
cp /home/kf27w/github/Juicebox/lib/jcuda/Archive.JCuda.0.7.0.zip ./
unzip Archive.JCuda.0.7.0.zip

##### modify Juicebox/build.xml
vim /home/kf27w/github/Juicebox/build.xml
# <zipfileset src="${basedir}/lib/jcuda/jcuda-0.7.0.jar"/>
# <zipfileset src="${basedir}/lib/jcuda/jcudaUtils-0.0.4.jar"/>

##### modify Juicebox/juicebox.properties
vim /home/kf27w/github/Juicebox/juicebox.properties
# jdk.home.1.8=/share/pkg/jdk/1.8.0_31




# 2. compile jar from source
cd /home/kf27w/github/Juicebox/
bsub -q interactive -W 6:00 -R rusage[mem=50000] -Is bash
module load jdk/1.8.0_31
export ANT_OPTS="-Xmx256M"
export ANT_HOME="/home/kf27w/github/apache-ant-1.9.13"
export PATH=${ANT_HOME}/bin:${JAVA_HOME}/bin:${PATH}
ant
cp out/artifacts/Juicebox_clt_jar/Juicebox.jar native_launcher/Juicebox_CLT.jar



# 3. test
cd /home/kf27w/github/Juicebox/native_launcher/
bsub -q interactive -q gpu_geforceGTX  -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]"
module load jdk/1.8.0_31
module load gcc/4.7.4
module load cuda/7.0.28

mkdir test
bash juicebox.sh hiccups --ignore_sparsity -r 10000 -f 0.1 -p 2 -i 5 -d 20000 -c 8 \
https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic test/


# 4. run HICCUPS
cd /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/
bsub -W 12:00 -q gpu_geforceGTX -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]" -n 1 /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/run_hiccups_cn_1.sh
bsub -W 12:00 -q gpu_geforceGTX -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]" -n 1 /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/run_hiccups_cn_2.sh
bsub -W 12:00 -q gpu_geforceGTX -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]" -n 1 /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/run_hiccups_cn_3.sh
bsub -W 12:00 -q gpu_geforceGTX -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]" -n 1 /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/run_hiccups_cn_4.sh

bsub -W 12:00 -q gpu_geforceGTX -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]" -n 1 /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/run_hiccups_npc_1.sh
bsub -W 12:00 -q gpu_geforceGTX -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]" -n 1 /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/run_hiccups_npc_2.sh
bsub -W 12:00 -q gpu_geforceGTX -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]" -n 1 /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/run_hiccups_npc_3.sh
bsub -W 12:00 -q gpu_geforceGTX -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]" -n 1 /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/run_hiccups_npc_4.sh


#### Sep25, 2018
# run merged .hic
cd /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/
bsub -W 12:00 -q gpu_geforceGTX -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]" -n 1 /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/run_hiccups_ch12.sh
bsub -W 12:00 -q gpu_geforceGTX -R "rusage[mem=10000,ngpus_excl_p=1] span[hosts=1]" -n 1 /project/umw_zhiping_weng/fankaili/for_Juicer/scripts/run_hiccups_heart.sh
