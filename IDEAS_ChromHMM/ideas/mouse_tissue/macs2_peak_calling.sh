#!/bin/bash

# -- Kaili
# This script is for using only rep2 to call peaks by MACS2.

bam="/data/zusers/fankaili/ideas/CTCF_impute/ctcf_samples/ENCFF676SED.bam"

mkdir macs2

macs2 callpeak -t TC1-H3K4-ST2-D0.GRCm38.p3.q30.bam \
               -c TC1-I-ST2-D0.GRCm38.p3.q30.bam \
               -f BAM -g mm -n TC1-ST2-H3K4-D0 -B -q 0.01 --outdir ./macs2
