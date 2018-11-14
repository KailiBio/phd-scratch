#!/usr/bin/env python

# -- Kaili
# This script is for working on the small gaps between DHSs.
### can be well packaged later!!!

import re, os, sys
import subprocess

ocr = "/data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs.bed_withLength.txt"
outfile = "/data/zusers/fankaili/ideas/dhs_bins/mm10-rOCRs_v3.txt"

dhsDic = {}
for line in open(ocr).readlines():
    id = line.rstrip().split("\t")[4]
    dhsDic[id] = line.rstrip()

for line in open(ocr).readlines():
    gap_length = int(line.rstrip().split("\t")[6])
    if gap_length <= 100:
        id = line.rstrip().split("\t")[4]
        id_0 = str(int(id)-1)
        #
        l0 = gap_length/2
        l1 = gap_length - l0
        #
        line0 = dhsDic[id_0]
        items_0 = line0.split("\t")
        new_0 = str(int(items_0[2])+l0)
        line0_new = ("\t").join([items_0[0], items_0[1], new_0, items_0[3], items_0[4], items_0[5], items_0[6]])
        dhsDic[id_0] = line0_new
        #
        line1 = dhsDic[id]
        items_1 = line1.split("\t")
        new_1 = str(int(items_1[1])-l1)
        line1_new = ("\t").join([items_1[0], new_1, items_1[2], items_1[3], items_1[4], items_1[5], items_1[6]])
        dhsDic[id] = line1_new

output = open(outfile, "w+")
for line in dhsDic.values():
    print >> output, line

output.close()
