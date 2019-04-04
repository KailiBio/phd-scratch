#!/usr/bin/env python

# -- Kaili
# This script is for getting ENCODE expID and info from given url.

import re, os, sys
import subprocess

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

##################################
if __name__ == "__main__":
    url = sys.argv[1]
    outDir = sys.argv[2]

    q = QueryDCC(auth=False)
    expID = []
    out = open(outDir, "w")
    for exp in q.getExps(url):
        expID = exp.encodeID
        assay = exp.assay_term_name
        mark = exp.label
        biosample = ("_").join(exp.biosample_summary.rstrip().split(" ")).replace('(','').replace(')','')
        print >> out, ("\t").join([expID, assay, mark,biosample]).encode('utf-8').strip()

    out.close()
