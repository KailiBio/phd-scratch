#!/bin/env python

# -- Kaili
# This script is for converting all signal file into compressed numpy file.

import numpy as np
import itertools

samples = ['forebrain_0', 'midbrain_0', 'hindbrain_0', 'heart_0', 'intestine_0', 'kidney_0', 'liver_0', 'liver_14.5', 'lung_0', 'lung_14.5', 'stomach_0']
marks = ['H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9me3', 'H3K9ac', 'H3K27me3', 'H3K27ac','H3K36me3', 'ATAC', 'CTCF', 'DNAme']

for sample, mark in itertools.product(samples, marks):
    inputfile = '/data/zusers/fankaili/ideas/signal/rep1_signal_chr19_25bp/{}_{}_chr19_25bp.txt'.format(sample, mark)
    outputfile = '/data/zusers/fankaili/ideas/signal/rep1_signal_chr19_25bp/{}_{}_chr19_25bp'.format(sample, mark)

    dat = np.loadtxt(inputfile)
    np.savez_compressed(outputfile, dat=dat)
