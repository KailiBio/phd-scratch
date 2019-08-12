#!/bin/env python

# -- Kaili
# This script is for running Avocado to do CTCF imputation.

import itertools
import numpy as np;
np.random.seed(0)

import keras
from keras.layers import Input, Embedding, Dense
from keras.layers import Multiply, Dot, Flatten, concatenate
from keras.models import Model
from keras.optimizers import Adam

from avocado import Avocado

import sys

signal_path=sys.argv[1]
output_path=sys.argv[2]

###############
samples = ['forebrain_0', 'midbrain_0', 'hindbrain_0', 'heart_0', 'intestine_0', 'kidney_0', 'liver_0', 'liver_14.5', 'lung_0', 'lung_14.5', 'stomach_0']
marks = ['H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9me3', 'H3K9ac', 'H3K27me3', 'H3K27ac','H3K36me3', 'ATAC', 'CTCF', 'DNAme']

data = {}
for sample, mark in itertools.product(samples, marks):
    if sample == 'liver_14.5' and mark == 'CTCF':
        continue
    if sample == 'lung_14.5' and mark == 'CTCF':
        continue
    filename= (signal_path+'{}_{}_chr19_25bp.npz').format(sample, mark)
    data[(sample, mark)] = np.load(filename)['dat']

model = Avocado(celltypes=samples, assays=marks, n_genomic_positions=2328192)
model.fit(data)

track1 = model.predict("liver_14.5", "CTCF")
track2 = model.predict("lung_14.5", "CTCF")

np.savetxt(output_path+'avocado_predict_liver_14.5_CTCF.txt', track1, fmt='%.4f')
np.savetxt(output_path+'avocado_predict_lung_14.5_CTCF.txt', track2, fmt='%.4f')
