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

#########
def build_model(n_celltypes, n_celltype_factors, n_assays, n_assay_factors,
    n_genomic_positions, n_25bp_factors, n_250bp_factors, n_5kbp_factors,
    n_layers, n_nodes, freeze_celltypes=False, freeze_assays=False,
    freeze_genome_25bp=False, freeze_genome_250bp=False,
    freeze_genome_5kbp=False, freeze_network=False):
    """This function builds a multi-scale deep tensor factorization model."""

    celltype_input = Input(shape=(1,), name="celltype_input")
    celltype_embedding = Embedding(n_celltypes, n_celltype_factors,
        input_length=1, name="celltype_embedding")
    celltype_embedding.trainable = not freeze_celltypes
    celltype = Flatten()(celltype_embedding(celltype_input))

    assay_input = Input(shape=(1,), name="assay_input")
    assay_embedding = Embedding(n_assays, n_assay_factors,
        input_length=1, name="assay_embedding")
    assay_embedding.trainable = not freeze_assays
    assay = Flatten()(assay_embedding(assay_input))

    genome_25bp_input = Input(shape=(1,), name="genome_25bp_input")
    genome_25bp_embedding = Embedding(n_genomic_positions, n_25bp_factors,
        input_length=1, name="genome_25bp_embedding")
    genome_25bp_embedding.trainable = not freeze_genome_25bp
    genome_25bp = Flatten()(genome_25bp_embedding(genome_25bp_input))

    genome_250bp_input = Input(shape=(1,), name="genome_250bp_input")
    genome_250bp_embedding = Embedding((n_genomic_positions / 10) + 1,
        n_250bp_factors, input_length=1, name="genome_250bp_embedding")
    genome_250bp_embedding.trainable = not freeze_genome_250bp
    genome_250bp = Flatten()(genome_250bp_embedding(genome_250bp_input))

    genome_5kbp_input = Input(shape=(1,), name="genome_5kbp_input")
    genome_5kbp_embedding = Embedding((n_genomic_positions / 200) + 1,
        n_5kbp_factors, input_length=1, name="genome_5kbp_embedding")
    genome_5kbp_embedding.trainable = not freeze_genome_5kbp
    genome_5kbp = Flatten()(genome_5kbp_embedding(genome_5kbp_input))

    layers = [celltype, assay, genome_25bp, genome_250bp, genome_5kbp]
    inputs = (celltype_input, assay_input, genome_25bp_input,
        genome_250bp_input, genome_5kbp_input)

    x = concatenate(layers)
    for i in range(n_layers):
        layer = Dense(n_nodes, activation='relu', name="dense_{}".format(i))
        layer.trainable = not freeze_network
        x = layer(x)

    layer = Dense(1, name="y_pred")
    layer.trainable = not freeze_network
    y = layer(x)

    model = Model(inputs=inputs, outputs=y)
    model.compile(optimizer='adam', loss='mse', metrics=['mse'])
    return model

###############
samples = ['forebrain_0', 'midbrain_0', 'hindbrain_0', 'heart_0', 'intestine_0', 'kidney_0', 'liver_0', 'liver_14.5', 'lung_0', 'lung_14.5', 'stomach_0']
marks = ['H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9me3', 'H3K9ac', 'H3K27me3', 'H3K27ac','H3K36me3', 'ATAC', 'CTCF', 'DNAme']

data = {}
for sample, mark in itertools.product(samples, marks):
    if sample == 'liver_14.5' and mark == 'CTCF':
        continue
    if sample == 'lung_14.5' and mark == 'CTCF':
        continue
    filename= '/data/zusers/fankaili/ideas/signal/rep1_signal_chr19_25bp/{}_{}_chr19_25bp.npz'.format(sample, mark)
    data[(sample, mark)] = np.load(filename)['dat']

model = Avocado(celltypes=samples, assays=marks, n_genomic_positions=2328192)
model.fit(data)

track1 = model.predict("liver_14.5", "CTCF")
track2 = model.predict("lung_14.5", "CTCF")

np.savetxt('/data/zusers/fankaili/ideas/imputation_comparison/avocado/avocado_predict_liver_14.5_CTCF.txt', track1, fmt='%.4f')
np.savetxt('/data/zusers/fankaili/ideas/imputation_comparison/avocado/avocado_predict_lung_14.5_CTCF.txt', track2, fmt='%.4f')
