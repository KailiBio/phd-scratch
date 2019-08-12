#!/bin/env python

# -- Kaili
# This script is for converting given file into compressed numpy file.

import numpy as np
import itertools
import sys



inputfile = sys.argv[1]
outputfile = sys.argv[2]

dat = np.loadtxt(inputfile)
np.savez_compressed(outputfile, dat=dat)
