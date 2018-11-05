#!/usr/bin/env python

# source code from Jill
# modify for adding different column as signal-input.
# no log for z-score
# -- Kaili


import numpy, sys, math

sig=[[],[]]
masterPeak=[]
calculate=[]
bigWig=open(sys.argv[1])
n = int(sys.argv[2])-1    # column number. refer to bigWigAverageOverBed --help

for line in open(sys.argv[1]):
    line=line.rstrip().split("\t")
    if float(line[n]) == 0:
	sig[1].append("Zero")
	sig[0].append((float(line[n])))
	masterPeak.append(line[0])
    else:
    	sig[1].append(float(line[n]))
	sig[0].append((float(line[n])))
	calculate.append(float(line[n]))
	masterPeak.append(line[0])

lmean=numpy.mean(calculate)
lstd=numpy.std(calculate)
i=0
for entry in sig[1]:
    if entry != "Zero":
	print masterPeak[i], "\t", (entry-lmean)/lstd, "\t", sig[0][i], "\t", sig[1][i]
    else:
	print masterPeak[i], "\t", -10, "\t", 0, "\t", -10
    i+=1
bigWig.close()
