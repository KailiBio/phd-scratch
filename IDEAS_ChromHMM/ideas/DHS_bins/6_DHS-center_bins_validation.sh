#!/bin/bash

# -- Kaili
# This script is for validating DHS-bins results.
# 1. aggregation plot
# 2. states transition rate

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/"
workDir="/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/"

cd ${workDir}

# 1. aggregation plot
bash ${scriptDir}make_DHSbins_marks_aggregation_plot.sh

# 2. states transition rate
bash ${scriptDir}calculate_state_transition_rate.sh
