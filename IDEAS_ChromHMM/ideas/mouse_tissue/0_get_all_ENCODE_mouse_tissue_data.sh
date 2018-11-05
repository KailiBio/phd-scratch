#!/bin/bash

# -- Kaili
# This is script for getting all the 870 mouse tissue file ID from ENCODE.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/"

python ${scriptDir}get_mouse_data_from_ENCODE_json.py
