#!/usr/bin/env python

# -- Kaili
# This script is for calculating transition probability.
# for normal bins and DHS-bins.abs, output state transition count&probability.
# INPUT: state_file. (in signle chromatin, sorted by loci)
#              type of bins. (normal or dhs)
#              outfile_count. Output a matrix of state transition counts.
#              outfile_prob. Output a matrix of state transition probability.
# OUTPUT: a matrix of state transition counts & a matrix of state transition probability.
# EXP: python calculate_transition_probability.py \
#           "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/transition_rate/dhs_bins_chr10_intestine_16.5_state.txt" \
#           "dhs" "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/transition_rate/intestine_16.5_chr10_dhs_transitionCount.txt"
#           "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/transition_rate/intestine_16.5_chr10_dhs_transitionProb.txt"
import re, os, sys

def read_state_annotation(file):
    dic = {}
    for line in open(file).readlines():
        line = line.rstrip().split("\t")
        dic[line[0]] = line[2]
    return dic

def initiate_transition_dic(state_list, initiation_value):
    dic = {}
    for i in state_list:
        for j in state_list:
            dic[i+"2"+j] = initiation_value
    return dic

def count_state_transition(state_file, state_dic, state_list):
    count = initiate_transition_dic(state_list, 0)
    end = ""
    end_state = ""
    for line in open(state_file).readlines():
        line = line.rstrip().split("\t")
        if end == line[1]:
            state = state_dic[line[3]]
            transition = end_state + "2" + state
            count[transition] +=1
        else:
            end = line[2]
            if state_dic.has_key(line[3]):
                end_state = state_dic[line[3]]
            else:
                end_state=""
    return count

def count_dic_2_matrix(state_transition_count_dic, state_list, outfile_count, outfile_prob):
    output_count = open(outfile_count, "w+")
    output_prob = open(outfile_prob, "w+")
    print >> output_count, ("\t").join(state_list)
    print >> output_prob, ("\t").join(state_list)
    #
    for i in state_list:
        outline = [i]
        num = list()
        for j in state_list:
            outline.append(str(state_transition_count_dic[i+"2"+j]))
            num.append(state_transition_count_dic[i+"2"+j])
        print >> output_count, ("\t").join(outline)
        total = float(sum(num))
        prob = [str(round(c/total,2)) for c in num]
        print >> output_prob, ("\t").join([i]+prob)
    output_count.close()
    output_prob.close()


if __name__ == "__main__":
    normal_state_anno = read_state_annotation("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/normal_bins_state_annotation.txt")
    dhs_state_anno = read_state_annotation("/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/dhs_bins_state_annotation.txt")
    #
    state_file = sys.argv[1]
    type_of_bins = sys.argv[2]
    outfile_count = sys.argv[3]
    outfile_prob = sys.argv[4]
    # choose the right state annotation file
    if type_of_bins=="normal":
        state_dic = normal_state_anno
    else:
        state_dic = dhs_state_anno
    # count state transition
    state_list = ["P", "O", "E", "B", "T", "Q", "H", "R"]
    state_transition_count_dic = count_state_transition(state_file, state_dic, state_list)
    # transform to transition count matrix
    count_dic_2_matrix(state_transition_count_dic, state_list, outfile_count, outfile_prob)
