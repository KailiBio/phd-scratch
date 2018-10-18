#!/usr/bin/env python

# -- Kaili
# This script is for calculating mean signal values for each state for one mark.
# INPUT:
# OUTPUT: signal mean value for all states.
# EXP: python get_state_signal_mean_mark.py
#       /data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/run_IDEAS_8hm_atac_dname_pvalue.chr
#       /data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input
#       /data/zusers/fankaili/ideas/compare_results/66samples_10marks/66samples_10marks_state_signal_mean_H3K4me2.txt
#       H3K4me2 66

import os,sys

state_file = sys.argv[1]
inputFile = sys.argv[2]
outFile = sys.argv[3]
mark = sys.argv[4]
n = sys.argv[5]

# state_file="/data/zusers/fankaili/ideas/run_ideas_p_value/IDEAS_8hm_atac_dname_pvalue_result/run_IDEAS_8hm_atac_dname_pvalue.chr"
# inputFile="/data/zusers/fankaili/ideas/run_ideas_p_value/run_IDEAS_8hm_atac_dname_pvalue.input"
# outFile="/data/zusers/fankaili/ideas/compare_results/66samples_10marks/66samples_10marks_state_signal_mean_H3K4me2.txt"
# mark="H3K4me2"
# n=66

################
def read_input_into_Dic(input_file_path):
    '''
    This function for
    1. read input file into dic.
    2. get all the marks into list.
    '''
    input = {}
    mark_list = []
    for line in open(input_file_path).readlines():
        biosample, mark, file = line.rstrip().split(" ")
        file = file.replace("txt","tab")
        # save marks
        if mark not in mark_list:
            mark_list.append(mark)
        # save input into dic
        if input.has_key(biosample):
            input[biosample][mark] = file
        else:
            d={}
            d[mark] = file
            input[biosample] = d
    return input, mark_list


def read_signal_file_into_dic(signal_file, mark):
    '''
    read signal file into dic.
    id as key and signal_value as value.
    '''
    signal_dic = {}
    if mark=="DNAme":
        for line in open(signal_file).readlines():
            id = line.rstrip().split("\t")[0]
            value = float(line.rstrip().split("\t")[5])
            signal_dic[id] = value
    else:
        for line in open(signal_file).readlines():
            id = line.rstrip().split("\t")[0]
            value = float(line.rstrip().split("\t")[4])
            signal_dic[id] = value
    return signal_dic


def get_state_signal_sum_num(chr, state_file, input, mark, n):
    ''' get the sum of all signal values for given state, also count numbers. '''
    state_signal_sum = {}
    state_num = {}
    for i in chr:
        file = state_file + i + ".state"
        for j in range(4,n+4):
            for line in open(file).readlines():
                # get header (biosample for each column)
                if line.startswith("#ID"):
                    header = line.rstrip().split(" ")
                    biosample = header[j]
                    if input[biosample].has_key(mark):
                        signal_file = input[biosample][mark]
                        signal_dic = read_signal_file_into_dic(signal_file, mark)
                    else:
                        break
                # sum all signal values
                else:
                    l = line.rstrip().split(" ")
                    id = l[0]
                    state = l[n]
                    if state_signal_sum.has_key(state):
                        state_signal_sum[state]+=signal_dic[id]
                        state_num[state]+=1
                    else:
                        state_signal_sum[state] = signal_dic[id]
                        state_num[state] = 1
    #
    return state_signal_sum, state_num, biosample


# def print_sum_signal_num(outFile, dic, mark):
#     '''print matrix.'''
#     output = open(outFile, "w")
#     print >> output, "state\t"+mark
#     for state, value in dic.items():
#         print >> output, state+"\t"+str(value)
#     #
#     output.close()


def calculate_state_mean_signal(outFile, state_signal_sum, state_num, mark):
    '''calculate mean of signal, print out'''
    output = open(outFile, "w")
    print >> output, "state\t"+mark
    for state, value in state_signal_sum.items():
        mean = round(value/ state_num[state],2)
        print >> output, state+"\t"+str(mean)
    #
    output.close()


def main():
    # read input file
    input, mark_list = read_input_into_Dic(inputFile)
    # get chromosome number. For mouse, here is chr1-19&X&Y.
    chr = [str(x+1) for x in range(19)]
    chr.append("X")
    chr.append("Y")
    # get state sum signal value & num
    state_signal_sum, state_num, biosample = get_state_signal_sum_num(chr, state_file, input, mark, int(n))
    # calculate mean and print
    calculate_state_mean_signal(outFile, state_signal_sum, state_num, mark)


if __name__ == "__main__":
    main()
