#!/usr/bin/env python

# -- Kaili
# This script is for geting data matrix for calculating r-square across cell type.
# INPUT: prefix (dhs or normal)
#               gene list file(only one column for geneID)
#               work directory
# OUTPUT: matrix of adjusted r-squared
# EXP: python get_matrix_for_regression_across_celltype.py dhs 0-0.5
#           /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/across_sample_matrix/

import re, os, sys
import subprocess
import numpy as np

def get_list(file):
    '''
    This function is for reading file into list.
    '''
    list = []
    for line in open(file).readlines():
        line = line.rstrip()
        list.append(line)
    return list

def get_exp_dic(expression_file):
    '''
    This function is for reading expression matix into dictionary.
    First key is geneID, second is sample.
    '''
    read_exp = np.genfromtxt(expression_file, dtype=None, encoding=None)
    header = read_exp[0]
    exp = read_exp[1:]
    #
    dic = {}
    for i in range(exp.shape[0]):
        gene = exp[i][0]
        dic[gene] = {}
        for j in range(1,exp.shape[1]):
            dic[gene][header[j]] = exp[i][j]
    return(dic)

def get_gene_state_proportion_dic(window, celltype_list, prefix):
    '''
    This function is for reading all the state proportion file into disctionary.
    First key is window number, second is geneID, third is sample.
    '''
    out_dic = {}
    #
    sample = celltype_list[0]
    state_file = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/" + prefix + "_state_proportion/mm10_gene_" + sample + "_" + prefix + "_state_count_window_" + str(window) + ".txt"
    for line in open(state_file).readlines():
        line = line.rstrip().split("\t")
        gene = line[0]
        states = line[2:]
        out_dic[gene] = {}
        out_dic[gene][sample] = ("\t").join(states)
    #
    for sample in celltype_list[1:]:
        state_file = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/" + prefix + "_state_proportion/mm10_gene_" + sample + "_" + prefix + "_state_count_window_" + str(window) + ".txt"
        for line in open(state_file).readlines():
            line = line.rstrip().split("\t")
            gene = line[0]
            states = line[2:]
            out_dic[gene][sample] = ("\t").join(states)
    #
    return out_dic

if __name__ == "__main__":
    prefix = sys.argv[1]
    sd_range = sys.argv[2]
    workDir = sys.argv[3]

    # get files:
    ## get celltype list
    celltype_file = "/data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt"
    celltype_list = get_list(celltype_file)
    ## get expression
    expression_file = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/mm10_RNA_protein-coding_tpm_matrix_matched.txt"
    exp = get_exp_dic(expression_file)
    ## get gene list
    gene_file = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list_sd_"+sd_range+".txt"
    gene_list = get_list(gene_file)
    ## get state proportion
    ### for each window, put all cell type state proportion into dic
    matrix = {}
    for i in range(1,21):
        # read state proportion file, put into dic make gene as key
        matrix[i] = get_gene_state_proportion_dic(i, celltype_list, prefix)

    # calculate R square for each gene
    for i in range(1,21):
        outfile = open(workDir+"regression_across_celltype_data_"+prefix+"_"+sd_range+"_window"+str(i)+".txt", "w+")
        #
        for gene in gene_list:
            if exp.has_key(gene) and matrix[i].has_key(gene):
                for sample in celltype_list:
                    y = exp[gene][sample]
                    x = matrix[i][gene][sample]
                    print >> outfile, gene+"\t"+sample+"\t"+y+"\t"+x
        outfile.close()
