#!/usr/bin/env python

# -- Kaili
# This script is for geting data matrix for calculating r-square across cell type.
# INPUT: prefix (dhs or normal)
#               gene list file(only one column for geneID)
#               work directory
# OUTPUT: matrix of adjusted r-squared
# EXP: python get_gene_across_celltype_rsquare.py dhs
#           /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list.txt
#           /data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/

import re, os, sys
import subprocess


def get_list(file):
    list = []
    for line in open(file).readlines():
        line = line.rstrip()
        list.append(line)
    return list

def get_exp_dic(expression_file):
    dic = {}
    for line in open(expression_file).readlines()[1:]

def get_gene_state_proportion_dic(window, celltype_list, prefix):
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
    num = sys.argv[4]

    # get celltype list
    celltype_file = "/data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt"
    celltype_list = get_list(celltype_file)
    # get expression
    expression_file = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/mm10_RNA_protein-coding_tpm_matrix_matched.txt"
    # get gene list
    gene_file = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/gene_list_sd_"+sd_range+".txt"
    gene_list = get_list(gene_file)

    # for each window, put all cell type state proportion into dic
    matrix = {}
    for i in range(1,21):
        # read state proportion file, put into dic make gene as key
        matrix[i] = get_gene_state_proportion_dic(i, celltype_list, prefix)

    # calculate R square for each gene
    outfile = open(workDir+"regression_across_celltype_data_"+prefix+"_"+sd_range+".txt", "w+")
    for gene in gene_list:







    n=0
    for gene in gene_list:
        n +=1
        print(n)
        outlist = [gene]
        p = Pool(5)
        r2_1 = p.map(partial(calculate_rsquare, workDir = workDir, gene = gene, prefix = prefix), range(1,6))
        r2_2 = p.map(partial(calculate_rsquare, workDir = workDir, gene = gene, prefix = prefix), range(6,11))
        r2_3 = p.map(partial(calculate_rsquare, workDir = workDir, gene = gene, prefix = prefix), range(11,16))
        r2_4 = p.map(partial(calculate_rsquare, workDir = workDir, gene = gene, prefix = prefix), range(16,21))
        print >> outfile, gene + "\t" + ("\t").join(r2_1) + "\t" + ("\t").join(r2_2) + "\t" + ("\t").join(r2_3) + "\t" + ("\t").join(r2_4)
        subprocess.call("rm "+workDir+"tmp_rsquare_"+gene+"_"+prefix+"*.txt", shell=True)
    outfile.close()

matrix = {}
for i in range(1,21):
    # read state proportion file, put into dic make gene as key
    matrix[i] = get_gene_state_proportion_dic(i, celltype_list, prefix)
