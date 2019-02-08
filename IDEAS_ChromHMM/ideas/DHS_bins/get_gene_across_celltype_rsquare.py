#!/usr/bin/env python

# -- Kaili
# This script is for geting r-squre for all genes across cell type.

import re, os, sys
from subprocess import Popen, PIPE, STDOUT

def get_gene_list(gene_file):
    gene_list = []
    for line in open(gene_file).readlines():
        line = line.rstrip().split("\t")
        gene_list.append(line[3])
    return gene_list

def get_gene_state_proportion_dic(state_file):
    out_dic = {}
    for line in open(state_file).readlines():
        line = line.rstrip().split("\t")
        gene = line[0]
        states = line[2:]
        out_dic[gene] = states
    return out_dic

if __name__ == "__main__":
    prefix = sys.argv[1]
    workDir = sys.argv[2]
    #
    gene_file = "/home/fankaili/genome/mm10_vM4_protein_coding.bed"
    celltype_file = "/data/zusers/fankaili/ideas/ENCODE_mouse_rep1_sample_list.txt"
    # get gene list
    gene_list = get_gene_list(gene_file)
    # for each window, put all cell type state proportion into dic
    matrix = {}
    for i in range(1,21):
        matrix[i] = {}
        for sample in open(celltype_file).readlines():
            sample = sample.rstrip()
            # read state proportion file, put into dic make gene as key
            state_file = "/data/zusers/fankaili/ideas/dhs_bins/v3_100_400bp/gene_expression/"+prefix+"_state_proportion/mm10_gene_"+sample+"_"+prefix+"_state_count_window_"+str(i)+".txt"
            matrix[i][sample] = get_gene_state_proportion_dic(state_file)
    # calculate R square for each gene
    outfile = open(workDir+"gene_rsquare_matrix_"+prefix+".txt", "w+")
    for gene in gene_list:
        outlist = [gene]
        for i in range(1,21):
            # get tmp file for gene in window i
            output = open("tmp_rsquare_"+prefix+".txt", "w+")
            for sample in open(celltype_file).readlines():
                sample = sample.rstrip()
                print >> output, gene+"\t"+("\t").join(matrix[i][sample][gene])+"\n"
            output.close()
            # calculate r-square
            cmd = "Rscript /data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/DHS_bins/calculate_rsqaure.R tmp_rsquare_"+prefix+".txt "+prefix+" | awk '{print $2}'"
            p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            outlist.append(p.stdout.read()[:-2])
        print >> outfile, ("\t").join(outlist) +"\n"
    outfile.close()
