#!/usr/bin/env python

# -- Kaili
# This script is for matching all data.

import re, os, sys
import subprocess

## function
def read_CTCF_file(filename):
    dic = {}
    for line in open(filename).readlines():
        id = line.rstrip().split("\t")[0]
        biosample = line.rstrip().split("\t")[3]
        if dic.has_key(biosample):
            new_id = dic[biosample] + ";" + id
            dic[biosample] = new_id
        else:
            dic[biosample] = id
    return dic

def read_OCR_file(filename):
    dic={}
    dic["DNase-seq"] = {}
    dic["ATAC-seq"] = {}
    for line in open(filename).readlines():
        id = line.rstrip().split("\t")[0]
        biosample = line.rstrip().split("\t")[3]
        assay = line.rstrip().split("\t")[1]
        if assay=="ATAC-seq":
            if dic["ATAC-seq"].has_key(biosample):
                new_id = dic["ATAC-seq"][biosample] + ";" +id
                dic["ATAC-seq"][biosample] = new_id
            else:
                dic["ATAC-seq"][biosample] = id
        else:
            if dic["DNase-seq"].has_key(biosample):
                new_id = dic["DNase-seq"][biosample] + ";" +id
                dic["DNase-seq"][biosample] = new_id
            else:
                dic["DNase-seq"][biosample] = id
    return dic

def read_3D_file(filename):
    dic={}
    dic["HiC"] = {}
    dic["ChIA-PET"] = {}
    for line in open(filename).readlines():
        id = line.rstrip().split("\t")[0]
        biosample = line.rstrip().split("\t")[3]
        assay = line.rstrip().split("\t")[1]
        if assay=="HiC":
            if dic["HiC"].has_key(biosample):
                new_id = dic["HiC"][biosample] + ";" +id
                dic["HiC"][biosample] = new_id
            else:
                dic["HiC"][biosample] = id
        else:
            if dic["ChIA-PET"].has_key(biosample):
                new_id = dic["ChIA-PET"][biosample] + ";" +id
                dic["ChIA-PET"][biosample] = new_id
            else:
                dic["ChIA-PET"][biosample] = id
    return dic

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

def get_id_from_dic(dic, key):
    if dic.has_key(key):
        out = dic[key]
    else:
        out = "---"
    return out

if __name__ == "__main__":
    # human
    # data
    encode_CTCF_file = "/data/zusers/fankaili/CTCF_and_OCR/encode_human_CTCF_explist.txt"
    encode_OCR_file = "/data/zusers/fankaili/CTCF_and_OCR/encode_human_OCR_explist.txt"
    cistrome_CTCF_file = "/data/zusers/fankaili/CTCF_and_OCR/Cistrome_human_CTCF_datalist.txt"
    cistrome_OCR_file = "/data/zusers/fankaili/CTCF_and_OCR/Cistrome_human_OCR_datalist.txt"
    encode_3D_file = "/data/zusers/fankaili/CTCF_and_OCR/encode_human_3D_explist.txt"
    #
    outDir = "human_matched_data_masterlist.txt"
    ## read data
    encode_CTCF = read_CTCF_file(encode_CTCF_file)
    encode_OCR = read_OCR_file(encode_OCR_file)
    cistrome_CTCF = read_CTCF_file(cistrome_CTCF_file)
    cistrome_OCR = read_OCR_file(cistrome_OCR_file)
    encode_3D = read_3D_file(encode_3D_file)
    ## find match, output
    biosample_list = union(encode_CTCF.keys(), cistrome_CTCF.keys())
    out = open(outDir, "w")
    print >> out, ("\t").join(["biosample", "ENCODE_CTCF", "ENCODE_DNase", "ENCODE_ATAC", "Cistrome_CTCF", "Cistrome_DNase", "Cistrome_ATAC", "ENCODE_Hi-C", "ENCODE_ChIA-PET"])
    for biosample in biosample_list:
        encode_ctcf_id = get_id_from_dic(encode_CTCF, biosample)
        encode_dnase_id = get_id_from_dic(encode_OCR["DNase-seq"], biosample)
        encode_atac_id = get_id_from_dic(encode_OCR["ATAC-seq"], biosample)
        cistrome_ctcf_id = get_id_from_dic(cistrome_CTCF, biosample)
        cistrome_dnase_id = get_id_from_dic(cistrome_OCR["DNase-seq"], biosample)
        cistrome_atac_id = get_id_from_dic(cistrome_OCR["ATAC-seq"], biosample)
        encode_hic_id = get_id_from_dic(encode_3D["HiC"], biosample)
        encode_cp_id = get_id_from_dic(encode_3D["ChIA-PET"], biosample)
        #
        print >> out, ("\t").join([biosample, encode_ctcf_id, encode_dnase_id, encode_atac_id, cistrome_ctcf_id, cistrome_dnase_id, cistrome_atac_id, encode_hic_id, encode_cp_id])

    out.close()

    # mouse
    # data
    encode_CTCF_file = "/data/zusers/fankaili/CTCF_and_OCR/encode_mouse_CTCF_explist.txt"
    encode_OCR_file = "/data/zusers/fankaili/CTCF_and_OCR/encode_mouse_OCR_explist.txt"
    cistrome_CTCF_file = "/data/zusers/fankaili/CTCF_and_OCR/Cistrome_mouse_CTCF_datalist.txt"
    cistrome_OCR_file = "/data/zusers/fankaili/CTCF_and_OCR/Cistrome_mouse_OCR_datalist.txt"
    encode_3D_file = "/data/zusers/fankaili/CTCF_and_OCR/encode_mouse_3D_explist2.txt"
    #
    outDir = "mouse_matched_data_masterlist.txt"
    ## read data
    encode_CTCF = read_CTCF_file(encode_CTCF_file)
    encode_OCR = read_OCR_file(encode_OCR_file)
    cistrome_CTCF = read_CTCF_file(cistrome_CTCF_file)
    cistrome_OCR = read_OCR_file(cistrome_OCR_file)
    encode_3D = read_3D_file(encode_3D_file)
    ## find match, output
    biosample_list = union(encode_CTCF.keys(), cistrome_CTCF.keys())
    out = open(outDir, "w")
    print >> out, ("\t").join(["biosample", "ENCODE_CTCF", "ENCODE_DNase", "ENCODE_ATAC", "Cistrome_CTCF", "Cistrome_DNase", "Cistrome_ATAC", "ENCODE_Hi-C", "ENCODE_ChIA-PET"])
    for biosample in biosample_list:
        encode_ctcf_id = get_id_from_dic(encode_CTCF, biosample)
        encode_dnase_id = get_id_from_dic(encode_OCR["DNase-seq"], biosample)
        encode_atac_id = get_id_from_dic(encode_OCR["ATAC-seq"], biosample)
        cistrome_ctcf_id = get_id_from_dic(cistrome_CTCF, biosample)
        cistrome_dnase_id = get_id_from_dic(cistrome_OCR["DNase-seq"], biosample)
        cistrome_atac_id = get_id_from_dic(cistrome_OCR["ATAC-seq"], biosample)
        encode_hic_id = get_id_from_dic(encode_3D["HiC"], biosample)
        encode_cp_id = get_id_from_dic(encode_3D["ChIA-PET"], biosample)
        #
        print >> out, ("\t").join([biosample, encode_ctcf_id, encode_dnase_id, encode_atac_id, cistrome_ctcf_id, cistrome_dnase_id, cistrome_atac_id, encode_hic_id, encode_cp_id])

    out.close()
