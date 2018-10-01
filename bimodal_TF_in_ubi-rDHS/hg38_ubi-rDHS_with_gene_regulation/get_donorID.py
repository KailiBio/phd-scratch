#!/usr/bin/env python

# -- Kaili
# This script is getting all the donor ID of the file from given expID and fileID.
# INPUT: inputFile: file that eacho line is expID, fileID and sample. (eg: ENCSR000AEC	ENCFF906LSJ	GM12878_NA)
#        outputFile
# OUTPUT:file that each line is expID, fileID, sample, file bio_rep and donorIDs.

# EXP: python get_donorID.py "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_closest_gene_exp_list.txt" "/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_closest_gene_exp_list_donor.txt"


import re, os, sys
import subprocess

sys.path.append('/home/fankaili/git/metadata/utils')
from querydcc import QueryDCC
from exp import Exp
from exp_file import ExpFile

def write_file(file, output_path):
	o= open(output_path, 'w')
	o.writelines(file)
	o.close()

inputFile = sys.argv[1]
outputFile = sys.argv[2]

# inputFile="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_closest_gene_exp_list.txt"
# outputFile="/data/zusers/fankaili/ccre/hg38_ubi-rDHS/hg38_closest_gene_exp_list_donor.txt"

##################################
if __name__ == "__main__":
    output = open(outputFile, "w")

    for line in open(inputFile).readlines():
        expID = line.rstrip().split("\t")[0]
        fileID = line.rstrip().split("\t")[1]
        sample = line.rstrip().split("\t")[2]

        myexp = myexp = Exp.fromJsonFile(expID)

        # get donorID
        if len(myexp.jsondata["replicates"])==2 :
            if myexp.jsondata["replicates"][0]["biological_replicate_number"]==1:
                donor1 = myexp.jsondata["replicates"][0]["library"]["biosample"]["accession"]
                donor2 = myexp.jsondata["replicates"][1]["library"]["biosample"]["accession"]
            else:
                donor1 = myexp.jsondata["replicates"][1]["library"]["biosample"]["accession"]
                donor2 = myexp.jsondata["replicates"][0]["library"]["biosample"]["accession"]
        elif len(myexp.jsondata["replicates"])==1 :
            donor1 = myexp.jsondata["replicates"][0]["library"]["biosample"]["accession"]
            donor2 = "---"
        else:
            for i in range(0,len(myexp.jsondata["replicates"])):
		if myexp.jsondata["replicates"][i]["biological_replicate_number"]==1:
		    donor1 = myexp.jsondata["replicates"][i]["library"]["biosample"]["accession"]
	        else:
                    donor2 += myexp.jsondata["replicates"][i]["library"]["biosample"]["accession"]+"_"

        # get file bio_replicate number
        for i in range(0,len(myexp.files)):
            if myexp.files[i].accession == fileID:
                bio_rep = str(myexp.files[i].bio_rep[0])
                break

        print >> output, "\t".join([expID, fileID, sample, bio_rep, donor1, donor2])

    output.close()
