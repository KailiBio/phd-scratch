#!/usr/bin/env python

# Kaili
# This script is for getting cell-type specific ccRE definition file list.
# INPUT: master file from Jill: /data/projects/screen/Version-4/ver10/hg19/hg19-Look-Up-Matrix.txt
# OUTPUT: cell type specific file list: /data/zusers/fankaili/ccre/tf/cell_type_specific/hg19_ccREs_cell_type_specific_definition_file_list.txt

def write_file(file, output_path):
	o= open(output_path, 'w')
	o.writelines(file)
	o.close()

if __name__ == "__main__":
	out = []
	file=open("/data/projects/screen/Version-4/ver10/hg19/hg19-Look-Up-Matrix.txt").readlines()
	for line in file:
		cell_line = line.rstrip().split("\t")[0]
		f1 = line.rstrip().split("\t")[2]
		f2 = line.rstrip().split("\t")[3]
		f3 = line.rstrip().split("\t")[4]
		f4 = line.rstrip().split("\t")[5]
		file_left = []
		if (f1!="NA"):
			file_left.append(f1)
		if (f2!="NA"):
			file_left.append(f2)
		if (f3!="NA"):
			file_left.append(f3)
		if (f4!="NA"):
			file_left.append(f4)
		filename = "_".join(file_left)+".cREs.bed.gz"
		out.append(cell_line+"\t"+filename+"\n")
	write_file(out, "/data/zusers/fankaili/ccre/tf/cell_type_specific/hg19_ccREs_cell_type_specific_definition_file_list.txt")
