#!/usr/bin/env python2.7

## Jul24.get_data_from_json.py

import json
import os,sys

def write_file(file, output_path):
	o= open(output_path, 'w')
	o.writelines(file)
	o.close()

jfile = json.load(open("/data/public_html_users/fankaili/IDEAS/allFiles.json"))
data = jfile["input data"]

path = "/data/zusers/fankaili/ideas/run_ideas/"

code = []
input = []
for i in range(len(data)):
    if data[i]['assay_term_name']=='whole-genome shotgun bisulfite sequencing':
        if (data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube"):
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_DNAme"
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test1_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group1.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test2_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group2.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test3_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group3.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test4_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group4.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test5_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group5.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test6_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group6.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test7_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group7.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test8_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group8.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test9_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group9.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test10_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group10.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test11_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group11.tab"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group1.tab > "+path+"signal/"+signal_file_name+"_200bin_group1.txt"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group2.tab > "+path+"signal/"+signal_file_name+"_200bin_group2.txt"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group3.tab > "+path+"signal/"+signal_file_name+"_200bin_group3.txt"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group4.tab > "+path+"signal/"+signal_file_name+"_200bin_group4.txt"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group5.tab > "+path+"signal/"+signal_file_name+"_200bin_group5.txt"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group6.tab > "+path+"signal/"+signal_file_name+"_200bin_group6.txt"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group7.tab > "+path+"signal/"+signal_file_name+"_200bin_group7.txt"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group8.tab > "+path+"signal/"+signal_file_name+"_200bin_group8.txt"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group9.tab > "+path+"signal/"+signal_file_name+"_200bin_group9.txt"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group10.tab > "+path+"signal/"+signal_file_name+"_200bin_group10.txt"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin_group11.tab > "+path+"signal/"+signal_file_name+"_200bin_group11.txt"+"\n")
        input.append(biosample+" DNAme "+ path + "signal/" + signal_file_name + "_200bin_group1.txt"+"\n")
    elif data[i]['assay_term_name']=='ATAC-seq':
        if data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube":
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_ATAC"
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test1_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group1.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test2_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group2.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test3_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group3.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test4_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group4.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test5_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group5.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test6_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group6.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test7_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group7.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test8_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group8.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test9_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group9.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test10_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group10.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test11_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group11.tab"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group1.tab > "+path+"signal/"+signal_file_name+"_200bin_group1.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group2.tab > "+path+"signal/"+signal_file_name+"_200bin_group2.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group3.tab > "+path+"signal/"+signal_file_name+"_200bin_group3.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group4.tab > "+path+"signal/"+signal_file_name+"_200bin_group4.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group5.tab > "+path+"signal/"+signal_file_name+"_200bin_group5.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group6.tab > "+path+"signal/"+signal_file_name+"_200bin_group6.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group7.tab > "+path+"signal/"+signal_file_name+"_200bin_group7.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group8.tab > "+path+"signal/"+signal_file_name+"_200bin_group8.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group9.tab > "+path+"signal/"+signal_file_name+"_200bin_group9.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group10.tab > "+path+"signal/"+signal_file_name+"_200bin_group10.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group11.tab > "+path+"signal/"+signal_file_name+"_200bin_group11.txt"+"\n")
        input.append(biosample+" ATAC "+ path + "signal/" + signal_file_name + "_200bin_group1.txt"+"\n")
    elif data[i]['assay_term_name']=='ChIP-seq':
        if data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube":
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_" + data[i]['label']
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test1_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group1.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test2_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group2.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test3_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group3.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test4_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group4.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test5_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group5.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test6_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group6.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test7_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group7.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test8_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group8.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test9_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group9.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test10_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group10.tab"+"\n")
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"tab_bed6/mm10_200bin_test11_tab.bed "+path+"signal/"+signal_file_name+"_200bin_group11.tab"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group1.tab > "+path+"signal/"+signal_file_name+"_200bin_group1.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group2.tab > "+path+"signal/"+signal_file_name+"_200bin_group2.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group3.tab > "+path+"signal/"+signal_file_name+"_200bin_group3.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group4.tab > "+path+"signal/"+signal_file_name+"_200bin_group4.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group5.tab > "+path+"signal/"+signal_file_name+"_200bin_group5.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group6.tab > "+path+"signal/"+signal_file_name+"_200bin_group6.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group7.tab > "+path+"signal/"+signal_file_name+"_200bin_group7.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group8.tab > "+path+"signal/"+signal_file_name+"_200bin_group8.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group9.tab > "+path+"signal/"+signal_file_name+"_200bin_group9.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group10.tab > "+path+"signal/"+signal_file_name+"_200bin_group10.txt"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin_group11.tab > "+path+"signal/"+signal_file_name+"_200bin_group11.txt"+"\n")
        input.append(biosample+" "+data[i]['label']+" "+ path + "signal/" + signal_file_name + "_200bin_group1.txt"+"\n")


write_file(code, path+"Jul24.get_signal_code.sh")
write_file(input, path+"IDEAS_group1.input")




code = []
input = []
for i in range(len(data)):
    if data[i]['assay_term_name']=='whole-genome shotgun bisulfite sequencing':
        if (data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube"):
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_DNAme"
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"mm10_tab.bed "+path+"signal/"+signal_file_name+"_200bin.tab"+"\n")
        code.append("""awk '{print $6}' """+path+"signal/"+signal_file_name+"_200bin.tab > "+path+"signal/"+signal_file_name+"_200bin.txt"+"\n")
        input.append(biosample+" DNAme "+ path + "signal/" + signal_file_name + "_200bin.txt"+"\n")
    elif data[i]['assay_term_name']=='ATAC-seq':
        if data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube":
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_ATAC"
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"mm10tab.bed "+path+"signal/"+signal_file_name+"_200bin.tab"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin.tab > "+path+"signal/"+signal_file_name+"_200bin.txt"+"\n")
        input.append(biosample+" ATAC "+ path + "signal/" + signal_file_name + "_200bin.txt"+"\n")
    elif data[i]['assay_term_name']=='ChIP-seq':
        if data[i]['biosample_term_name']=="embryonic facial prominence" or data[i]['biosample_term_name']=="neural tube":
            a = data[i]['biosample_term_name'].split(" ")
            biosample = "-".join(a) + "_" + data[i]['age']
        else:
            biosample = data[i]['biosample_term_name'] + "_" + data[i]['age']
        signal_file_name = biosample + "_" + data[i]['label']
        code.append("bigWigAverageOverBed "+data[i]['fn']+" "+path+"mm10_200bin_tab.bed "+path+"signal/"+signal_file_name+"_200bin.tab"+"\n")
        code.append("""awk '{print $5}' """+path+"signal/"+signal_file_name+"_200bin.tab > "+path+"signal/"+signal_file_name+"_200bin.txt"+"\n")
        input.append(biosample+" "+data[i]['label']+" "+ path + "signal/" + signal_file_name + "_200bin.txt"+"\n")


write_file(code, path+"Jul24.get_signal_code_for_all_chr.sh")
write_file(input, path+"IDEAS_8hm_atac_dname.input")
