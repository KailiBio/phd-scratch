import requests
import re
import sys

# functions
def getRefNumber(gene_id):
    url = "https://www.ncbi.nlm.nih.gov/gene/" + gene_id
    response = requests.get(url)
    #refN = re.findall('See all GeneRIFs \(([0-9]+)\)', response.text)
    refN = re.findall("See all \(([0-9]+)\) citations", response.text)
    refLink = "https://www.ncbi.nlm.nih.gov/pubmed?LinkName=gene_pubmed&from_uid=%s"%gene_id
    if len(refN) == 0:
        return [0, "NA"]
    else:
        return [int(refN[0]), refLink]

# run here #

# help
if len(sys.argv) < 2:
    print "python " + sys.argv[0] + "in.gene.list(entrez_id in the first column)"
    sys.exit(0)
# search refence number for each gene and print to stdout
dict_refN = []
count = 0
for l in open(sys.argv[1]):
    dict_refN.append(l.strip())
for g in dict_refN:
    ret = getRefNumber(g.split("\t")[0])
    sys.stdout.write("%s\t%s\t%s\n"%(g, ret[0], ret[1]))
    count += 1
    sys.stderr.write("{0:<10}genes processed\r".format(count))