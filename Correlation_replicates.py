import sys,re,fileinput,os


Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 2:
        print "Usage: Allelic_ratio_file1 Allelic_ratio_file2" 
        sys.exit()

File1 = Argument[0]
File2 = Argument[1]
output = open(Argument[2],'w')

def numeric_compare(x, y):
        x1 = int(x)
        y1 = int(y)
        return x1 - y1

import math

def average(x):
	assert len(x) > 0
    	return float(sum(x)) / len(x)

def pearson(x, y):
	assert len(x) == len(y)
    	n = len(x)
    	assert n > 0
    	avg_x = average(x)
    	avg_y = average(y)
    	diffprod = 0
    	xdiff2 = 0
    	ydiff2 = 0
    	for idx in range(n):
        	xdiff = x[idx] - avg_x
        	ydiff = y[idx] - avg_y
        	diffprod += xdiff * ydiff
        	xdiff2 += xdiff * xdiff
        	ydiff2 += ydiff * ydiff

    	return diffprod / math.sqrt(xdiff2 * ydiff2)

SNP1 = {}

for line in fileinput.input([File1]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')

        if line.startswith("#") or rowlist[0] == "chrX":
                continue
        else:
                if rowlist[0]+"\t"+rowlist[1] not in SNP1 and int(rowlist[4])+int(rowlist[6]) >= 30:
			SNP1[rowlist[0]+"\t"+rowlist[1]] = 0.0
			SNP1[rowlist[0]+"\t"+rowlist[1]] = float(rowlist[4])/float(int(rowlist[4])+int(rowlist[6]))

SNP2 = {}

for line in fileinput.input([File2]):
	rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')

        if line.startswith("#") or rowlist[0] == "chrX":
                continue
        else:
        	if rowlist[0]+"\t"+rowlist[1] in SNP1 and int(rowlist[4])+int(rowlist[6]) >= 30:
                	SNP2[rowlist[0]+"\t"+rowlist[1]] = 0.0
                        SNP2[rowlist[0]+"\t"+rowlist[1]] = float(rowlist[4])/float(int(rowlist[4])+int(rowlist[6]))

Array1 = []
Array2 = []

from statlib import stats
reload(stats)

for snp in SNP2:
	output.write(str(SNP2[snp])+"\t"+str(SNP1[snp])+"\n")
	Array2.append(SNP2[snp])
	Array1.append(SNP1[snp])

	
print pearson(Array1,Array2)
print stats.pearsonr(Array1,Array2)
print len(Array2)

output.close()
