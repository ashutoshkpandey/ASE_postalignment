"""
This function implements False Discovery rate (FDR) procedure developed by Yoav Benjamini and Yosef Hochberg.
Arguments: List containing the P-values , significance level(Alpha) 
Output: List containing P-values qualifying as significant P-values after FDR correction
Reference: Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing by Yoav Benjamini and Yosef Hochberg.
"""

def FDRAnalysis(Pvaluelist,signlevel):
	Pvaluelistsort = Pvaluelist[:]
	Pvaluelistsort.sort(lambda a,b: cmp(float(a), float(b))) 
	for i in range(len(Pvaluelistsort)-1,-1,-1):
		if Pvaluelistsort[i] <= ((i*signlevel)/(len(Pvaluelistsort))):
			return float(Pvaluelistsort[i])
			
import sys,re,fileinput,os,math

Argument = []
Argument = sys.argv[1:]

if (len(Argument)) < 3:
        print "Usage: Chisquare_indel_flagged alpha-pvalue Output"
        sys.exit()

Pvalues = []
adjpvalue = 0.0 

for line in fileinput.input([Argument[0]]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')
	if rowlist[-2] == "PASS":
		Pvalues.append(float(rowlist[8]))

adjpvalue = FDRAnalysis(Pvalues,float(Argument[1]))
print adjpvalue 

output = open(Argument[2],'w')

for line in fileinput.input([Argument[0]]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')
	line = line.rstrip("\n")
	if rowlist[-2] != "PASS":
		output.write(str(line)+"\tMARKER_FAIL\n")
		continue
        if float(rowlist[8]) <= adjpvalue:
		output.write(str(line)+"\tPASS\n")
	else:
		output.write(str(line)+"\tPVALUE_FAIL\n")		

output.close()
