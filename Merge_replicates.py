import sys,re,fileinput,os
from statlib import stats
reload(stats)

#print stats.chisqprob(250.0,1)

Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 2:
        print "Usage: Path_to_folder_containing_alleleic_ratio_files New_file" 
        sys.exit()

dpath = Argument[0]
output = open(str(Argument[1]),"w")

Listoffile = []

def dir(patharray):
	Listdir = []
    	for infile in patharray:
        	Listdir.append(os.path.join(dpath,infile))
    	return Listdir    

Listoffile = dir(os.listdir(dpath))

Allelic_ratio = {}

for Filepath in Listoffile:
	for (path, dirs, files) in os.walk(Filepath):
        	for file in files:
                	if file.endswith("_allelic_ratio.txt") and not os.path.getsize(path+"/"+file) == 0:
                        	if path+"/"+file not in Allelic_ratio:
                                	Allelic_ratio[path+"/"+file] = (path+"/"+file)

print Allelic_ratio

def numeric_compare(x, y):
        x1 = int(x)
        y1 = int(y)
        return x1 - y1

import math

def average(x):
	assert len(x) > 0
    	return float(sum(x)) / len(x)

Q_crit = {'3':['0.941','0.970','0.994'],
'4':['0.765','0.829','0.926'],
'5':['0.642','0.710','0.821'],
'6':['0.560','0.625','0.740']}

def median(mylist):
        sorts = sorted(mylist,key=float)
        length = len(sorts)

	if length == 1:
		return sorts[0]

	if length == 2:
		return float(sorts[0]+sorts[1])/2

        if not length % 2:
                median = (sorts[length/2] + sorts[length/2-1])/2.0
        else: 
                median =  sorts[length/2]

	return median

def interquartile(x):
	
	count_dict = x
        ratio_dict = {}
        for file in count_dict:
                ratio_dict[file] = float(count_dict[file][0]+count_dict[file][3])/float(count_dict[file][0]+count_dict[file][3]+count_dict[file][1]+count_dict[file][2])

        sorts = sorted(ratio_dict.values(),key=float)
        Outlier = []

	lower_quartile = 0
	upper_quartile = 0

    	length = len(sorts)

    	if not length % 2:
        	lower_quartile = median(sorts[0:length/2])
		upper_quartile = median(sorts[length/2:])
    	else:
		lower_quartile = median(sorts[0:(length-1)/2])
		upper_quartile = median(sorts[(length+1)/2:]) 

	IQR = 0
	IQR = upper_quartile-lower_quartile
	
	lower_limit = lower_quartile-(1.5*IQR)
	upper_limit = upper_quartile+(1.5*IQR)

	for a in sorts:
		if a < lower_limit or a > upper_limit:
			Outlier.append(a)

	#print sorts
	#print Outlier 

	return_list = []
	for a in ratio_dict:
                if ratio_dict[a] not in Outlier:
                        return_list.append(a)

        return return_list

def dixon_outlier(x):
	count_dict = x
	ratio_dict = {}
	for file in count_dict:
		ratio_dict[file] = float(count_dict[file][0]+count_dict[file][3])/float(count_dict[file][0]+count_dict[file][3]+count_dict[file][1]+count_dict[file][2])

	ratio_sort = sorted(ratio_dict.values(),key=float)
	Outlier = []
	
	print ratio_sort

	for i in range(0,len(ratio_sort)-1):
		Q_value = 0.0
		Q_value = (ratio_sort[i+1]-ratio_sort[i])/(ratio_sort[-1]-ratio_sort[0])
		print Q_value,Q_crit[str(len(ratio_sort))][0]
				
		if Q_value > Q_crit[str(len(ratio_sort))][0]:
			Outlier.append(ratio_sort[i])			 

	if len(Outlier) >= 1:
		print Outlier

	return_list = []
	for a in ratio_dict:
		if ratio_dict[a] not in Outlier:
			return_list.append(a)
	
	return return_list

Replicates = {}
SNPs = {}
Unsort_SNPs = {}
Sort_SNPs = {}
Map_SNPs = {}

for filepath in Allelic_ratio:
	for line in fileinput.input([filepath]):
        	rowlist = []
        	rowlist = (line.rstrip("\n")).split('\t')

        	if line.startswith("#") or rowlist[0] == "chrX":
                	continue
        	else:
			if int(rowlist[4])+int(rowlist[5])+int(rowlist[6])+int(rowlist[7]) < 10:
				continue

			if rowlist[0] not in Unsort_SNPs:
                                Unsort_SNPs[rowlist[0]] = []
                                Unsort_SNPs[rowlist[0]].append(rowlist[1])

                                Map_SNPs[rowlist[0]] = {}
                                Map_SNPs[rowlist[0]][rowlist[1]] = rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]
                        else:
				if rowlist[1] not in Unsort_SNPs[rowlist[0]]:
                                	Unsort_SNPs[rowlist[0]].append(rowlist[1])
                                	Map_SNPs[rowlist[0]][rowlist[1]] = rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]

				
                	if rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]  not in SNPs:
				SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]] = {}
				SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]][filepath] = []
				SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]][filepath].append(int(rowlist[4]))
				SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]][filepath].append(int(rowlist[5]))
				SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]][filepath].append(int(rowlist[6]))
				SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]][filepath].append(int(rowlist[7]))
			else:
				SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]][filepath] = []
                                SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]][filepath].append(int(rowlist[4]))
                                SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]][filepath].append(int(rowlist[5]))
                                SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]][filepath].append(int(rowlist[6]))
                                SNPs[rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[2]+"\t"+rowlist[3]][filepath].append(int(rowlist[7]))
	
Chisqr_saved = {}

for chrom in Unsort_SNPs:
	Sort_SNPs[chrom] = []
	Sort_SNPs[chrom] = sorted(Unsort_SNPs[chrom],key=int)	

chromosome = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"]

for chrome in chromosome:
	for position in Sort_SNPs[chrome]:
		snp = ""
		snp = Map_SNPs[chrome][position]

		if len(SNPs[snp]) >= 3:
			passed = []
			#print snp
			passed = interquartile(SNPs[snp])

			if len(passed) >= 3: 
				countA = 0
				countB = 0

				for file in SNPs[snp]:
					count = []
					count = SNPs[snp][file]
					if file in passed:
						countA = countA + count[0]+count[3]
						countB = countB + count[1]+count[2]

				sum = countA+countB
	
				o1 = float(countA)
				o2 = float(countB)
                		e1 = float(sum)/2.0
                		e2 = float(sum)/2.0
		
                		chisqr_value = 0.0
                		chisqr_value  = (((o1-e1)*(o1-e1))/e1)+(((o2-e2)*(o2-e2))/e2)

    
                		if chisqr_value not in Chisqr_saved:
                        		p_value = 1.00
                        		p_value = stats.chisqprob(chisqr_value,1)
                        		Chisqr_saved[chisqr_value] = p_value
				else:
					p_value = Chisqr_saved[chisqr_value]
	
				output.write(str(snp)+"\t"+str(countA)+"\t"+str(countB)+"\t"+str(sum)+"\t"+str(chisqr_value)+"\t"+str(p_value)+"\t"+str(float(o1/sum))+"\t"+str(len(passed))+"\n")

output.close()
