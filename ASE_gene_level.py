import sys,re,fileinput,os,math

Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 2:
        print "Usage: ChiSquare_file_annotated_FDRflagged RefseqgeneID Output_file"  
        sys.exit()

File_Chi = Argument[0]
ref2gene = Argument[1]
output = open(Argument[2],"w")

Ref2gene = {}

for line in fileinput.input([ref2gene]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')

        if rowlist[3].split(".")[0] not in Ref2gene:
        	Ref2gene[rowlist[3].split(".")[0]] = rowlist[1]+"\t"+rowlist[-1]

SNP_done = []
Gene_ASE = {}
Gene_order = []
Gene_Total = {}
Gene2ref = {}
Gene2ref_many = {}
Allele_ratio_B6 = {}
Allele_ratio_D2 = {}
SNP_ratio = {}

for line in fileinput.input([File_Chi]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')

        if line.startswith("#") or rowlist[0] == "chrX" or rowlist[-3] != "PASS":
		continue
        else:
		if rowlist[0]+"\t"+rowlist[1] in SNP_done:
                        continue

		if rowlist[11] not in Gene_order:
			Gene_order.append(rowlist[11])
		
		if rowlist[11] not in Gene2ref_many:
			Gene2ref_many[rowlist[11]] = {}
			Gene2ref_many[rowlist[11]][rowlist[12]] = 1
		else:
			if rowlist[12] in Gene2ref_many[rowlist[11]]:
				Gene2ref_many[rowlist[11]][rowlist[12]] = Gene2ref_many[rowlist[11]][rowlist[12]] + 1
			else:
				Gene2ref_many[rowlist[11]][rowlist[12]] = 1
		
		#if rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[11] in SNP_done:
		#	continue

                if rowlist[11] not in Gene_Total:
			Gene2ref[rowlist[11]] = rowlist[12]
			Gene_ASE[rowlist[11]] = 0
			Gene_Total[rowlist[11]] = 0
			Allele_ratio_B6[rowlist[11]] = 0
			Allele_ratio_D2[rowlist[11]] = 0
			SNP_ratio[rowlist[11]] = []	
	
			if rowlist[-1] == "PASS":
				SNP_ratio[rowlist[11]].append(rowlist[9])
				Gene_ASE[rowlist[11]] = 1
			
				if float(rowlist[9]) >= 0.55:
					Allele_ratio_B6[rowlist[11]] = 1 
				
				if float(rowlist[9]) <= 0.45:
					Allele_ratio_D2[rowlist[11]] = 1

			Gene_Total[rowlist[11]] = 1

		else:

			if rowlist[-1] == "PASS":
				SNP_ratio[rowlist[11]].append(rowlist[9])
                        	Gene_ASE[rowlist[11]] = Gene_ASE[rowlist[11]] + 1

				if float(rowlist[9]) >= 0.55:
                                        Allele_ratio_B6[rowlist[11]] = Allele_ratio_B6[rowlist[11]] + 1

				if float(rowlist[9]) <= 0.45:
                                        Allele_ratio_D2[rowlist[11]] = Allele_ratio_D2[rowlist[11]] + 1
				
                        Gene_Total[rowlist[11]] =  Gene_Total[rowlist[11]] + 1
	
		#SNP_done.append(rowlist[0]+"\t"+rowlist[1]+"\t"+rowlist[11])
		SNP_done.append(rowlist[0]+"\t"+rowlist[1])

#print Gene_order
output.write("GENE\tGENEID\tGENE\tREFSEQID\tTESTABLE_MARKERS\tASE_MARKERS\tASE_MARKERS_B6_0.55\tASE_MARKERS_D2_0.45\tASE_MARKERS_VALUES\tREFSEQ_IDS_markers\n")

for gene in Gene_order:
	print gene
	print Gene2ref[gene]
	print Ref2gene[Gene2ref[gene]]
	output.write(str(gene)+"\t"+str(Ref2gene[Gene2ref[gene]])+"\t"+str(Gene2ref[gene])+"\t"+str(Gene_Total[gene])+"\t"+str(Gene_ASE[gene])+"\t"+str(Allele_ratio_B6[gene])+"\t"+str(Allele_ratio_D2[gene])+"\t"+str(";".join(SNP_ratio[gene]))+"\t")
	for refseq in Gene2ref_many[gene]:
		output.write(str(refseq)+"["+str(Gene2ref_many[gene][refseq])+"]"+";")
	output.write("\n")
output.close()

