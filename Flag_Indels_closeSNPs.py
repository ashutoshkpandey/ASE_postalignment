import sys,re,fileinput,os


Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 3:
        print "Usage: ChiSquare_file_annotated VCF_file Output_file" 
        sys.exit()

File_Chi = Argument[0]
File_Indel = Argument[1]
output = open(Argument[2],"w")

Indel_info = {}

for line in fileinput.input([File_Indel]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')

        if line.startswith("#") or rowlist[0] == "chrX" or not re.search('INDEL;',line):
                continue

	if len(rowlist[3]) < 4 and len(rowlist[4]) < 4:
		continue

        else:
                if rowlist[0] not in Indel_info:
			Indel_info[rowlist[0]] = []
			Indel_info[rowlist[0]].append(rowlist[1])
		else:
			Indel_info[rowlist[0]].append(rowlist[1])

snp_location = 0

SNPs = {}

for line in fileinput.input([File_Chi]):
	rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')
	Status = "NotIndel"
	
        if line.startswith("#") or rowlist[0] == "chrX":
                continue
        else:
        	if rowlist[0] in Indel_info:
			for indel in Indel_info[rowlist[0]]:
			
				if int(indel) + 30 > int(rowlist[1]) and int(indel) < int(rowlist[1]):  
					output.write(str(line.rstrip("\n"))+"\tINDEL\t"+str(indel)+"\n")
					Status = "Indel"
					break

				if int(indel) - 30 < int(rowlist[1]) and int(indel) > int(rowlist[1]):  
                                        output.write(str(line.rstrip("\n"))+"\tINDEL\t"+str(indel)+"\n")
                                        Status = "Indel"
                                        break

				if rowlist[0]+"\t"+rowlist[1] in SNPs:
					output.write(str(line.rstrip("\n"))+"\tADJ_SNP\t"+str(snp_location)+"\n")
					Status = "Indel"
					break

				if int(rowlist[1])-snp_location < 60 and int(rowlist[1]) != snp_location:
					output.write(str(line.rstrip("\n"))+"\tADJ_SNP\t"+str(snp_location)+"\n")
					Status = "Indel"
					SNPs[rowlist[0]+"\t"+rowlist[1]] = "SNP"
                                        break
				

			if Status != "Indel":
				output.write(str(line.rstrip("\n"))+"\tPASS\tPASS"+"\n")

			snp_location = int(rowlist[1])
output.close()


