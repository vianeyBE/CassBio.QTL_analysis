import sys
import getopt
import vcf
import os
import time
import fileinput
import re 

###############################################################################
##### GET ARGUMENTS FROM THE COMMAND LINE ####
# To provide the vcf file name, parental genotypes, 
# and the output file name in the command line
# -i path/to/inputFile (--input=)
# -m mother name (--mother-name=)
# -f father name (--father-name=)
# -o outputFile (--output=)
# --three-allele (optional, to consider locus with three alleles)
# call the script like this:
# python encoding-snp.py -i path/to/input -m mother_name -f father_name -o outputFile (--three-allele)
###############################################################################
arguments = sys.argv[1:]
print arguments
if len(arguments) < 8:
	print "incorrect number of arguments"
	sys.exit(1)

optList, arg = getopt.getopt(arguments,
	'i:m:f:o:',
	["input-file=", "mother-name=", "father-name=", "output=", "three-allele", "four-allele"])

#optList = [('-i', 'D:/OneDrive - CGIAR/Documents/Vianey/04 Whitefly - Genetic Map/snp-encoding/cassava.filtered.vcf'), ('-m', 'ECU72_FEMALE'), ('-f', 'COL2246_MALE'), ('-o', 'cassava.filtered'), ('--three-allele', '')]
#### SAVE ARGUMENTS AS VARIABLES #### 
three_allele = "False"
four_allele = "False"
for name, value in optList:
	if name in ('-i', '--input'):
		inputFile = value
		print "VCF file:", value
	if name in ('-m', '--mother-name'):
		motherName = value
		print "Mother genotype:", value
	if name in ('-f', '--father-name'):
		fatherName = value
		print "Father genotype:", value
	if name in ('-o', '--output'):
		outputFile = value
		print "Output file:", value	
	if name in ("--three-allele"): 
		three_allele = "True"
		print "Locus with three alleles (segregation <fexeg>) will be consider"
	if name in ("--four-allele"):
		four_allele = "True"
		print "Locus with four alleles (segregation <abxcd>) will be consider"        

###############################################################################
#### CREATE OUPUT FILE ####
outputFile1 = ''.join([outputFile, '.loc'])
date = time.strftime("%d/%m/%Y")
f = open(outputFile1, "w")
f.write("""; locus genotype file for JoinMap5
; created on %s

""" % date)
f.close()

outputFile2 = ''.join([outputFile, '.positions.txt'])
f2 = open(outputFile2, "w")
f2.write("""Chromosome	Positions	JoinMapName\n""")
f2.close()

#### READ VCF FILE AND READ DATA ####
vcfFile = vcf.Reader(open(inputFile, 'r'))
name = os.path.splitext(os.path.basename(inputFile))[0]
nind = len(vcfFile.samples)-2
print "Number of individuals:", nind
nloc = "pendiente"
popt = "CP" 

### WRITE HEADER ###
header_template = """name = {name}
popt = {popt}
nloc = {nloc}
nind = {nind}

"""
header_data = {
 "name":name,
 "popt":popt,
 "nloc":nloc,
 "nind":nind
 }

f = open(outputFile1, "a")
f.write(header_template.format(**header_data))
f.close()

print "Proccessing... It will take some time, be patient"

#### ENCODING AND WRITE GENOTYPE DATA ####
genodata_template = """{locus}  <{segre}> 
  {geno_code} 
"""
postitions_template = "{chromosome}	{position}	{joinmapName}\n"
nloc = 0
###############################################################################
#### FUNCTION TO WRITE GENOTYPE DATA ####
def savegeno(record, segre, geno_list):
	chromosome = record.CHROM.replace("Scaffold", "S").replace("Chromosome", "C")
	position = str(record.POS)
	
	locus = chromosome + "_" + position + "_"  + ''.join(str(alt)).replace('[', '').replace(']', '').replace(',','')
	position_data = {
			"chromosome": chromosome,
			"position": position,
			"joinmapName": locus}
	genodata_data = {
			"locus":locus,
			"segre":segre, 
			"geno_code":' '.join(geno_list)
			}
	f = open(outputFile1, "a")
	f.write(genodata_template.format(**genodata_data))
	f.close()
	f2 = open(outputFile2, "a")
	f2.write(postitions_template.format(**position_data)) 
	f2.close()


#### FUNCTION TO ENCODE THREE ALLELES ####
def threeallele(mgeno_allele, fgeno_allele): 
	segre = "efxeg"
	for mother_allele in mgeno_allele:
		for father_allele in fgeno_allele:
			if mother_allele == father_allele:
				e_allele = mother_allele
	mgeno_allele.remove(e_allele)
	f_allele = mgeno_allele[0]			
	fgeno_allele.remove(e_allele)
	g_allele = fgeno_allele[0]			

	for sample in vcfFile.samples:		
		if sample not in parents:
			childgeno = record.genotype(sample)["GT"]
			if childgeno == '/'.join([e_allele, e_allele]):
				geno_list.append("ee")
			elif childgeno == '/'.join([e_allele, f_allele]):
				geno_list.append("ef")
			elif childgeno == '/'.join([e_allele, g_allele]):
				geno_list.append("eg")
			elif childgeno == '/'.join([f_allele, g_allele]):
				geno_list.append("fg")
			else:
				geno_list.append("--")
	savegeno(record, segre, geno_list) 
	
###############################################################################
for record in vcfFile:
	
	#### READ PARENTS' DATA AND LOCUS
	parents = [motherName, fatherName]
	mgeno = record.genotype(motherName)["GT"]
	fgeno = record.genotype(fatherName)["GT"]
	alt = record.ALT
	ref = record.REF
	geno_list = []
	 
	#### CHECK IF REF IS SNP, PARENTS'GENOTYPES ARE KNOWN, AND THERE ARE AT LEAST TWO ALLELES
	check1 = "False"
	if len(ref) == 1 and mgeno not in ("./.") and fgeno not in ("./.") and len(set(re.split("/|\|", mgeno) + re.split("/|\|",fgeno)))>1:
		check1 = "True"
	
	#### CHECK IF ALT IS SNP
	check2 = "True"
	for alternative in alt: 
		if len(alternative) != 1 or alternative == '*': 
			check2 = "False"

	#### FOUR ALTERNATIVE ALLELES
	if len(alt) > 2 and four_allele == "True" and check1 == "True" and check2 == "True":
		#print "four alleles"
		a = re.split("/|\|", mgeno)[0]
		b = re.split("/|\|", mgeno)[1]
		c = re.split("/|\|", fgeno)[0]
		d = re.split("/|\|", fgeno)[1]
		nloc = nloc + 1
		segre = "abxcd"
		for sample in vcfFile.samples:
			if sample not in parents:
				childgeno = record.genotype(sample)["GT"]
				if childgeno in ( a+'/'+c , a+"|"+c ):
					geno_list.append("ac")
				elif childgeno in ( a+'/'+d , a+"|"+d ):
					geno_list.append("ad")
				elif childgeno in ( b+'/'+c , b+"|"+c ):
					geno_list.append("bc")
				elif childgeno in ( b+'/'+d , b+"|"+d ):
					geno_list.append("bd")
				else:
					geno_list.append("--")
		savegeno(record, segre, geno_list)
		
	#### THREE ALTERNATIVE ALLELES
	if (len(alt) == 2 and check1 == "True" and check2 == "True" and three_allele == "True"):
		mgeno_allele = re.split("/|\|", mgeno)
		fgeno_allele = re.split("/|\|", fgeno)	  
		parents_alleles = mgeno_allele + fgeno_allele
		
		#### THREE AND HETEROZYGOUS ALLELES 
		if len(set(parents_alleles))==3 and len(set(mgeno_allele))==2 and len(set(fgeno_allele))==2:
			nloc = nloc + 1
			threeallele(mgeno_allele, fgeno_allele)
			
		#### THREE ALLELES BUT SEGREGATION AS TWO ALLELES: nnxnp
		if len(set(parents_alleles))==3 and len(set(mgeno_allele))==1 and len(set(fgeno_allele))==2:
			a = mgeno_allele[0]
			b = fgeno_allele[0]
			c = fgeno_allele[1]
			#print "3 as 2", mgeno_allele, fgeno_allele, a, b, c			
			nloc = nloc + 1
			segre = "nnxnp"
			for sample in vcfFile.samples:
				if sample not in parents:
					childgeno = record.genotype(sample)["GT"]
					if childgeno in ( a+'/'+a , a+"|"+a ):
						geno_list.append("nn")
					elif childgeno in ( b+'/'+c , b+"|"+c ):
						geno_list.append("np")
					else:
						geno_list.append("--")
			savegeno(record, segre, geno_list)

		#### THREE ALLELES BUT SEGREGATION AS TWO ALLELES: lmxll
		if len(set(parents_alleles))==3 and len(set(mgeno_allele))==2 and len(set(fgeno_allele))==1:
			a = mgeno_allele[0]
			b = mgeno_allele[1]
			c = fgeno_allele[0]
			print "3 as 2", mgeno_allele, fgeno_allele, a, b, c			
			nloc = nloc + 1
			segre = "lmxll"
			for sample in vcfFile.samples:
				if sample not in parents:
					childgeno = record.genotype(sample)["GT"]
					if childgeno in ( a+'/'+b , a+"|"+b ):
						geno_list.append("lm")
					elif childgeno in ( c+'/'+c , c+"|"+c ):
						geno_list.append("ll")
					else:
						geno_list.append("--")
			savegeno(record, segre, geno_list)
			   			   
	#### TWO ALLELES
	if (len(alt) == 1 and check1 == 'True' and check2 == 'True'):
		#twoallele(mgeno, fgeno, nloc)
		geno_list = []
		if mgeno in ("0/1", "0|1") and fgeno in ("0/1", "0|1"):
			nloc = nloc + 1
			segre = "hkxhk"
			for sample in vcfFile.samples:		
				if sample not in parents:
					childgeno = record.genotype(sample)["GT"]
					if childgeno in ("0/0", "0|0"):
						geno_list.append("hh")
					elif childgeno in ("0/1", "0|1"):
						geno_list.append("hk")
					elif childgeno in ("1/1", "1|1"):
						geno_list.append("kk")
					else:
						geno_list.append("--")
			savegeno(record, segre, geno_list)
						
		elif mgeno in ("0/1", "0|1") and fgeno in ("0/0", "0|0","1/1", "1|1"):
			nloc = nloc + 1
			segre = "lmxll"
			for sample in vcfFile.samples:
				if sample not in parents:
					childgeno = record.genotype(sample)["GT"]
					if childgeno in ("0/0", "0|0","1/1", "1|1"):
						geno_list.append("ll")
					elif childgeno in ("0/1", "0|1"):
						geno_list.append("lm")
					else:
						geno_list.append("--")
			savegeno(record, segre, geno_list)
										 
		elif mgeno in ("0/0", "0|0","1/1", "1|1") and fgeno in ("0/1", "0|1"):
			nloc = nloc + 1
			segre = "nnxnp"
			for sample in vcfFile.samples:
				if sample not in parents:
					childgeno = record.genotype(sample)["GT"]
					if childgeno in ("0/0", "0|0","1/1", "1|1"):
						geno_list.append("nn")
					elif childgeno in ("0/1", "0|1"):
						geno_list.append("np")
					else:
						geno_list.append("--")
			savegeno(record, segre, geno_list)
		else:
			#print "check this"
			#print "position", record.CHROM, "-" , record.POS		  
			#print "REF", ref
			#print "ALT", alt
			#print "mother", mgeno
			#print "father", fgeno
			pass 
						
print "Number of SNPs processed:", nloc

f = open(outputFile1, "a")
f.write("""
individual names:
	
""" )

for sample in vcfFile.samples:
	if sample not in parents:
		f.write(sample + '\n')
f.close()

for line in fileinput.input(outputFile1, inplace = 1):
	print line.replace("nloc = pendiente", ("nloc = " + str(nloc))).strip()

print "DONE!"