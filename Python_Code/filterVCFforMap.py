# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 07:29:45 2021
Script to replace FilterVCFforMap.groovy. This new script will perform the same task as the previous groovy version, 
@author: VBarrera
"""

#import sys
#import getopt
import vcf 
import os
from collections import namedtuple
from scipy.stats import chisquare
#from model import _Call, _Record, make_calldata_tuple
#import time
#import fileinput
#import re #regular expression

#pkg_resources.get_distribution('pyvcf').version
directory = "D:/OneDrive - CGIAR/Documents/Vianey/05 Betacarotene/filterVCFforMap/"
inputFile = "CM8996_original.vcf"
#list_indiv = "F1.txt"
motherName = "ECU72_FEMALE"
fatherName = "COL2246_MALE"
p_threshold = 0.05

vcfFile = vcf.Reader(open(inputFile, 'r'))
name = os.path.splitext(os.path.basename(inputFile))[0]
output = name+".initial2.vcf"
nind = len(vcfFile.samples)

#f1 = open(list_indiv, 'r')
#lines = f1.readlines()
#f1_family = []
#for line in lines:
#	f1_family.append(line)
#	#print(line)
#f1.close()

# Open new file to write VCF output
vcf_writer = vcf.Writer(open(directory+output, 'w'), vcfFile)


for record in vcfFile:

	#### read information
	parents = [motherName, fatherName]
	mgeno = record.genotype(motherName)["GT"]
	fgeno = record.genotype(fatherName)["GT"]
	alt = record.ALT
	ref = record.REF
	geno_list = []
	skip=False

	#### Check if there are missing genotypes in the parents
	if mgeno in ("./.") or fgeno in ("./."):
		#print(mgeno, fgeno)
		skip = True

	#### Check if both parents are homozygous
	if (mgeno in ("0/0") and fgeno in ("0/0")) or\
		(mgeno in ("1/1") and fgeno in ("1/1")) or\
		(mgeno in ("1/1") and fgeno in ("0/0")) or\
		(mgeno in ("0/0") and fgeno in ("1/1")):
		skip = True
		#print(mgeno, fgeno)

	calls = 0
	no_calls = 0
	AA = 0
	Aa = 0
	aa = 0
	chi_val = 0
	p_val = 0

	for sample in vcfFile.samples:

		if sample not in parents:
			#print(sample)
			childgeno = record.genotype(sample)["GT"]

			#### Check progeny genotype inconsistencies
			if childgeno in ("0/0") and ((mgeno in ("1/1") and fgeno in ("0/1")) or (mgeno in ("0/1") and fgeno in ("1/1"))):
				#record.genotype(sample)["GT"] = "./."
				#new_gt = "./."
				#print("original", record.genotype(sample).data)
				new_CallData = namedtuple('CallData', record.genotype(sample).data._fields)
				calldata = ["./."] + list(record.genotype(sample).data[1:])
				record.genotype(sample).data = new_CallData(*calldata)
				#print("update", record.genotype(sample).data)
			elif childgeno in ("1/1") and ((mgeno in ("0/0") and fgeno in ("0/1")) or (mgeno in ("0/1") and fgeno in ("0/0"))):
				#record.genotype(sample)["GT"] = "./."
				#print(sample, mgeno, fgeno, childgeno)
				#new_gt = "./."
				new_CallData = namedtuple('CallData', record.genotype(sample).data._fields)
				calldata = ["./."] + list(record.genotype(sample).data[1:])
				record.genotype(sample).data = new_CallData(*calldata)
			else:
				pass
			# end if child genotype

			# Segregation Distortion
			if childgeno in ("0/1"): Aa+=1
			elif (childgeno in ("0/0") and mgeno in ("0/0") and fgeno in ("0/1")): AA+=1
			elif (childgeno in ("1/1") and mgeno in ("1/1") and fgeno in ("0/1")): AA+=1
			elif (childgeno in ("0/0") and mgeno in ("0/1") and fgeno in ("0/0")): aa+=1
			elif (childgeno in ("1/1") and mgeno in ("0/1") and fgeno in ("1/1")): aa+=1
			elif (childgeno in ("0/0") and mgeno in ("0/1") and fgeno in ("0/1")): AA+=1
			elif (childgeno in ("1/1") and mgeno in ("0/1") and fgeno in ("0/1")): aa+=1
			else: no_calls+=1
			# end if segregation distorion
		# end if sample 
	# end for samples

	calls = AA + Aa + aa
	if calls>0:
		## hkxhk AA-Aa-aa
		if (mgeno in ("0/1") and fgeno in ("0/1")):
			chi_val = chisquare(f_obs=[AA, Aa, aa], f_exp=[0.25*calls, 0.50*calls, 0.25*calls])[0]
			p_val = chisquare(f_obs=[AA, Aa, aa], f_exp=[0.25*calls, 0.50*calls, 0.25*calls])[1]
		## nnxnp male backcross AA-Aa-0
		elif (AA!=0 and Aa!=0 and aa==0):
			chi_val = chisquare(f_obs=[AA, Aa], f_exp=[0.5*calls, 0.5*calls])[0]
			p_val = chisquare(f_obs=[AA, Aa], f_exp=[0.50*calls, 0.50*calls])[1]
		## lmxll female backcross 0-Aa-aa
		elif (AA==0 and Aa!=0 and aa!=0):
			chi_val = chisquare(f_obs=[Aa, aa], f_exp=[0.5*calls, 0.5*calls])[0]
			p_val = chisquare(f_obs=[Aa, aa], f_exp=[0.50*calls, 0.50*calls])[1]
		## O-Aa-0, when parents are 0/0 x 1/1 - It is not informative in a F1 population 
		elif (AA==0 and Aa!=0 and aa==0):
			p_val = 1.0
			skip = True
		## AA-0-aa, dominant genotypes h- or k-
		elif (AA!=0 and Aa==0 and aa!=0):
			p_val = 0.0
			skip = True
		# end if chi square test
	else:
		skip = True
	# end if calls - chi square test
		
	## Write VCF - save variants with skip false
	if (skip==False and p_val>=p_threshold):
		vcf_writer.write_record(record)
	## end if filter skip
# end for records-variants
vcf_writer.close()
vcfFile.close()
print("Done")
