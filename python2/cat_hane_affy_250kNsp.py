#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
import glob

parser = OptionParser()
(options, args) = parser.parse_args()

def get_files(names):
	return (file for file in names if os.path.isfile(file))

in_glob = "/Users/myourshaw/lab/MMS/genotyping_redo/hane/*.TXT"
out = open("/Users/myourshaw/lab/MMS/genotyping_redo/work/hane_chr9_genotypes2.txt","w")
sep = "\t"
header = """sample	SNP_ID	Chromosome	Physical_Position	dbSNP_RS_ID	Call	Confidence"""
header_fields = header.split(sep)


out.write(header+'\n')
for file in get_files(glob.iglob(in_glob)):
	print file
	line_number = 0
	for line in open(file):
		line_number += 1
		line_fields=line.rstrip('\n').split(sep)
		if line.startswith('Dynamic'): continue
		if line_fields[1] == 'SNP ID':
			sample_indices = []
			samples = []
			confidence_indices = []
			for i in range(5,len(line_fields)):
				if  line_fields[i].endswith('_Call'):
					sample_indices.append(i)
					samples.append(re.sub(r'_Call','',line_fields[i]))
					if len(line_fields) > i+1 and re.sub(r'_Call','',line_fields[i]) == re.sub(r'_Confidence','',line_fields[i+1]):
						confidence_indices.append(i+1)
		if line_fields[2] == '9':
			for i in range(len(sample_indices)):
				ix = sample_indices[i]
				sample = re.sub(r'_Call','',samples[i])
				if len(confidence_indices) > i:
					confidence = line_fields[confidence_indices[i]]
				else:
					confidence = ''
				out.write("%s\t%s\t%s\t%s\n" % (sample,sep.join(line_fields[1:5]),line_fields[ix],confidence))
