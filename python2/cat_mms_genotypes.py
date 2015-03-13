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

in_glob = "/Users/myourshaw/lab/MMS/genotyping_redo/work/*.txt"
out = open("/Users/myourshaw/lab/MMS/genotyping_redo/work/all_genotypes","w")
sep = "\t"
header = """sample	SNP_ID	Chromosome	Physical_Position	dbSNP_RS_ID	Call	Confidence"""
header_fields = header.split(sep)


out.write(header+'\n')
for file in get_files(glob.iglob(in_glob)):
	line_number = 0
	for line in open(file):
		line_number += 1
		line_fields=line.rstrip('\n').split(sep)
		if line.startswith('Dynamic'): continue
		if line_fields[1] == 'SNP ID':
			sample = re.sub(r'_Call','',line_fields[5])
			continue
		out.write("%s\t%s\n" % (sample,sep.join(line_fields[1:])))
