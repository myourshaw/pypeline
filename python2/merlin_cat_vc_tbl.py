#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
import glob
import math

#-o /Users/myourshaw/lab/ADHD/Paper/link/merlin_2011-03-15/vc/merlin_20110315.vc.cat.txt /Users/myourshaw/lab/ADHD/Paper/link/merlin_2011-03-15/vc/*-vc-chr*.tbl

usage = "usage: %prog [options] <list of merlin vc-chr*.tbl files>"
parser = OptionParser(usage=usage)
parser.add_option("-o", "--output",default="STDOUT", metavar="OUTPUTFILE", dest="out", help="write results to OUTPUTFILE [%default]")
(options, args) = parser.parse_args()
outputFile = options.out

def get_files(names):
	return (file for file in names if os.path.isfile(file))

def is_number(s):
	try:
		n = float(s)
		return False if (math.isinf(n) or math.isnan(n)) else True
	except ValueError:
		return False

header_in = "CHR\tPOS\tLABEL\tTRAIT\tH2\tLOD\tPVALUE\n"
header_out = "chr\tcM\tmarker\tvariable\th2\tlod\tpvalue"
if options.out != "STDOUT":
	saveout = sys.stdout                                     
	out = open(options.out, 'w')                             
	sys.stdout = out                                       
print header_out

vc_file_re = re.compile(r".*\/(?P<contig>.+)-vc-chr(?P<chrom>.+)\.tbl$",re.I)

for arg in args:
	for file in get_files(glob.iglob(arg)):
		m = vc_file_re.match(file)
		contig = m.group('contig')
		chrom = m.group('chrom')
		for line in open(file):
			if (line == header_in or line.strip() == ''): continue
			fields = line.rstrip('\n').split('\t')
			chr = fields[0]
			cm = fields[1]
			marker = fields[2]
			variable = fields[3]
			h2 = fields[4]
			lod = fields[5]
			pvalue = fields[6]
			if(is_number(cm) and is_number(h2) and is_number(lod) and is_number(pvalue)):
				if chr == '999':
					if chrom != '999': chr = chrom
					elif contig.lower() == 'autosome': chr = 'XY'
					elif contig.lower() == 'x': chr = 'X'
				print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chr,cm,marker,variable,h2,lod,pvalue))
			else:
				print (file,line)
