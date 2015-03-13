#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
from IndentedHelpFormatterWithNL import *

description = """
This script takes a bed file and converts it into the format expected by samtools (contig\tposition\n).
"""
parser = OptionParser( description=description, usage="usage: %prog [options] INPUT-FILE", formatter=IndentedHelpFormatterWithNL())
parser.add_option("-o", "--output-filename",help="Output file path [Default: %default]", default="stdout")
parser.add_option("-i", "--input-filename",help="Input file path [Default: %default]", default="stdin")
(options, args) = parser.parse_args()
bedfile=options.input_filename
targetfile = options.output_filename
try:
	out = open(targetfile,"w")
	for line in open(bedfile):
		fields = line.split('\t')
		for pos in range(int(fields[1])+1,int(fields[2])):
			target_line = "\t".join((fields[0],str(pos)))
			out.write (target_line+'\n')
except Exception as ex:
	print ex
