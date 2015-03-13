#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup

# Init cmd-line args
description = """
This script splits SOLiD .csfastq paired end files into a number of smaller files.
"""
parser = OptionParser(description=description)

#parser.add_option("-f", "--files", nargs=2, metavar="FILES", dest="files", help="Specifies two paired end .csfastq files to split")
parser.add_option("-r", "--reads", type="int", metavar="READS", dest="reads", help="Maximum number of reads per output file")
parser.add_option("-d", "--dir", metavar="DIR", dest="dir", help="Directory for output files")

(options, args) = parser.parse_args()

# -r 100 -d "/Volumes/raid/gmd/runs/temp" "/Volumes/raid/gmd/runs/python/GMD19_F3.csfasta" "/Volumes/raid/gmd/runs/python/GMD19_F5-P2.csfasta"
# -r 1000000 "/data/storage-1-02/abi/Solid_UCLA153/SOLiD_ucla153_20101102_PE_GMD/GMD19/results/primary.20101107000722537/reads/SOLiD_ucla153_20101102_PE_GMD_GMD19_F3.csfasta" "/data/storage-1-02/abi/Solid_UCLA153/SOLiD_ucla153_20101102_PE_GMD/GMD19/results/primary.20101114004221888/reads/SOLiD_ucla153_20101102_PE_GMD_GMD19_F5-P2.csfasta"
# -r 100 -d "/scratch1/tmp/myourshaw/gmd/runs/python" "/scratch1/tmp/myourshaw/gmd/runs/python/GMD19_F3.csfasta"

if(not os.path.isdir(options.dir)): os.makedirs(options.dir)
for file in args:
	inputCount = 0 #line count for input file
	outputReads = -1 #output read count
	file = os.path.realpath(file)
	(dir,base) = os.path.split(file)
	(name,ext)=os.path.splitext(base)
	outputFileCount = 0 #count of output files
	#get header
	for line in open(file):
		inputCount += 1
		if (line.startswith("#") or line.strip() == ''): continue
		elif line.startswith(">"):
			if (outputReads >= options.reads or outputReads == -1):
				if (outputReads != -1): out.close()
				outputFile = os.path.join(options.dir,name+"."+'%06d' % outputFileCount+ext)
				out = open(outputFile,"write")
				outputFileCount += 1
				outputReads = 0
				#out.write(header)
			out.write(line)
			outputReads += 1
		else: out.write(line)
