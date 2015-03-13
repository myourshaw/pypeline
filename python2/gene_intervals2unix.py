#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
import MySQLdb
import codecs

description = """
This script creates bed and interval_list files from ucsc gene tables.
Each row represents an exon of a gene transcript. exonStart and exonEnd are
expanded by 2 bases to include essential splice site locations.
"""
parser = OptionParser()
(options, args) = parser.parse_args()

#/Volumes/raid/resources/intervals/b36.ccdsGene.exons.bed.dos
for file in args:
	if not file.endswith(".bed.dos"):
		continue
	bed = open(file.rsplit('.',2)[0]+'.bed','w',)
	interval_list = open(file.rsplit('.',2)[0]+'.interval_list','w')
	for line in codecs.open(file,encoding='utf_8_sig'):
		line = line.strip()
		if line.startswith('#') or len(line) == 0:
			continue
		fields = line.strip().split("\t")
		bedline = "\t".join(fields)+"\n"
		bed.write(bedline)
		intervalline = "\t".join((fields[0],str(int(fields[1])+1),fields[2],'+',fields[3]))+"\n"
		interval_list.write(intervalline)
