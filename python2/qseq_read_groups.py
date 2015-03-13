#!/usr/bin/env python

import os
from optparse import OptionParser, OptionGroup
import re

#path to metadata
#/scratch1/tmp/myourshaw/gmd/new/qseq/gmd.metadata.txt

description = """
foo bar
"""

parser = OptionParser(description=description)
(options, args) = parser.parse_args()
metadata = args[0]

line_count = 0

for line in open (metadata):
	line_count +=1
	line = line.strip()
	#header
	if line.startswith("#"):
		line = line.lstrip("#")
		col_names = line.split("\t")
		DIR = col_names.index("dir")
		LANE = col_names.index("lane")
		BARCODE = col_names.index("barcode")
		RG_ID = col_names.index("read_group_id")
	else:
		fields = line.split("\t")
		dir = fields[DIR]
		file_re = re.compile(r"s_"+fields[LANE]+r"_(?P<read>\d)_(?P<tile>\d{4})_qseq.txt" )
		files = os.listdir(dir)
		for file in files:
			base = os.path.basename(file)
			m = file_re.match(base)
			if not m: continue
			new_file = fields[RG_ID]+'.'+base
			os.rename(os.path.join(dir,file), os.path.join(dir,new_file))
			print os.path.join(dir,file), os.path.join(dir,new_file)
print "done"

