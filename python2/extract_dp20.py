#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup

parser = OptionParser()
(options, args) = parser.parse_args()

FORMAT = 8
DAD = 9
MOM = 10
KID = 11
for file in args:
	path_ext = os.path.splitext(os.path.abspath(file))
	dad_vcf = open(path_ext[0]+'.dad_dp20'+path_ext[1], 'w')
	mom_vcf = open(path_ext[0]+'.mom_dp20'+path_ext[1], 'w')
	kid_vcf = open(path_ext[0]+'.kid_dp20'+path_ext[1], 'w')
	for line in open(file):
		if line.startswith('#'):
			dad_vcf.write(line)
			mom_vcf.write(line)
			kid_vcf.write(line)
			continue
		line_strip = line.strip()
		if line_strip == '': contiune
		fields = line_strip.split('\t')
		DP = fields[FORMAT].split(':').index('DP')
		if fields[DAD] != './.':
			dad_dp = int(fields[DAD].split(':')[DP])
			if dad_dp <= 20: dad_vcf.write(line)
		if fields[MOM] != './.':
			mom_dp = int(fields[MOM].split(':')[DP])
			if mom_dp <= 20: mom_vcf.write(line)
		if fields[KID] != './.':
			kid_dp = int(fields[KID].split(':')[DP])
			if kid_dp <= 20: kid_vcf.write(line)
		
