#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import time
import warnings
import my

class VcfIds2PosError(Exception): pass

#-i /scratch0/tmp/myourshaw/gmd/analysis/vcf/GMD*.gatk.snpFiltered.indelFiltered.recalibrated.vcf /scratch1/tmp/myourshaw/resources//00-All.vcf

def main():

	#command line arguments
	parser = argparse.ArgumentParser(
		description = 'put id=chrom:pos:ref:alt in ID column for consumption by VEP and database',
		epilog = 'pypeline.vcfids2pos version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--input', '-i', nargs='+',
											help='input vcf files')
	args = parser.parse_args()
	
	CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GT = range(10)
	
	for vcf in sorted(my.flatten([glob.glob(v) for v in args.input])):
		print vcf
		with open(vcf) as input:
			with open(vcf.rstrip('.vcf')+'.pos.vcf','w') as output:
				for line in input:
					if line.startswith('#'):
						output.write(line)
					else:
						fields = line.split('\t')
						fields[ID] = '{}={}:{}:{}:{}'.format(fields[ID],fields[CHROM],fields[POS],fields[REF],fields[ALT])
						output.write('\t'.join(fields))
	print 'done'
	

	
if __name__ == "__main__": sys.exit(main())
