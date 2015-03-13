#!/usr/bin/env python

#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
import glob
import array
import numpy

usage = "usage: %prog [options] <list of pileup files in consensus format>"
parser = OptionParser(usage=usage)
parser.add_option("-o", "--output",default="STDOUT", metavar="OUTPUTFILE", dest="out", help="write results to OUTPUTFILE [%default]")
(options, args) = parser.parse_args()
outputFile = options.out

def get_files(names):
	return (file for file in names if os.path.isfile(file))

header = "file\tcoverage\thom_ref\thom_alt\thom_%_alt\thet_ref\thet_alt\thet_%_alt"
if options.out != "STDOUT":
	saveout = sys.stdout                                     
	out = open(options.out, 'w')                             
	sys.stdout = out                                       
print header

acgt = re.compile(r"[ACGT]",re.I)
(hom_ref_count,hom_alt_count,het_ref_count,het_alt_count) = range(4)
iupac = dict(A=('A'),C=('C'),G=('G'),T=('T'),U=('U'),R=('A','G'),Y=('C','T'),S=('C','G'),W=('A','T'),K=('G','T'),M=('A','C'),B=('C','G','T'),D=('A','G','T'),H=('A','C','T'),V=('A','C','G'),N=('A','C','G','T'))
counters_merge = numpy.zeros((101,4))

for arg in args:
	for file in get_files(glob.iglob(arg)):
		counters = numpy.zeros((101,4))
		for line in open(file):
			fields = line.strip().split()
			if fields[2] == "*" or fields[3] == "N" or fields[7] == "0": #ignore indel lines, uncallable consensus, and no reads
				continue
			(chrom,pos,ref,consensus,consensus_quality,SNP_quality,RMS_mapping_quality,read_bases,read_qualities,alignment_qualities) = fields
			coverage = int(read_bases)
			ref_count = read_qualities.count(".") + read_qualities.count(",")
			alt_count = coverage - ref_count
			coverage = min(coverage,101)-1
			if acgt.match(consensus) != None: #homozygous
				counters[coverage,hom_ref_count] += ref_count
				counters[coverage,hom_alt_count] += alt_count
				counters_merge[coverage,hom_ref_count] += ref_count
				counters_merge[coverage,hom_alt_count] += alt_count
			else: #IUPAC code for heterozygous
				counters[coverage,het_ref_count] += ref_count
				counters[coverage,het_alt_count] += alt_count
				counters_merge[coverage,het_ref_count] += ref_count
				counters_merge[coverage,het_alt_count] += alt_count
		file_name = os.path.basename(file)
		for i in range(101):
			hom_ref = counters[i,hom_ref_count]
			hom_alt = counters[i,hom_alt_count]
			hom_pct_alt = hom_alt/(hom_ref+hom_alt)
			het_ref = counters[i,het_ref_count]
			het_alt = counters[i,het_alt_count]
			het_pct_alt = het_alt/(het_ref+het_alt)
			line_out = "%s\t%u\t%u\t%u\t%f\t%u\t%u\t%f"%(file_name,i+1,hom_ref,hom_alt,hom_pct_alt,het_ref,het_alt,het_pct_alt)
			print (line_out)
file_name = "*"
for i in range(101):
	hom_ref = counters_merge[i,hom_ref_count]
	hom_alt = counters_merge[i,hom_alt_count]
	hom_pct_alt = hom_alt/(hom_ref+hom_alt)
	het_ref = counters_merge[i,het_ref_count]
	het_alt = counters_merge[i,het_alt_count]
	het_pct_alt = het_alt/(het_ref+het_alt)
	line_out = "%s\t%u\t%u\t%u\t%f\t%u\t%u\t%f"%(file_name,i+1,hom_ref,hom_alt,hom_pct_alt,het_ref,het_alt,het_pct_alt)
	print (line_out)
if options.out != "STDOUT":
	sys.stdout = saveout                                     
	out.close()                                            
