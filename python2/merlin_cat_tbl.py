#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
import glob
import math

#-o /Users/myourshaw/lab/ADHD/Paper/link/merlin_2011-03-15/merlin_20110315.regress_vc_covar.cat.txt /Volumes/rocks5.scratch1/link/merlin_2011-03-15/vc_covar/chr*/chr*-vc-chr*.tbl /Volumes/rocks5.scratch1/link/merlin_2011-03-15/regress_covar/chr*/chr*-regress-chr*.tbl
#-o /Users/myourshaw/lab/ADHD/Paper/link/merlin_2011-03-15/merlin_20110315.regress_vc.cat.txt /Volumes/rocks5.scratch1/link/merlin_2011-03-15/vc/chr*/chr*-vc-chr*.tbl /Volumes/rocks5.scratch1/link/merlin_2011-03-15/regress/chr*/chr*-regress-chr*.tbl

usage = "usage: %prog [options] <list of merlin *.tbl files>"
#*-nonparametric.tbl
#CHR	POS	LABEL	ANALYSIS	ZSCORE	DELTA	LOD	PVALUE
#na	na	min	adhd [ALL]	-22.612	-0.271	-40.132	1
#na	na	max	adhd [ALL]	28.276	0.707	100.475	6.204e-103
#22	0.082	rs7288876	adhd [ALL]	1.243	0.111	0.528	0.0594

#without covariates
#chr*-regress-chr*.tbl
#CHR	POS	PHENOTYPE	H2	SD	INFO	LOD	PVALUE
#22	rs7288876	ADHD_sum	0.000	0.125	48.032	0.000	0.5
#with covariates:
#chr*-regress-chr*.tbl
#CHR	POS	PHENOTYPE	H2	SD	INFO	LOD	PVALUE
#22	rs7288876	Trait: ADHD_sum	0.000	0.105	49.754	0.000	0.5

#chr*-vc-chr*.tbl
#CHR	POS	LABEL	TRAIT	H2	LOD	PVALUE
#22	0.082	rs7288876	ADHD_sum	0.000	0.005	0.4374


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

header_npl = "CHR\tPOS\tLABEL\tANALYSIS\tZSCORE\tDELTA\tLOD\tPVALUE\n"
header_regress = "CHR\tPOS\tPHENOTYPE\tH2\tSD\tINFO\tLOD\tPVALUE\n"
header_vc = "CHR\tPOS\tLABEL\tTRAIT\tH2\tLOD\tPVALUE\n"
header_out = "analysis\tvariable\tmarker\tlod\tpvalue"
if options.out != "STDOUT":
	saveout = sys.stdout                                     
	out = open(options.out, 'w')                             
	sys.stdout = out                                       
print header_out

#/Volumes/rocks5.scratch1/link/merlin_2011-03-15/vc_covar/chr01/chr01-vc-chr01.tbl
file_re = re.compile(r".*\/chr(?P<chr>[0-9XYMT]{1,3})-(?P<analysis>\S+)-chr(?P<chrnum>[0-9XYMT]{1,3})\.tbl$",re.I)

for arg in args:
	for file in get_files(glob.iglob(arg)):
		m = file_re.match(file)
		if m == None: raise Exception("Unrecognized file name format: %s" % (file))
		filechr = m.group('chr')
		fileanalysis = m.group('analysis')
		filechrnum = m.group('chrnum')
		linecount = 0
		analysis_type = ''
		for line in open(file):
			linecount += 1
			if (line.strip() == ''): continue
			if analysis_type == '':
				if (line == header_npl): analysis_type = 'n'
				elif (line == header_regress): analysis_type = 'r'
				elif (line == header_vc): analysis_type = 'v'
				else: raise Exception("Unrecognized header: %s in file %s" % (line,file))
				continue
			fields = line.rstrip('\n').split('\t')
			chr = fields[0]
			if(analysis_type == 'n'):
				foo = fields[3].split(' ')
				bar = foo[1].strip(' []')
				if bar == 'ALL': analysis = 'a'
				elif bar == 'Pairs': analysis = 'p'
				elif bar == 'QTL': analysis = 'q'
				else: raise Exception("Unknown nonparametric analysis %s in file %s line %u %s" % (bar,file,linecount,line))
				marker = fields[2]
				variable = foo[0].strip()
				lod = '0' if fields[6] == 'na' else fields[6]
				pvalue = '1' if fields[7] == 'na' else fields[7]
			if(analysis_type == 'r'):
				analysis = 'r'
				marker = fields[1]
				variable = fields[2].split(' ')[1] if fields[2].startswith("Trait") else fields[2]
				lod = '0' if fields[6] == 'na' else fields[6]
				pvalue = '1' if fields[7] == 'na' else fields[7]
			if(analysis_type == 'v'):
				analysis = 'v'
				marker = fields[2]
				variable = fields[3]
				lod = '0' if fields[5] == 'na' else fields[5]
				pvalue = '1' if fields[6] == 'na' else fields[6]
			if (not is_number(lod) or not is_number(pvalue)):
				raise Exception("lod or pvalue not a number in file %s line %u %s", (file,linecount, line))
			print("%s\t%s\t%s\t%s\t%s" % (analysis,variable,marker,lod,pvalue))
