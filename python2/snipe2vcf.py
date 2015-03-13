#!/usr/bin/env python
import sys,os
import re
import traceback
from optparse import OptionParser, OptionGroup
import glob
import time,datetime

# Init cmd-line args
description = '''
Converts somatic sniper output to vcf.
'''
parser = OptionParser(description=description,usage="%prog -d <output dir> -i <snipe results file> -o <vcf output file> -s <sample ID>", version="%prog 1.0.beta.1")
parser.add_option('-i', '--input', action='store', type='string', dest='input', help='Input file = snipe results file [default = STDIN')
parser.add_option('-o', '--output', action='store', type='string', dest='output', help='Output vcf file [default: = STDOUT')
parser.add_option('-s', '--sample', action='store', type='string', dest='sample', default='sample1', help='Sample file [default: sample1]')
parser.add_option('-q', '--quiet', action='store_true', dest='quiet', default=False, help='Quiet mode [default: off]' )
(options, args) = parser.parse_args()

def makedir(dir):
	try:
	  os.makedirs(dir)
	except OSError:
		if os.path.isdir(dir):
			pass
		else:
			raise

def get_files(names):
	return (file for file in names if os.path.isfile(file))

def is_number(s):
	try:
		n = float(s)
		if (math.isinf(n) or math.isnan(n)): return False
		return True
	except ValueError:
		return False

def out(file,string):
	if (not string.endswith('\n')): string += '\n'
	if (not options.quiet):
		print(string)
	file.write(string)

IUPAC = dict(A=('A'),C=('C'),G=('G'),T=('T'),U=('U'),R=('A','G'),Y=('C','T'),S=('C','G'),W=('A','T'),K=('G','T'),M=('A','C'),B=('C','G','T'),D=('A','G','T'),H=('A','C','T'),V=('A','C','G'),N=('A','C','G','T'))

if options.input == None: snipe_file = sys.stdin
else: snipe_file = open(options.input,'r')
if options.output == None:
	options.quiet = True
	vcf_file = sys.stdout
else:
	makedir(os.path.dirname(options.output))
	vcf_file = open(options.output,'w')

CHR=0
POS=1
REF=2
TUMOR_GT=3
SOMATIC_SCORE=4
TUMOR_CONSENSUS_QUAL=5
TUMOR_SNV_QUAL=6
TUMOR_MAP_QUAL=7
DEPTH_TUMOR=8
DEPTH_NORMAL=9
VCF_HEADER=('''##fileformat=VCFv4.0
##fileDate=%s
##source=somaticsniper
##reference=NCBI37
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Somatic score">
##FORMAT=<ID=TCQ,Number=1,Type=Integer,Description="Tumor consensus quality">
##FORMAT=<ID=TSQ,Number=1,Type=Integer,Description="Tumor SNV quality">
##FORMAT=<ID=TMQ,Number=1,Type=Integer,Description="Tumor map quality">
##FORMAT=<ID=DT,Number=1,Type=Integer,Description="Depth tumor">
##FORMAT=<ID=DN,Number=1,Type=Integer,Description="Depth normal">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s
''' % (time.strftime('%Y%m%d',time.localtime()),options.sample))
out(vcf_file, VCF_HEADER)
line_number=0
for line in snipe_file:
	line_number+=1
	fields=line.strip().split('\t')
	ref = fields[REF]
	tgt = fields[TUMOR_GT]
	alt = list(IUPAC[tgt])
	if len(alt) == 1:
		if tgt == ref: gt = '0/0'
		else: gt = '1/1'
	elif len(alt) == 2:
		if ref in alt:
			alt.remove(ref)
			gt= '0/1'
		else:
			gt= '1/2'		
	else:
		sys.stderr.write ('line %u genotype %s has more than two alleles' % (line_number, alt))
		#gt is ref/alt = 0/1
	out(vcf_file, '%s\t%u\t.\t%s\t%s\t.\t.\t.\tGT:SS:TCQ:TSQ:TMQ:DT:DN\t%s:%u:%u:%u:%u:%u:%u\n' % (
		fields[CHR],
		int(fields[POS]),
		fields[REF],
		','.join(alt),
		gt,
		int(fields[SOMATIC_SCORE]),
		int(fields[TUMOR_CONSENSUS_QUAL]),
		int(fields[TUMOR_SNV_QUAL]),
		int(fields[TUMOR_MAP_QUAL]),
		int(fields[DEPTH_TUMOR]),
		int(fields[DEPTH_NORMAL])))
	
	