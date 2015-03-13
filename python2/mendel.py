#!/usr/bin/env python

import sys
import os
from optparse import OptionParser, OptionGroup

description = """
This script calculates Mendellian inconsistency statistics from a VCF file.
The script expects three genotypes: father, mother, child in that order.
The INFO field of the VCF file can contain a tag "ALN=<alignment option>"
"""

parser = OptionParser(description=description)
(options, args) = parser.parse_args()
file = args[0]

proband_sex = 'm'
line_count = 0
#VCF fields
CHROM = 0
POS = 1
ID = 2
REF = 3
ALT = 4
QUAL = 5
FILTER = 6
INFO = 7
FORMAT = 8
FATHER = 9
MOTHER = 10
PROBAND = 11
#totals is a key, value dictionary
#key is alignment name, e.g. "bwa_snp"
#value is a list of counters
totals = {}
ok = 0
error = 1
indeterminate = 2
total = 3
ok_known = 4
error_known = 5
indeterminate_known = 6
total_known = 7
ok_novel = 8
error_novel = 9
indeterminate_novel = 10
total_novel = 11
missing_dad_geno = 12
missing_mom_geno = 13
missing_kid_geno = 14
missing_any_geno = 15
#valid genotypes (dad,mom,kid)
valid = (
	'0/0,0/0,0/0',
	'0/0,0/1,0/0',
	'0/0,0/1,0/1',
	'0/0,1/1,0/1',
	'0/1,0/0,0/0',
	'0/1,0/0,0/1',
	'0/1,0/1,0/0',
	'0/1,0/1,0/1',
	'0/1,0/1,1/1',
	'0/1,1/1,0/1',
	'0/1,1/1,1/1',
	'1/1,0/0,0/1',
	'1/1,0/1,0/1',
	'1/1,0/1,1/1',
	'1/1,1/1,1/1'
	)
#genotypes that are invalid even though a parent is missing
invalid_missing = (
	'0/0,./.,1/1',
	'1/1,./.,0/0',
	'./.,0/0,1/1',
	'./.,1/1,0/0',
)

for line in open (file):
	line_count +=1
	line = line.strip()
	if line.startswith("#") or line == "":
		continue
	#get the fields of the VCF file
	fields = line.split("\t")
	chrom = fields[CHROM].upper()
	pos = int(fields[POS])
	#skip Y and funny contigs
	if chrom == 'Y' or len(chrom) > 2:
		continue
	#put the INFO tags in a list, then check for some flags
	info = fields[INFO].split(';')
	PASS = True if 'PASS' in info else False
	REFSEQ = True if 'REFSEQ' in info else False
	CONCORDANT = True if 'CONCORDANT' in info else False
	NOVEL = True if 'NOVEL' in info else False
	#make a dictionary of the key-value pairs in INFO
	tags = {}
	for tag in info:
		kv = tag.split('=')
		if len(kv) == 2:
			tags[kv[0]] = kv[1]
	DBSNP132 = tags['DBSNP132'] if 'DBSNP132' in tags else '.'
	VC = tags['VC'] if 'VC' in tags else '?'
	ALN = tags['ALN'] if 'ALN' in tags else '.'
	#skip records unless they have a PASS tag and a REFSEQ tag in the INFO field
	if not PASS or not REFSEQ:
		continue
	#if this alignment hasn't been seen before, make an entry in the totals dictionary with counters initially set to zero
	if ALN not in totals:
		totals[ALN] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	#treat known as CONCORDANT and novel as not CONCORDANT
	known = CONCORDANT
	#extract the genotypes from the sample fields
	dad_gt = fields[FATHER].split(":")[0]
	mom_gt = fields[MOTHER].split(":")[0]
	kid_gt = fields[PROBAND].split(":")[0]
	#accumulate totals depending on various tags and whether the pedigree is consistent
	totals[ALN][total] += 1
	if known:
		totals[ALN][total_known] += 1
	else:
		totals[ALN][total_novel] += 1
	if chrom == 'X' and (dad_gt == '0/1' or (proband_sex == 'm' and (kid_gt == '0/1' or kid_gt.split('/')[0] not in mom_gt))):
			print ALN,chrom,pos,dad_gt,mom_gt,kid_gt
			totals[ALN][error] += 1
			if known:
				totals[ALN][error_known] += 1
			else:
				totals[ALN][error_novel] += 1
	elif chrom == 'X' and proband_sex == 'm':
		totals[ALN][ok] += 1
		if known:
			totals[ALN][ok_known] += 1
		else:
			totals[ALN][ok_novel] += 1
	else:
		test = ",".join((dad_gt,mom_gt,kid_gt))
		if test in valid:
			#mendellian consistent and all samples have genotypes
			totals[ALN][ok] += 1
			if known:
				totals[ALN][ok_known] += 1
			else:
				totals[ALN][ok_novel] += 1
		else:
			#mendellian inconsistent or missing genotype(s)
			if'.' not in test or test in invalid_missing:
				#mendellian inconsistent
				print ALN,chrom,pos,dad_gt,mom_gt,kid_gt
				totals[ALN][error] += 1
				if known:
					totals[ALN][error_known] += 1
				else:
					totals[ALN][error_novel] += 1
			else:
				#indeterminate: missing genotype(s) and not inconsistent
				totals[ALN][indeterminate] += 1
				if known:
					totals[ALN][indeterminate_known] += 1
				else:
					totals[ALN][indeterminate_novel] += 1
	#count missing genotypes
	if dad_gt == './.':
		totals[ALN][missing_dad_geno] += 1
	if mom_gt == './.':
		totals[ALN][missing_mom_geno] += 1
	if kid_gt == './.':
		totals[ALN][missing_kid_geno] += 1
	if dad_gt == './.' or mom_gt == './.' or  kid_gt == './.':
		totals[ALN][missing_any_geno] += 1
#print results to STDOUT
print "Alignment\tOK\terror\tindeterminate\ttotal\tok_known\terror_known\tindeterminate_known\ttotal_known\tok_novel\terror_novel\tindeterminate_novel\ttotal_novel\tmissing_dad_geno\tmissing_mom_geno\tmissing_kid_geno\tmissing_any_geno"
for k,v in sorted(totals.iteritems()):
    print "%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u" % \
    (k,v[ok],v[error],v[indeterminate],v[total],v[ok_known],v[error_known],v[indeterminate_known],v[total_known],v[ok_novel],v[error_novel],v[indeterminate_novel],v[total_novel],v[missing_dad_geno],v[missing_mom_geno],v[missing_kid_geno],v[missing_any_geno])
