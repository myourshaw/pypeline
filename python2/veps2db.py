#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import time
import my

class Veps2DbError(Exception): pass

#-i /scratch0/tmp/myourshaw/gmd/analysis/vcf/*.gatk.snpFiltered.indelFiltered.recalibrated.pos.vcf.vep -o /scratch0/tmp/myourshaw/gmd/analysis/vcf/gmddb2/GMD

def main():

	#command line arguments
	parser = argparse.ArgumentParser(
		description = 'merge multiple Variant Effect Predictor 2.1 files into a txt file for database entry',
		epilog = 'pypeline.veps2db version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--input', '-i', nargs='+',
											help='input vep files')
	parser.add_argument('--output', '-o', required=True,
											help='path/prefix for output files')
	args = parser.parse_args()
	
	Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,Extra = range(14)
	
	#regular expressions for vep header lines and data fields
	extra_re = re.compile(r'##\s+(?P<col>\S+)\s+:\s+(?P<description>.+)', re.I)
	fileformat_re = re.compile(r'##\s+ENSEMBL VARIANT EFFECT PREDICTOR\s+v(?P<version>[0-9.]+)', re.I)
	location_re = re.compile(r'(?P<chrom>\S+):(?P<chromStart>\d+)-?(?P<chromEnd>\d+)?', re.I)
	
	#get all possible extra columns
	extra_cols = {}
	in_cols = False
	veps = sorted(my.flatten([glob.glob(v) for v in args.input]))
	for vep in veps:
		line_count = 0
		for line in open(vep, 'r'):
			line_count += 1
			line = line.rstrip('\n')
			if line.startswith('#Uploaded_variation'):
				break
			elif line.startswith('## Extra column keys:'):
				in_cols = True
				continue
			elif in_cols and line.startswith('##'):
				m = extra_re.match(line)
				if m:
					extra_cols[m.group('col')] = m.group('description')
	extra_cols = sorted(extra_cols.keys())
	
	
	out_base = args.output.rstrip('.txt')
	my.makedir(os.path.dirname(out_base))
	vep_out = open(out_base+'vep.db.txt', 'w')
	vep_out.write('VEP\tTIME\tCHROM\tPOS\tREF\tALT\tUploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\t{}\n'.format('\t'.join(extra_cols)))
	consequence_out = open(out_base+'vep.db.consequence.txt', 'w')
	consequence_out.write('VEP\tTIME\tCHROM\tPOS\tREF\tALT\tUploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\t{}\n'.format('\t'.join(extra_cols)))

	for vep in veps:
		print vep
		vep_time = time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(os.path.getmtime(vep)))
		format_ok = False
		line_count = 0
		for line in open(vep, 'r'):
			line_count += 1
			line = line.rstrip('\n')
			if not bool(line.strip()):
				continue
			if not format_ok:
				m = fileformat_re.match(line)
				if not m:
					raise Veps2DbError('missing fileformat in {} at line {}'.format(vep, line_count))
				elif not float(m.group('version')) >= 2.1:
					raise Veps2DbError('obsolete fileformat in {} at line {}, {}'.format(vep, line_count, m.group('version')))
				else:
					format_ok = True
					continue
			elif line.startswith('#'):
				continue
			else:
				fields = line.split('\t')
				values = fields[:Extra]
				uploaded_variation, location, allele, gene, feature, feature_type, consequence, cdna_position, cds_position, protein_position, amino_acids, codons, existing_variation = values
				m = location_re.match(location)
				chrom, chromStart, chromEnd = m.groups(m.group('chromStart')) if m else ('','','')
				uv = uploaded_variation.split('=')
				uploaded_variation = uv[0]
				if len(uv) > 1:
					CHROM,POS,REF,ALT = uv[1].split(':')
					alt_list = ALT.split(',')
				extra = fields[Extra] if len(fields) > Extra else ''
				extras = {c: '' for c in extra_cols}
				if extra and extra != '-':
					exs = [x.split('=') for x in extra.split(';')]
					for e in exs:
						if len(e)>1:
							extras[e[0]] = e[1]
				extra_values = [extras[e] for e in sorted(extras.keys())]
				for this_alt in alt_list:
					vep_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(vep, vep_time, CHROM, POS, REF, this_alt, '\t'.join(fields[:Extra]), '\t'.join(extra_values)))
					consequence_list = consequence.split(',')
					for this_consequence in consequence_list:
						consequence_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(vep, vep_time, CHROM, POS, REF, this_alt, '\t'.join(fields[:Consequence]), this_consequence, '\t'.join(fields[cDNA_position:Extra]), '\t'.join(extra_values)))

	vep_out.close()
	consequence_out.close()
	print 'done'
	
	
if __name__ == "__main__": sys.exit(main())
