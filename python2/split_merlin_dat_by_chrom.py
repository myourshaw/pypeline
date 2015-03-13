#!/usr/bin/env python

import sys
import os

out_prefix = "/Users/myourshaw/lab/ADHD/Paper/link/merlin_2011-03-15/merlin.dat.chrom+ord.all_variables.dat"
header = ''
that_chrom = ''
#CHROMOSOME	ord	locus_type	variable_marker
for chrom in ('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','X','XY','Y'):
	dat = open(os.path.join(out_prefix,'chr'+chrom+'.dat'),"w")
	for line in open("/Users/myourshaw/lab/ADHD/Paper/link/merlin_2011-03-15/merlin.dat.chrom+ord.all_variables.txt"):
		if line.startswith("#") or line.strip() == '': continue
		fields = line.rstrip().split("\t")
		chr = fields[0]
		ord = fields[1]
		locus_type = fields[2]
		variable_marker = fields[3]
		if (locus_type == 'M' and chr != chrom): locus_type = 'S2'
		if ((chrom == "X" and (chr == "X" or chr == "")) or (chrom != "X" and chr != "X")):
			dat.write("%s\t%s\n" % (locus_type,variable_marker))
			