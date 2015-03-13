#!/usr/bin/env python

import sys
import os

count = 0
out = open("/Users/myourshaw/lab/ADHD/Paper/link/merlin_2011-03-15/merlin-00123456/cat-sim-npl_qtl_regress_LOD1.txt","w")
for line in open("/Users/myourshaw/lab/ADHD/Paper/link/merlin_2011-03-15/merlin-00123456/cat-sim-npl_qtl_regress.txt"):
	fields = line.split("\t")
	if (fields[0] == 'sim' or float(fields[4]) >= 1 ):
		out.write(line)