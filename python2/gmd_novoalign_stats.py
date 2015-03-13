#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup

f = os.path.join(os.path.dirname(os.path.abspath(__file__)),"foo")
foo = sys.path
bar = foo[0]

parser = OptionParser()
(options, args) = parser.parse_args()

read_stats_re = re.compile(r"#\s+(?P<key>Paired Reads|Pairs Aligned|Read Sequences|Aligned|Unique Alignment|Gapped Alignment):\s*(?P<value>\d+)")
fragment_stats_re = re.compile(r"#\s+(?P<from>\d+)\s+(?P<to>\d+)\s+(?P<count>\d+)$")
mean_sd_re = re.compile(r"#\s+Mean\s+(?P<mean>[0-9.-eE]+),\s+Std Dev\s+(?P<sd>[0-9.-eE]+)")

#/Volumes/raid/gmd/qout/novoalign_SOLiD_ucla153_20100603_1_GMD18_GMD18_F3.000000.sam.e3023252 /Volumes/raid/gmd/qout/novoalign_SOLiD_ucla153_20100603_1_GMD18_GMD18_F3.000001.sam.e3023255 /Volumes/raid/gmd/qout/novoalign_SOLiD_ucla153_20101102_PE_GMD_GMD98_F5-P2.000054.sam.e3026009 /Volumes/raid/gmd/qout/novoalign_SOLiD_ucla153_20101102_PE_GMD_GMD98_F5-P2.000055.sam.e3026012
read_stats = open("/Volumes/raid/gmd/qout/novoalign_read_stats.txt","w")
read_stats.write("sample\tslice\tkey\tvalue\n")
fragment_stats = open("/Volumes/raid/gmd/qout/novoalign_fragment_stats.txt","w")
fragment_stats.write("sample\tslice\tfrom\tto\tcount\n")
for file in args:
	file = os.path.realpath(file)
	(dir,base) = os.path.split(file) #/Volumes/raid/gmd/qout,novoalign_SOLiD_ucla153_20100603_1_GMD18_GMD18_F3.000000.sam.e3023252
	(name,ext)=os.path.splitext(base) #novoalign_SOLiD_ucla153_20100603_1_GMD18_GMD18_F3.000000.sam,.e3023252
	slice = name.split(".")[1]
	name_split = name.split("_")
	sample = name_split[len(name_split)-2]
	for line in open(file):
		line = line.strip()
		match = read_stats_re.match(line)
		if match != None:
			read_stats.write("%s\t%s\t%s\t%s\n" % (sample,slice,match.group("key"),match.group("value")))
		else:
			match = fragment_stats_re.match(line)
			if match != None:
				fragment_stats.write("%s\t%s\t%u\t%u\t%u\n" % (sample,slice,int(match.group("from")),int(match.group("to")),int(match.group("count"))))
			else:
				match = mean_sd_re.match(line)
				if match != None:
					read_stats.write("%s\t%s\t%s\t%s\n" % (sample,slice,"Mean Fragment Length",match.group("mean")))
					read_stats.write("%s\t%s\t%s\t%s\n" % (sample,slice,"Std Dev Fragment Length",match.group("sd")))
