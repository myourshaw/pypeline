#!/usr/bin/env python
import sys
import os
import re
#>hg18_ct_UserTrack_3545_chr1:58962[T/C] range=chr1:58962-58962 5'pad=0 3'pad=0 strand=+ repeatMasking=none
#T
#>hg18_ct_UserTrack_3545_chr1:58989[T/C] range=chr1:58989-58989 5'pad=0 3'pad=0 strand=+ repeatMasking=none
#T
out = open("/Volumes/raid/tibhan/hg18.refncbi.txt","w")
for line in open("/Volumes/raid/tibhan/hg18.tibhan.seq"):
	if line.startswith(">"):
		fields = line.strip().split(" ")
		range = fields[1].split("=")[1]
		coord = range.split(":")
		chrom = coord[0]
		pos = coord[1].split("-")[0]
	else:
		ref = line.strip()[0]
		out.write("%s\t%s\t%s\n" % (chrom,pos,ref))