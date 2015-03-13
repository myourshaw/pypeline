#!/usr/bin/env python
import sys
import os
import re
from string import maketrans   # Required to call maketrans function.
bed = open("/Volumes/raid/tibhan/tibhan.hg18.bed","w")
vcf = open("/Volumes/raid/tibhan/tibhan.hg18.vcf","w")
vcf.write("""##fileformat=VCFv4.0
##fileDate=201011203
##source=Tibetan/Han exomes majorminor.txt from Emilia Huerta <emilia.huertasanchez@gmail.com>
##reference=hg18
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
""")
table = maketrans("0123","ACGT")
for line in open("/Volumes/raid/tibhan/majorminor.txt"):
	if line.startswith("chromosome"):
		continue
	(chrom,position,major,minor) = line.strip().split("\t")
	ref = major.translate(table)
	alt = minor.translate(table)
	vcf.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n" % (chrom,position,ref,alt))
	bed.write("%s\t%u\t%u\t%s\n" % (chrom, int(position)-1, int(position), chrom+":"+position+"["+ref+"/"+alt+"]"))