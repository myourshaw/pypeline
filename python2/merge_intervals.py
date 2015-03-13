#!/usr/bin/env python
import sys
import os
import re
#input is bed file sorted by reference genome sort order, start, end
#containing all exons or all cds parts of exons (including 2 base essential splice site within each intron)
#name is comma-separated list of gene models that have the exon
#this script merges overlapping exons and gene model names
# 1	17912	18063	geneid,genscan,sgpGene,xenoRefGene
# 1	17913	18063	geneid,genscan,sgpGene,xenoRefGene	+
interval = open("/Volumes/raid/geneIntervals2/beds/hg19.allGenes14ModelsMerged.cds.interval_list","w")
bed = open("/Volumes/raid/geneIntervals2/beds/hg19.allGenes14ModelsMerged.cds.bed","w")
inputIsBed = True
thisChrom = ""
thisStart = -1
thisEnd = -1
thisModels = set()
for line in open("/Volumes/raid/geneIntervals2/beds/hg19.allGenes14ModelsCollapsed.cds.bed"):
	line = line.strip()
	if line.startswith("#")or len(line) == 0:
		continue
	elif line.startswith("@"):
		interval.write(line+"\n")
	else:
		fields = line.split("\t")
		chrom = fields[0]
		start = int(fields[1]) if inputIsBed else int(fields[1]) -1
		end = int(fields[2])
		name = fields[3]
		if thisChrom != "" and (chrom != thisChrom or start > thisEnd):
			models = ",".join(sorted(thisModels))
			interval.write("%s\t%u\t%u\t%s\t+\n" % (thisChrom,thisStart+1,thisEnd,models))
			bed.write("%s\t%u\t%u\t%s\t1000\t+\n" % (thisChrom,thisStart,thisEnd,models))
			thisChrom = chrom
			thisStart = start
			thisEnd = end
			thisModels = set(name.split(","))
		else:
			thisChrom = chrom
			if thisStart == -1 or start < thisStart:
				thisStart = start
			if end >= thisEnd:
				thisEnd = end
			thisModels = thisModels | set(name.split(","))
models = ",".join(sorted(thisModels))
interval.write("%s\t%u\t%u\t%s\t+\n" % (thisChrom,thisStart+1,thisEnd,models))
bed.write("%s\t%u\t%u\t%s\t1000\t+\n" % (thisChrom,thisStart,thisEnd,models))
