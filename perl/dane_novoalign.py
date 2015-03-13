#!/usr/bin/env python

import sys
import os
import re

#try:
for line in open("/Volumes/raid/dane/SRX012142/SRR029011.fastq","r"):
	# @SRR029011.1 FC427EKAAXX:6:1:0:85
	sra_header_regex = re.compile(r"@SRR(?P<run>\d+)\.(?P<read>\d+)\s(?P<flowcell>\S+):(?P<lane>\d+):(?P<tile>\d+):(?P<x>\d+):(?P<y>\d+)")
	match = sra_header_regex.search(line)
	if match != None:
		run = int(match.group("run"))
		flowcell = match.group("flowcell")
		lane = match.group("lane")
		sample = run-22463
		break
# '@RG\tID:readGroupId\tSM:sample\tLB:library\tDS:description\tPU:platformUnit(lane|slide)\tPI:predictedMedianInsertSize\tCN:sequencingCenter\tDT:runDate(ISO8601)\tPL:platform(ILLUMINA|SOLID)'
#  @RG\tID:SRR029011\tSM:SRS006548\tPU:BGI-FC427EKAAXX-6\tCN:BGI\tPL:ILLUMINA
rg_header = "@RG\\tID:SRR{0:06}\\tSM:SRS{1:06}\\tPU:BGI-{2}-{3}\\tCN:BGI\\tPL:ILLUMINA".format(run,sample,flowcell,lane)
print rg_header