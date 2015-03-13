#!/usr/bin/env python

out = open('refgene_miscounts.txt','w')
out.write('line\texonCount\ttxCount\tcdsCount\tphaseCount\tbin\tid\tgene\n')
lineN = 0
for line in open('/share/apps/myourshaw/resources/refgene.b37.sorted.txt'):
	lineN += 1
	fields = line.rstrip('\n').split('\t')
	exonCount = int(fields[8])
	txStarts = fields[9].rstrip(',').split(',')
	cdsStarts = fields[10].rstrip(',').split(',')
	phases = fields[15].rstrip(',').split(',')
	txCount = len(txStarts)
	cdsCount = len(cdsStarts)
	phaseCount = len(phases)
	if txCount != exonCount or cdsCount != exonCount or phaseCount != exonCount:
		out.write("%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\n" % (lineN,exonCount,txCount,cdsCount,phaseCount,int(fields[0]),fields[1],fields[12]))
