#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

#-i /Volumes/scratch/ucsc/ucsc.GRCh37.refGene
#-i /share/apps/myourshaw/resources/refgene/refgene.b37_ucsc_20130609_header_sorted.txt

parser = argparse.ArgumentParser(
    description = 'create an exome cds + essential splice site interval list from ucsc.refGene table',
    epilog = 'pypeline.refGene2interval_list version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
parser.add_argument('--input', '-i', required=True,
    help='input UCSC refGene file',)
parser.add_argument('--output', '-o',
    help='output interval_list')

args = parser.parse_args()

args.output = args.output if args.output else args.input+'.interval_list'

intervals = {}

with open(args.input) as input:
    for line in input:
        line = line.rstrip('\n')
        if line.startswith('#'):
            cols = line.lstrip('#').split('\t')
            continue
        data = line.split('\t')
        record = {cols[i]:data[i] for i in range(len(cols))}
        exonCount = int(record['exonCount'])
        exonStarts = [int(x) for x in record['exonStarts'].rstrip(',').split(',')]
        exonEnds = [int(x) for x in record['exonEnds'].rstrip(',').split(',')]
        cdsStart = int(record['cdsStart'])
        cdsEnd = int(record['cdsEnd'])
        for i in range(exonCount):
            exonStart = exonStarts[i]
            exonEnd = exonEnds[i]
            if cdsStart > exonEnd or cdsEnd < exonStart or cdsStart == cdsEnd:
                #print '{}\t{}\t{}-{}\t{}-{}'.format(i,record['name'],exonStart,exonEnd,cdsStart,cdsEnd)
                continue
            else:
                intervalStart = max(exonStart,cdsStart)+1 - (2 if exonStart >= cdsStart and i != 0 else 0)
                intervalEnd = min(exonEnd,cdsEnd) + (2 if exonEnd <= cdsEnd and i != exonCount-1 else 0)
                #print '{}\t{}\t{}-{}\t{}-{}\t{}-{}'.format(i,record['name'],exonStart,exonEnd,cdsStart,cdsEnd,intervalStart,intervalEnd)
                intervals[(record['chrom'],intervalStart,intervalEnd)] = (record['name'],record['name2'])
                
collapsed = {}
that_chrom = None
this_chrom = None
that_start = None
that_end = None
names = set()
name2s = set()
with open(args.output,'w') as output:
    for k in sorted(intervals.keys()):
        interval = intervals[k]
        chrom,start,end = k
        if chrom == that_chrom:
            if start > that_end:
                output.write('{}\t{}\t{}\t+\t{};{}\n'.format(that_chrom,that_start,that_end,','.join(names),','.join(name2s)))
                that_start,that_end = start,end
                names = {interval[0]}
                name2s = {interval[1]}
            else:
                that_end = max(end,that_end)
                names.add(interval[0])
                name2s.add(interval[1])
        else:
            if that_chrom:
                output.write('{}\t{}\t{}\t+\t{};{}\n'.format(that_chrom,that_start,that_end,','.join(names),','.join(name2s)))
            that_chrom,that_start,that_end = chrom,start,end
            names = {interval[0]}
            name2s = {interval[1]}

