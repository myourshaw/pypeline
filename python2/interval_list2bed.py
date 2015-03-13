#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from ConfigParser import SafeConfigParser
import my

class IntervalList2BedError(Exception): pass

#-i /scratch1/tmp/myourshaw/resources/intervals/TruSeq/TruSeq_exome_targeted_regions.txt.interval_list

def main():
   
    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'convert an interval_list to  a BED file',
        epilog = 'pypeline.verify_OK_file version 1.0β1 ©2011-2013 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True, nargs='+',
        help='input file(s)')
    parser.add_argument('--header',
        help='header line(s) for BED file(s)')
    args = parser.parse_args()

    interval_lists=my.unglob(args.input)
    for interval_list in interval_lists:
        with open(interval_list) as il, open(interval_list.rstrip('.interval_list')+'.bed', 'w') as bed:
            if args.header:
                if not args.header.endswith('\n'):
                    args.header += '\n'
                bed.write(args.header)
            for line in il:
                if line.strip() == '' or line.startswith('@'):
                    continue
                chrom, start,end, strand, name = line.rstrip('\n').split('\t')[:5]
                bedline = '\t'.join([chrom, start,end, name, '1000', strand])+'\n'
                bed.write(bedline)

if __name__ == "__main__": sys.exit(main())
