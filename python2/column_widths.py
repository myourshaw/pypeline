#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import time
import warnings
import copy
import my

# DEPRECATED: use sql_columns.py instead

class Vcfs2ColumnWidthsError(Exception): pass

#dbsnp135
#-i /scratch1/tmp/myourshaw/resources/dbsnp135/00-All.vcf

#-i /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulationNoGeno/YRI-1412-Y-nogeno.vcf.gz
#-i /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulationNoGeno/ByPopulationNoGeno.variants.txt
#-i /scratch1/tmp/myourshaw/resources/dbsnp135/00-All.vcf.flat.txt
#--skip ## --header_line #CHROM -i /scratch0/tmp/myourshaw/1000genomes/phase1_integrated_release_version3_20120430/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
#-i /scratch0/tmp/myourshaw/hgnc/hgnc_xref_20120902.txt --header_line 1

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'calculate maximum column widths for tab-delimited file',
        epilog = 'pypeline.column_widths version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required = True, nargs='+',
        help='input  file(s)[.gz]')
    parser.add_argument('--skip',
        help='skip lines that start with character(s)')
    parser.add_argument('--header_line', default = '#',
        help='1-based row number of column names or # if the first row that starts with a single # character')
    args = parser.parse_args()
    
    column_widths = {}
    column_names = {}
    line_count = 0
    header_line_number = int(args.header_line) if my.is_int(args.header_line) and int(args.header_line) > 0 else None
    header_chars = args.header_line if not my.is_int(args.header_line) else None
    header_found = False if header_line_number or header_chars else True
    
    files = my.unglob(args.input)
    for file in files:
        print file
        header_found = False if header_line_number or header_chars else True
        with my.open_gz_or_text(file) as input, open(file+'.column_widths.txt', 'w') as output:
            for line in input:
                line_count += 1
                line = line.rstrip('\n')
                if not bool(line.strip()):
                    continue
                if not header_found:
                    if (header_line_number and header_line_number == line_count) or (header_chars and line.startswith(header_chars)):
                        if header_chars == '#' and line[1] == '#':
                            continue
                        fields = line.split('\t')
                        column_names = {i:fields[i] for i in range(len(fields))}
                        header_found = True
                    continue
                elif args.skip and line.startswith(args.skip):
                    continue
                else: #data line
                    fields = line.split('\t')
                    indices = [i for i in range(len(fields))]
                    column_widths = {i: max(len(fields[i]),column_widths.setdefault(i,0)) for i in indices}
                    pass
            output.write('col_number\tcol_name\tcol_width\n')
            output.write(''.join(['{}\t{}\t{}\n'.format(col+1, column_names[col], column_widths[col]) for col in sorted(column_widths)]))


	
if __name__ == "__main__": sys.exit(main())
