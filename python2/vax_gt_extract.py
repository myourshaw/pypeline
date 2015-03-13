#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import time
import my

class VaxGt2VcfError(Exception): pass

#-i /scratch1/tmp/myourshaw/hlee/vcfs/hlee55gmd28bipolar20ocdtic18.analysis_ready.vcf.vax.gt.txt.gz

#-i /scratch1/tmp/myourshaw/hlee/vcfs/hlee55gmd28bipolar20ocdtic18.analysis_ready.vcf.vax.gt.txt.gz -o /scratch1/tmp/myourshaw/hlee/vcfs/hlee55gmd28bipolar20ocdtic18.analysis_ready.vcf.vax.gt.txt.gz.extracted2.txt -c #CHROM POS ID REF ALT FORMAT SAMPLE GT ALLELE1 PHASE ALLELE2 ZYGOSITY

def main():

    columns_out = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE', 'GT', 'ALLELE1', 'PHASE', 'ALLELE2', 'ZYGOSITY']

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'extract minimal vcf from vax.gt file',
        epilog = 'pypeline.vax_gt_extract version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='input vax.gt file')
    parser.add_argument('--output', '-o',
        help='output file (default:<input>.extracted.txt)')
    parser.add_argument('--headers',
        help='file containing vcf headers')
    parser.add_argument('--columns_to_extract', '-c', nargs='*',
        help='columns to extract (note: quote \'#CHROM\' on command line) (default: {})'.format(columns_out))
    args = parser.parse_args()
    
    if not args.output:
        args.output = args.input+'.extracted.txt'
    if args.columns_to_extract:
        columns_out =  args.columns_to_extract
        
    header = ''
    if args.headers:
        with open(args.headers) as h:
            for line in h:
                if line.startswith('##'):
                    header += line
                else:
                    break
    header += '\t'.join(columns_out)+'\n'
        
    columns_in = {}
    
    print 'extracting [{}] to {}'.format(','.join(columns_out), args.output)
    with my.open_gz_or_text(args.input) as i, open(args.output, 'w') as o:
        o.write(header)
        line_count = 0
        foo_count = 0
        for line in i:
            line_count += 1
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                cols = [c.upper() for c in line.rstrip('\n').split('\t')]
                missing_cols = list(set(columns_out) - set(cols))
                if missing_cols:
                    raise VaxGt2VcfError('input file does not have required columns [{}]'.format(','.join(missing_cols)))
                columns_in = {cols[x]: x for x in range(len(cols))}
                continue
            elif not columns_in:
                raise VaxGt2VcfError('data encountered before column names at line {}'.format(line_count))
            else:
                data_in = line.rstrip('\n').split('\t')
                if not data_in[0]:
                    foo_count += 1
                data_out = [data_in[columns_in[c]] for c in columns_out]
                o.write('\t'.join(data_out)+'\n')
    print 'done, wrote {} lines'.format(line_count)
    print '{} lines missing CHROM'.format(foo_count)


if __name__ == "__main__": sys.exit(main())
