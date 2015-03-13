#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import gzip
import my

class MergeHeaderFilesError(Exception): pass

#-i /scratch0/tmp/myourshaw/gmd/vcfs/tmp_Vcqu14/gmd32.analysis_ready.vcf.part.*.vax -o /scratch0/tmp/myourshaw/gmd/vcfs/gmd32.analysis_ready.vcf.vax
#--unique_key_data_columns 0 1 3 4 -i /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/tmp_0h8tut/ALL.chr14.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.part.00000000.vax /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/tmp_0h8tut/ALL.chr14.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.part.00000001.vax -o /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/test_tmp/ALL.chr14.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.vax.gz --must_have_data --all_data_records_must_have_all_columns --compress

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'merge files with single set of #-prefixed header lines taken from the first file',
        epilog = 'pypeline.merge_header_files version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', nargs='+',
        help='input file(s)')
    parser.add_argument('--output', '-o', required=True,
        help='output file')
    parser.add_argument('--separator', '-s', default='\t',
        help='field separator (default: tab)')
    parser.add_argument('--compress', action='store_true', default=False,
        help='gzip output (default: False)')
    parser.add_argument('--must_have_data', action='store_true', default=False,
        help='raise error if any file has no lines that do not start with # (default: False)')
    parser.add_argument('--all_data_records_must_have_all_columns', action='store_true', default=False,
        help='raise error if any file has no lines that do not have the same number of tab-separated columns as the last line that starts with a single # (default: False)')
    parser.add_argument('--unique_key_data_columns', nargs='*', default=None, type=int,
        help='list of zero-based indices to data columns for which to print to STDOUT a count of unique values (default: None)')
    args = parser.parse_args()

    skip_header = False
    input_files = my.unglob(args.input,sort_unique=False)
    column_count = 0
    unique_keys = set()
    output_line_count = 0
    output_record_count = 0
    if args.separator != '\t':
        if args.separator == '\\t':
            args.separator = '\t'
        if args.separator == 'tab':
            args.separator = '\t'
        if args.separator == 'space':
            args.separator = ' '
        
    
    my.makedir(os.path.dirname(args.output))
    with gzip.open(args.output, 'w') if args.compress else open(args.output, 'w') as out:
        for input in input_files:
            if not my.file_exists(input):
               raise MergeHeaderFilesError("merge failed: file {} doesn't exist".format(input))
            line_count = 0
            record_count = 0
            for line in my.open_gz_or_text(input):
                line_count += 1
                if line.startswith('#'):
                    if line[1] != '#':
                        column_count = len(line.split(args.separator))
                    if skip_header:
                        continue
                else:
                    record_count += 1
                    output_record_count += 1
                    if args.all_data_records_must_have_all_columns and len(line.split(args.separator)) != column_count:
                        raise MergeHeaderFilesError("merge failed: file {} line {} does not have {} columns".format(input, line_count, column_count))
                    if args.unique_key_data_columns:
                        fields = line.rstrip('\n').split(args.separator)
                        unique_keys.add(tuple([fields[c] for c in args.unique_key_data_columns]))
                output_line_count += 1
                out.write(line)
            skip_header = True
            if args.must_have_data and record_count == 0:
                raise MergeHeaderFilesError("merge failed: file {} has no data".format(input))
    print 'merged {} input files into {} with {} lines and {} data records'.format(len(input_files), args.output, output_line_count, output_record_count)
    if unique_keys and args.unique_key_data_columns:
        print 'having {} unique values in columns {}'.format(len(unique_keys), ','.join([str(c) for c in args.unique_key_data_columns]))

if __name__ == "__main__": sys.exit(main())
