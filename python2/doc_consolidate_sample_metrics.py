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

class SqlColumnsError(Exception): pass

name = 'pypeline.doc_consolidate_sample_metrics'
version = 1
copyright = '©2013 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()

#-i /scratch1/tmp/myourshaw/mmjj_20130514/metrics/depth_of_coverage/mmjj_20130514.GMD15*.* -d /scratch1/tmp/myourshaw/mmjj_20130514/metrics/depth_of_coverage/consolidated -p mmjj_20130514.sample.bam.depth

def run (input, output_dir, prefix):
    files = my.unglob(input)
    my.makedir(output_dir)
    #dict of [file handle, header]
    file_handles = {'.sample_summary': [None,None], '.sample_gene_summary': [None,None], '.sample_interval_summary': [None,None]}
    #file_handles = {'.sample_summary': [None,None]}
    for f in files:
        (root, ext) = os.path.splitext(f)
        if ext in file_handles.keys() and file_handles[ext][0] == None:
            file_handles[ext][0] = open(os.path.join(output_dir,prefix+ext), 'w')
    for f in files:
        (root, ext) = os.path.splitext(f)
        file_name = os.path.basename (f)
        if ext in file_handles.keys():
            line_count = 0
            for line in open(f):
                line_count += 1
                line = line.rstrip('\n')
                if line_count == 1:
                    if file_handles[ext][1] == None:
                        file_handles[ext][1] = line.split('\t')
                        file_handles[ext][0].write('#file\t'+'\t'.join(file_handles[ext][1])+'\n')
                    else:
                        continue
                else:
                    fields = line.split()
                    if ext == '.sample_summary' and fields[0] == 'Total':
                        continue
                    file_handles[ext][0].write(file_name+'\t'+'\t'.join(fields)+'\n')


def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'consolidate GATK DepthOfCoverage sample metrics',
        epilog = '{} version {} {}'.format(name, version, copyright))
    parser.add_argument('--input', '-i', required = True, nargs='+',
        help='input summary file(s) produced by GATK')
    parser.add_argument('--output_dir', '-d', required = True,
        help='/path/to/output/consolidated')
    parser.add_argument('--prefix', '-p', default = '',
        help='prefix for output files')
    args = parser.parse_args()

    run(input=args.input, output_dir=args.output_dir, prefix=args.prefix)


if __name__ == "__main__": sys.exit(main())
