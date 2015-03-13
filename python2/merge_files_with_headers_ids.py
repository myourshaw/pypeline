#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import gzip
from collections import OrderedDict
import my

class MergeFilesWithHeadersIDsError(Exception): pass

#--id_column_header sample --header_starts_with tracking_id -i "HB8n_P5_EGF_10_High,/scratch1/tmp/myourshaw/rnaseq_20140227/cufflinks/cufflinks_no-upper-quartile-norm/HB8n_P5_EGF_10_High.CAGATC.2/isoforms.fpkm_tracking" "HB8n_P5_EGF_10_Low,/scratch1/tmp/myourshaw/rnaseq_20140227/cufflinks/cufflinks_no-upper-quartile-norm/HB8n_P5_EGF_10_Low.GTGAAA.2/isoforms.fpkm_tracking" -o /scratch1/tmp/myourshaw/rnaseq_20140227/cufflinks/cufflinks_no-upper-quartile-norm/merged_HB8_isoforms.fpkm_tracking_test

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'merge files that have headers and add ids',
        epilog = 'pypeline.merge_files_with_headers_ids version 1.0β1 ©2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', nargs='+',
        help='list of id,file1,file2,...')
    parser.add_argument('--output', '-o', required=True,
        help='output file')
    parser.add_argument('--id_column_header', default='ID',
        help='field name of id column that will be added (default: ID)')
    parser.add_argument('--id_last', action='store_true', default=False,
        help='put id as last column (default: id in first column)')
    parser.add_argument('--header_starts_with', default='#',
        help='unique string that when first encountered marks the header (default: #)')
    parser.add_argument('--skip', default=None,
        help='skip lines that start with this string (default: None)')
    parser.add_argument('--separator', '-s', default='\t',
        help='field separator (default: tab)')
    parser.add_argument('--compress', action='store_true', default=False,
        help='gzip output (default: False)')
    parser.add_argument('--must_have_data', action='store_false', default=True,
        help='raise error if any file has no lines that do not start with header_starts_with string (default: True)')
    parser.add_argument('--all_data_records_must_have_all_columns', action='store_false', default=True,
        help='raise error if any file has no lines that do not have the same number of tab-separated columns as the last line that starts with a single # (default: False)')
    args = parser.parse_args()

    skip_header = False
    header_written = False
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
    
    #OrderedDict to preserve order of samples and their files
    ids = [i.split(',')[0] for i in args.input]
    file_lists = [i.split(',')[1:] for i in args.input]
    #key=file_path, value=sample_id
    id_files = OrderedDict(zip(ids, file_lists))
    input_file_count = 0
    
    with gzip.open(args.output, 'w') if args.compress else open(args.output, 'w') as out:
        for id in id_files.keys():
            print "ID: {}".format(id)
            input_files = my.unglob(id_files[id], sort_unique=False)
        
            for input in input_files:
                input_file_count += 1
                print "\t{}:".format(input)
                if not my.file_exists(input):
                   raise MergeFilesWithHeadersIDsError("merge failed: file {} doesn't exist".format(input))
                line_count = 0
                record_count = 0
                header_found = False
                for line in my.open_gz_or_text(input):
                    line_count += 1
                    if (args.skip and line.startswith(args.skip)) or line.strip() == '':
                        continue
                    elif line.startswith(args.header_starts_with):
                        header_found = True
                        if skip_header:
                            continue
                        else:
                            skip_header = True
                            column_count = len(line.split(args.separator))
                    else:
                        record_count += 1
                        output_record_count += 1
                        if args.all_data_records_must_have_all_columns and len(line.split(args.separator)) != column_count:
                            raise MergeFilesWithHeadersIDsError("merge failed: file {} line {} does not have {} columns".format(input, line_count, column_count))
                    output_line_count += 1
                    if args.id_last:
                        line_out = line+args.separator+(args.id_column_header if output_line_count == 1 else id)
                    else:
                        line_out = (args.id_column_header if output_line_count == 1 else id)+args.separator+line
                    out.write(line_out)
                    header_written = True
                if args.must_have_data and record_count == 0:
                    raise MergeFilesWithHeadersIDsError("merge failed: file {} has no data".format(input))
                print '\t\t{} lines, {} data records'.format(line_count, record_count)
    print '\ninput files: {}\n\noutput file: {}\n\t{} lines, {} data records'.format(input_file_count, args.output, output_line_count, output_record_count)
    print '\ndone'
    
if __name__ == "__main__": sys.exit(main())
