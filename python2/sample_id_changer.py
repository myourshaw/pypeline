#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import my

class GtxSampleIdChangerError(Exception): pass

#-i /scratch1/tmp/myourshaw/mmjj_20130514/vcfs/vax_tmp/mmjj_20130514.hc.analysis_ready.vcf.gtx --old_sample GMD16A --new_sample CIPO16A
#-i /scratch1/tmp/myourshaw/mmjj_20130514/vcfs/vax/mmjj_20130514.hc.analysis_ready.vcf.vex.vcf --old_sample GMD16A --new_sample CIPO16A

#-i /scratch1/tmp/myourshaw/mmjj_20130514/vcfs/mmjj_20130514.ug.recalibrated.*.vcf.gtx --old_sample GMD16A --new_sample CIPO16A
#-i /scratch1/tmp/myourshaw/mmjj_20130514/vcfs/mmjj_20130514.hc.raw.indels.vcf --old_sample GMD16A --new_sample CIPO16A


def main():

    name = 'sample_id_changer'
    version = 1
    copyright = '©2011-2013 Michael Yourshaw all rights reserved'


    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'change sample ids in a vcf or gtx file',
        epilog = '{} version {} {}'.format(name, version, copyright))
    #input VCF parameter
    parser.add_argument('--input', '-i', nargs='+', required=True,
        help='input vcf or gtx file(s); output will be <input>.<old_sample>_<new_sample>.<input_extension>')
    parser.add_argument('--old_sample', required=True,
        help='old sample id')
    parser.add_argument('--new_sample', required=True,
        help='new sample id')
    parser.add_argument('--no_change_column_header', action='store_true', default=False,
        help='do not change sample name in the column header')
    args = parser.parse_args()
    
    for input_file in my.unglob(args.input):
        header_found = False
        info_col, sample_col = -1, -1
        output_file = input_file + '.' + args.old_sample + '_' + args.new_sample + os.path.splitext(input_file)[1]
        with my.open_gz_or_text(input_file) as i, open(output_file, 'w') as o:
            for line in i:
                if not header_found:
                    if line.startswith('#CHROM'):
                        exact_columns=line.rstrip('\n').split('\t')
                        columns = line.rstrip('\n').upper().split('\t')
                        if 'INFO' in columns:
                            info_col = columns.index('INFO')
                        if 'SAMPLE' in columns:
                            sample_col = columns.index('SAMPLE')
                        if (not args.no_change_column_header) and args.old_sample.upper() in columns:
                            exact_columns[columns.index(args.old_sample.upper())] = args.new_sample
                        header_found = True
                        o.write('\t'.join(exact_columns)+'\n')
                    else:
                        o.write(line)
                else:
                    fields = line.rstrip().split('\t')
                    ##temp fix
                    #if info_col >= 0 and fields[info_col].startswith('ALT='):
                    #    fields[info_col] = 'o' + fields[info_col]
                    if sample_col >= 0 and fields[sample_col].upper() == args.old_sample.upper():
                        fields[sample_col] = args.new_sample
                        if info_col >= 0:
                            pos = fields[info_col].upper().find('SAMPLE='+args.old_sample.upper())
                            if pos >= 0:
                                new_info = fields[info_col][:pos] + 'SAMPLE='+args.new_sample + fields[info_col][pos+len('SAMPLE='+args.old_sample):]
                                fields[info_col] = new_info
                    new_line = '\t'.join(fields) + '\n'
                    o.write(new_line)



if __name__ == "__main__": sys.exit(main())

