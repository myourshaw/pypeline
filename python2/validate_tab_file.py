#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from warnings import warn
import re
import my

class ValidateTabError(Exception): pass

#-i /home/myourshaw/lab/pypeline/vax23/gmd_test.vcf.vax
#-i /home/myourshaw/lab/pypeline/vax23/gmd_test.vcf.vax.gt.txt

def write_error(err_count, max_errors, err_msg, output_file):
    err_count+=1
    if err_count > max_errors:
        raise ValidateTabError('recal.csv file has > {} errors; see {}'.format(max_errors, output_file.name))
    else:
        output_file.write(err_msg+'\n')
        return err_count

def run(input, output=None, max_errors=1):
    if not output:
        output = input+'.validate'

    err_count = 0
    line_count = 0
    header_seen = False
    column_count = 0
    data_seen = False
    
    output_file = open(output, 'w')
    try:
        with open(input) as v:
            for line in v:
                line_count += 1
                line = line.rstrip('\n')
                if line.startswith('##'):
                    continue
                elif not header_seen and line.startswith('#'):
                    header_seen = True
                    column_count = len(line.split('\t'))
                    continue
                elif not header_seen:
                    err_msg = 'line {} missing column name header {}'.format(line_count, line)
                    err_count = write_error(err_count, max_errors, err_msg, output_file)
                else:
                    if len(line.split('\t')) != column_count:
                        err_msg = 'line {} {} data fields instead of {} {}'.format(line_count, len(line.split('\t')), column_count, line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    else:
                        data_seen = True
        if header_seen and data_seen:
            output_file.write('OK')
        else:
            err_msg = 'line {} column header and/or data missing {}'.format(line_count, line)
            err_count = write_error(err_count, max_errors, err_msg, output_file)
    except Exception as e:
        err_msg = 'error {} validating input file {}.format(e, input)'
        err_count = write_error(err_count, max_errors, err_msg, output_file)
    else:
        return 0 if err_count == 0 else 100
    finally:
        output_file.close()

def main():
   
    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'validate tab delimited file',
        epilog = 'pypeline.validate_tab_file version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='tab-delimited file to be validated, with optional "##..." header lines and mandatory "#..." column names line and no empty lines')
    parser.add_argument('--output', '-o',
        help='output file(default: <input file>.validate; if input is valid file output contains "OK")')
    parser.add_argument('--max_errors', type=int, default=1,
        help='maximum number of errors to report (default: 1)')
    args = parser.parse_args()
    
    status = run(args.input, args.output, args.max_errors)
    
    return status


if __name__ == "__main__": sys.exit(main())
