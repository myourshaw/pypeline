#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from warnings import warn
import re
import my

class ValidateBrlmmError(Exception): pass

#-i /home/myourshaw/lab/pypeline/vax23/gmd_test.vcf.vax
#-i /home/myourshaw/lab/pypeline/vax23/gmd_test.vcf.vax.gt.txt

def write_error(err_count, max_errors, err_msg, output_file):
    err_count+=1
    if err_count > max_errors:
        raise ValidateBrlmmError('brlmm file has > {} errors; see {}'.format(max_errors, output_file.name))
    else:
        output_file.write(err_msg+'\n')
        return err_count

def run(input, output=None, max_errors=1):
    if not output:
        output = input+'.validate'

    err_count = 0
    line_count = 0
    
    output_file = open(output, 'w')
    try:
        with open(input) as v:
            for line in v:
                line_count += 1
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 2 or fields[1] not in ('0','1','2'):
                    err_msg = 'line {} "{}" is wrong format'.format(line_count, line.rstrip('\n'))
                    err_count = write_error(err_count, max_errors, err_msg, output_file)
        if err_count == 0:
            output_file.write('OK')
        else:
            raise ValidateBrlmmError('brlmm file has {} errors; see {}'.format(err_count, output_file.name))
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
        description = 'validate linkdatagen brmm file',
        epilog = 'pypeline.validate_brmm_file version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='linkdatagen brmm file to be validated')
    parser.add_argument('--output', '-o',
        help='output file(default: <input file>.validate; if input is valid file output contains "OK")')
    parser.add_argument('--max_errors', type=int, default=1,
        help='maximum number of errors to report (default: 1)')
    args = parser.parse_args()
    
    status = run(args.input, args.output, args.max_errors)
    
    return status


if __name__ == "__main__": sys.exit(main())
