#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from warnings import warn
import re
import my

class ValidateRecalCsvError(Exception): pass

#-i /scratch0/tmp/myourshaw/gmd/bams/sample_bams/zila_temp/GMD108.sample.recal.csv
#-i /scratch0/tmp/myourshaw/gmd/bams/sample_bams/zila_temp/GMD108.sample.recal2.csv

def write_error(err_count, max_errors, err_msg, output_file):
    err_count+=1
    if err_count > max_errors:
        raise ValidateRecalCsvError('recal.csv file has > {} errors; see {}'.format(max_errors, output_file.name))
    else:
        output_file.write(err_msg+'\n')
        return err_count

def run(input, output=None, max_errors=1):
    if not output:
        output = input+'.validate'

    ReadGroup,QualityScore,Cycle,Dinuc,nObservations,nMismatches,Qempirical = range(7)

    recal_csv_re = re.compile(r'(?P<ReadGroup>[^,]+),(?P<QualityScore>[-.0-9eE]+),(?P<Cycle>[-0-9]+),(?P<Dinuc>[ACGTN]{2}),(?P<nObservations>[0-9]+),(?P<nMismatches>[0-9]+),(?P<Qempirical>[-.0-9eE]+)', re.I)

    err_count = 0
    line_count = 0
    header_seen = False
    data_seen = False
    eof_seen = False
    
    output_file = open(output, 'w')
    try:
        with open(input) as v:
            for line in v:
                line_count += 1
                line = line.rstrip('\n')
                if line.startswith('#'):
                    continue
                elif line == 'ReadGroup,QualityScore,Cycle,Dinuc,nObservations,nMismatches,Qempirical':
                    header_seen = True
                    continue
                elif not header_seen:
                    err_msg = 'line {} missing header {}'.format(line_count, line)
                    err_count = write_error(err_count, max_errors, err_msg, output_file)
                elif line == 'EOF':
                    eof_seen = True
                else:
                    if eof_seen:
                        err_msg = 'line {} data after EOF {}'.format(line_count, line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    else:
                        m = recal_csv_re.match(line)
                        if not m:
                            err_msg = 'line {} invalid data {}'.format(line_count, line)
                            err_count = write_error(err_count, max_errors, err_msg, output_file)
                        else:
                            data_seen = True
        if header_seen and data_seen and eof_seen:
            output_file.write('OK')
        else:
            err_msg = 'line {} header and/or data and/or EOF missing {}'.format(line_count, line)
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
        description = 'validate recal.csv file',
        epilog = 'pypeline.validate_recal_csv version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='input recal.csv file to be validated')
    parser.add_argument('--output', '-o',
        help='output file(default: <input recal.csv>.validate; if input is valid file output contains "OK")')
    parser.add_argument('--max_errors', type=int, default=1,
        help='maximum number of errors to report (default: 1)')
    args = parser.parse_args()
    
    status = run(args.input, args.output, args.max_errors)
    
    return status


if __name__ == "__main__": sys.exit(main())
