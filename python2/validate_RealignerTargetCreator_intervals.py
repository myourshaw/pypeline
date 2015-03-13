#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from warnings import warn
import re
import my

class ValidateRecalCsvError(Exception): pass

#-i /scratch0/tmp/myourshaw/gmd/bams/sample_bams/zila_temp/GMD108.sample.realignertargetcreator.intervals

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

    interval_re = re.compile(r'(?P<chrom>[^:]+):(?P<startPos>[0-9]+)(?:-(?P<endPos>[0-9]+))?', re.I)

    err_count = 0
    line_count = 0
    that_chrom = None
    that_pos = 0
    data_seen = False
    
    output_file = open(output, 'w')
    try:
        with open(input) as v:
            for line in v:
                line_count += 1
                line = line.rstrip('\n')
                if line.startswith('#') or not line.strip():
                    continue
                else:
                    m = interval_re.match(line)
                    if not m:
                        err_msg = 'line {} invalid data {}'.format(line_count, line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    else:
                        foo = m.groupdict(None)
                        if m.group('chrom') == that_chrom:
                            if int(m.group('startPos')) < that_pos:
                                err_msg = 'line {} startPos {} < {} {}'.format(line_count, m.group('startPos'), that_pos, line)
                                err_count = write_error(err_count, max_errors, err_msg, output_file)
                            else:
                                that_pos = int(m.group('endPos')) if m.lastgroup == 'endPos' else int(m.group('startPos'))
                        else:
                            that_chrom = m.group('chrom')
                            that_pos = int(m.group('endPos')) if m.lastgroup == 'endPos' else int(m.group('startPos'))
                            data_seen = True
        if data_seen:
            output_file.write('OK')
        else:
            err_msg = 'line {} data missing {}'.format(line_count, line)
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
        description = 'validate realignertargetcreator.intervals file',
        epilog = 'pypeline.validate_RealignerTargetCreator_intervals version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='input realignertargetcreator.intervals file to be validated')
    parser.add_argument('--output', '-o',
        help='output file(default: <input realignertargetcreator.intervals>.validate; if input is valid file output contains "OK")')
    parser.add_argument('--max_errors', type=int, default=1,
        help='maximum number of errors to report (default: 1)')
    args = parser.parse_args()
    
    status = run(args.input, args.output, args.max_errors)
    
    return status


if __name__ == "__main__": sys.exit(main())
