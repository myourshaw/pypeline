#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import pipes
import shutil
import my

#-q /scratch0/tmp/myourshaw/hlee/merge_test/s_1_1_110[123]_head_qseq.txt.gz -o /scratch0/tmp/myourshaw/hlee/merge_test/s_1_1_1101-3_qseq.txt.gz

class MergeValidateQseqsError(Exception): pass

def run(qseqs, output, validate=None, overwrite=False):
    
    my.makedir(os.path.dirname(output))
    if not output.endswith('.gz'):
        output = output + '.gz'

    if not validate:
        validate = output + '.validate'

    #do nothing if valid output exists unless overwrite
    if overwrite or not my.check_files([output], [(validate,'OK')]):
        
        qseqs = my.unglob(qseqs)
        
        #validate input qseqs
        for q in qseqs:
            if not my.qseq_all_reads_valid(q):
                raise MergeValidateQseqsError('invalid qseq input file ' + q)

        #merge input qseq files and compress the output
        t = pipes.Template()
        t.append('gzip','--')
        with t.open(output, 'w') as o:
            for q in qseqs:
                shutil.copyfileobj(my.open_gz_or_text(q), o)

        #validate output qseq
        with open(validate,'w') as v:
            output_valid = my.qseq_all_reads_valid(output)
            v.write('OK' if output_valid else 'ERROR')
            if not output_valid:
                raise MergeValidateQseqsError('invalid qseq output file ' + output)

    return 0
    

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'validate, merge, compress, and validate qseq files',
        epilog = 'pypeline.merge_validate_qseqs version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--qseqs', '-q', nargs='+', required=True,
        help='list of compressed or uncompressed qseq files to be merged')
    parser.add_argument('--output', '-o', required=True,
        help='output file, which will be compressed and ".gz" appened if nexessary')
    parser.add_argument('--validate', '-v',
        help='validation results file; contains "OK" if valid, "ERROR" if invalid (default: <output>.gz.validate)')
    parser.add_argument('--overwrite', action='store_true', default=False,
        help='overwrite existing valid output')
    args = parser.parse_args()
    
    results = run(qseqs=args.qseqs, output=args.output, validate=args.validate, overwrite=args.overwrite)
    

if __name__ == "__main__": sys.exit(main())
