#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from ConfigParser import SafeConfigParser
import my

class ValidateSamFileError(Exception): pass

#-v /scratch0/tmp/myourshaw/mms/bams/sample_bams/MMS69.sample.recalibrated.bam.validate

def main():
   
    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'check output of picard ValidateSamFile',
        epilog = 'pypeline.verify_sam_file_validation version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--validate_file', '-v', required=True,
        help='picard ValidateSamFile ouput file')
    args = parser.parse_args()

    if not my.check_files(files_contain=[(args.validate_file, 'No errors found\n')]):
        raise ValidateSamFileError('{} had validation errrors\n{}'.format(args.validate_file,))

    return 0

if __name__ == "__main__": sys.exit(main())
