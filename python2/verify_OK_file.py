#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from ConfigParser import SafeConfigParser
import my

class VerifyOKFileError(Exception): pass

#-v /scratch0/tmp/myourshaw/mms/bams/sample_bams/MMS69.sample.recalibrated.bam.validate

def main():
   
    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'check that a file exists and contains only the two characters "OK" without an end of line marker',
        epilog = 'pypeline.verify_OK_file version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-v', '-i', required=True,
        help='input file')
    args = parser.parse_args()

    if not my.check_files(files_contain=[(args.input, 'OK')]):
        raise ValidateSamFileError('{} had validation errrors\n{}'.format(args.input, val_msg))

    return 0

if __name__ == "__main__": sys.exit(main())
