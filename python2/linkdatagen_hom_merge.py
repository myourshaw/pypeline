#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from ConfigParser import SafeConfigParser
import re
import my

class LinkdatagenHomMergeError(Exception): pass

#-i /scratch1/tmp/myourshaw/vf_20140828/plink/KMC*_plink/plink.hom -o /scratch1/tmp/myourshaw/vf_20140828/plink/vf_homozygous_intervals.txt

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'merge plink.hom files produced by linkdatagen',
        epilog = 'pypeline.linkdatagen_hom_merge version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', nargs='+',
        help='input plink.hom files')
    parser.add_argument('--output', '-o', required=True,
        help='merged ouput file')
    args = parser.parse_args()

    files = my.unglob(args.input)
    in_col_headings = 'FID  IID      PHE  CHR        SNP1        SNP2         POS1         POS2         KB     NSNP  DENSITY     PHOM     PHET'.split()
    out_col_headings = ['file', 'sample'] + in_col_headings
    with open(args.output, 'w') as o:
        o.write('\t'.join(out_col_headings)+'\n')
        for file in files:
            sample = os.path.basename(os.path.dirname(file))[:-6]
            with open(file) as f:
                for line in f:
                    fields = line.rstrip('\n').split()
                    if not fields or fields == in_col_headings:
                        continue
                    if fields[3] == '23':
                        fields[3] = 'X'
                    elif fields[3] == '24':
                        fields[3] = 'Y'
                    out_data = [file,sample] + fields
                    o.write('\t'.join(out_data)+'\n')
                    

if __name__ == "__main__": sys.exit(main())
