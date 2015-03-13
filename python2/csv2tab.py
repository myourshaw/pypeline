#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import csv
import argparse

class Csv2TabError(Exception): pass

#-i /Users/myourshaw/apps/vax/downloads/hpa/normal_tissue.csv

def run(input, output=None, delimiter=',', quotechar='"', output_delimiter='\t', output_quotechar=None):
    if not input:
        raise Csv2TabError('argument --input/-i is required')
    if not output:
        if output_delimiter == '\t':
            output = input + '.tab'
        elif output_delimiter == ',':
            output = input + '.csv'
        else:
            output = input + '.converted'
    with open(input, 'rb') as csvfile, open(output, 'wb') as tabfile:
        reader = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
        writer = csv.writer(tabfile, delimiter=output_delimiter, quotechar=output_quotechar) if output_quotechar \
        else csv.writer(tabfile, delimiter=output_delimiter, quoting=csv.QUOTE_NONE)
        for row in reader:
            writer.writerow(row)

def main():   

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'convert csv file to tab-delimited',
        epilog = 'pypeline.csv2tab version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required = True,
        help='input  file')
    parser.add_argument('--output', '-o',
        help='output  file (default: <input>.tab)')
    parser.add_argument('--delimiter', '-d', default=',',
        help='input field delimiter (default: ,')
    parser.add_argument('--quotechar', '-q', default='"',
        help='input quote character (default: ")')
    parser.add_argument('--output_delimiter', default = '\t',
        help='output field delimiter (default: ,')
    parser.add_argument('--output_quotechar', default=None,
        help='input quote character (default: csv.QUOTE_NONE)')
    args = parser.parse_args()
    
    run(input=args.input, output=args.output, delimiter=args.delimiter, quotechar=args.quotechar, output_delimiter=args.output_delimiter, output_quotechar=args.output_quotechar)

if __name__ == "__main__": sys.exit(main())
    
    
