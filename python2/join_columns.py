#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import itertools
from math import isinf, isnan
from pprint import pprint as pp 
import my

def AddColumnsError(): pass

#-i /scratch1/tmp/myourshaw/mouse_phenoypes/mus/gut_hormone_table_09_25_2014_mLgr5-GFP.txt -c /scratch1/tmp/myourshaw/mouse_phenoypes/mus/mus_growth_factor.txt -o /scratch1/tmp/myourshaw/mouse_phenoypes/mus/gut_hormone_table_09_25_2014_mLgr5-GFP_growth-factors.txt --input_foreign_key_columns 2 --added_primary_key_columns 2 --columns_to_add 1 3
#-i /scratch1/tmp/myourshaw/mouse_phenoypes/mus/gut_hormone_table_09_25_2014_mLgr5-GFP_growth-factors.txt -c /scratch1/tmp/myourshaw/mouse_phenoypes/mus/mus_growth_factor_receptor.txt -o /scratch1/tmp/myourshaw/mouse_phenoypes/mus/gut_hormone_table_09_25_2014_mLgr5-GFP_growth-factors+receptors.txt --input_foreign_key_columns 2 --added_primary_key_columns 2 --columns_to_add 1 3

def is_int(s):
    try:
        n = int(s)
        if (isinf(n) or isnan(n)):
            return False
        return True
    except ValueError:
        return False


def run(input, columns, output=None,
        input_skip_string=None, columns_skip_string=None,
        ignore_input_lines=None, ignore_columns_lines=None,
        input_header_line=1, columns_header_line=1,
        input_foreign_key_columns=[1], added_primary_key_columns=[1], columns_to_add=[2],):

    if not output:
        output=my.swap_ext('.txt', '.added-columns.txt')
    input_header_line_number = int(input_header_line) if is_int(input_header_line) and int(input_header_line) > 0 else None
    columns_header_line_number = int(columns_header_line) if is_int(columns_header_line) and int(columns_header_line) > 0 else None
    input_header_chars = input_header_line if not is_int(input_header_line) else None    
    columns_header_chars = columns_header_line if not is_int(columns_header_line) else None    
    input_foreign_key_columns = [int(c)-1 for c in input_foreign_key_columns]
    added_primary_key_columns = [int(c)-1 for c in added_primary_key_columns]
    columns_to_add = [int(c)-1 for c in columns_to_add]

    #create dictionary of new columns' data
    columns_dict = {}
    with open(columns,'r') as columns_fh:
        print('Creating database of new columns from {}.'.format(columns))
        line_count = 0
        first_data_line = None
        columns_header_found = False if columns_header_line_number or columns_header_chars else True
        for line in columns_fh:
            line_count += 1
            line = line.rstrip('\r\n')
            if not bool(line.strip()):
                continue
            if not columns_header_found:
                if ((columns_header_line_number and columns_header_line_number == line_count)
                    or (columns_header_chars and columns_header_chars != '#' and line.startswith(columns_header_chars))
                    or (columns_header_chars and columns_header_chars == '#' and line.startswith(columns_header_chars) and line[1] != '#')):
                    header_found = True
                    header_columns = line.split('\t')
                    added_header_column_names = [header_columns[i] for i in columns_to_add]
                    continue
                elif (columns_skip_string and line.startswith(columns_skip_string)) or (ignore_columns_lines and line_count <= ignore_columns_lines):
                    continue
                else: #column data line
                    fields = line.split('\t')
                    key = tuple([fields[i] for i in added_primary_key_columns])
                    value = tuple([fields[i] for i in columns_to_add])
                    columns_dict[key] = value
    
    #add new columns to input data and write to output file
    with open(input,'r') as input_fh, open(output,'w') as output_fh:
        print('Adding new columns to {}.'.format(input))
        line_count = 0
        first_data_line = None
        input_header_found = False if input_header_line_number or input_header_chars else True
        header_found = False
        header_written = False
        for line in input_fh:
            line_count += 1
            line = line.rstrip('\r\n')
            if not bool(line.strip()):
                output_fh.write(line+'\n')
                continue
            if not header_found:
                if ((input_header_line_number and input_header_line_number == line_count)
                      or (input_header_chars and input_header_chars != '#' and line.startswith(input_header_chars))
                      or (input_header_chars and input_header_chars == '#' and line.startswith(input_header_chars) and line[1] != '#')):
                    header_found = True
                    input_column_names = line.split('\t')
                continue
            elif (input_skip_string and line.startswith(input_skip_string)) or (ignore_input_lines and line_count <= ignore_input_lines):
                output_fh.write(line+'\n')
                continue
            else: #input data line
                if first_data_line == None:
                    first_data_line = line_count
                if not header_written:
                    output_column_names = input_column_names + added_header_column_names
                    output_fh.write('{}\n'.format('\t'.join(output_column_names)))
                    header_written = True
                input_fields = line.split('\t')
                foreign_key = tuple([input_fields[i] for i in input_foreign_key_columns if i < len(fields)])
                added_fields = '\t'.join(columns_dict.get(foreign_key,['' for i in columns_to_add]))
                output_line = line + ('\t' if added_fields else '') + added_fields + '\n'
                output_fh.write(output_line)
    print('Done.')



def main():
    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'add pathway annotations to a file',
        epilog = """1.0β1 ©2014 Michael Yourshaw all rights reserved.""")
    parser.add_argument('--input', '-i', type=str, metavar='FILE', required=True,
                        help='Tab-delimited input file')
    parser.add_argument('--columns', '-c', type=str, metavar='FILE', required=True,
                        help='Tab-delimited input file containing columns to be added')
    parser.add_argument('--output', '-o', type=str, metavar='FILE',
                        help='Ouput file (default: <input>.added-columns.txt')
    parser.add_argument('--input_skip_string',
                        help='skip input lines that start with this string')
    parser.add_argument('--columns_skip_string',
                        help='skip columns file lines that start with this string')
    parser.add_argument('--ignore_input_lines', type=int,
                        help='ignore the first n lines of input file')
    parser.add_argument('--ignore_columns_lines', type=int,
                        help='ignore the first n lines of columns file')
    parser.add_argument('--input_header_line', default = '1',
                        help='1-based row number of column names in input file or # if the first row that starts with a single # character (default: 1')
    parser.add_argument('--columns_header_line', default = '1',
                        help='1-based row number of column names in input file or # if the first row that starts with a single # character (default: 1')
    parser.add_argument('--input_foreign_key_columns',  default=[1], nargs='*',
                        help='input file key column numbers (1-based) (default: 1)')
    parser.add_argument('--added_primary_key_columns',  default=[1], nargs='*',
                        help='columns file key column numbers (1-based) (default: 1)')
    parser.add_argument('--columns_to_add',  default=[2], nargs='*',
                        help='columns file column numbers of columns to add (1-based) (default: 2)')
    args = parser.parse_args()    

    run(input=args.input, columns=args.columns, output=args.output,
        input_skip_string=args.input_skip_string, columns_skip_string=args.columns_skip_string,
        ignore_input_lines=args.ignore_input_lines, ignore_columns_lines=args.ignore_columns_lines,
        input_header_line=args.input_header_line, columns_header_line=args.columns_header_line,
        input_foreign_key_columns=args.input_foreign_key_columns, added_primary_key_columns=args.added_primary_key_columns, columns_to_add=args.columns_to_add,)

if __name__ == "__main__": sys.exit(main())