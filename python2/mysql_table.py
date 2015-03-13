#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import time
import warnings
import copy
import my

class SqlColumnsError(Exception): pass

#--header_line 28 -x 'OBJECT_SYMBOL' 'TERM_NAME' -i /scratch1/tmp/myourshaw/resources/rgd/20140328/rgd_*_terms.txt

def run(input, table=None, primary_key=None, indexes=None, skip=None, header_line='#', column_names=None, no_column_names=False):
    
    header_line_number = int(header_line) if my.is_int(header_line) and int(header_line) > 0 else None
    header_chars = header_line if not my.is_int(header_line) else None
    
    if not table:
        table = my.r_strip(os.path.basename(input), '.txt')
    
    output = os.path.join(os.path.dirname(input),table+'.mysql')
    with my.open_gz_or_text(input) as input_file, open(output, 'w') as mysql:
        column_widths = {}
        column_names = {}
        column_names_list = None
        line_count = 0
        first_data_line = None
        header_found = False if header_line_number or header_chars else True
        if column_names:
            header_found = True
            column_names_list = column_names
            column_names = {i:column_names_list[i] for i in range(len(column_names_list))}
            sql_column_spec = my.get_sql_column_spec(column_names_list)
        for line in input_file:
            line_count += 1
            line = line.rstrip('\n')
            if not bool(line.strip()):
                continue
            if not header_found:
                if no_column_names:
                    header_found = True
                    column_names_list = ['col'+str(i+1) for i in range(len(line.split('\t')))]
                    column_names = {i:column_names_list[i] for i in range(len(column_names_list))}
                    sql_column_spec = my.get_sql_column_spec(column_names_list)
                elif (header_line_number and header_line_number == line_count) or (header_chars and header_chars != '#' and line.startswith(header_chars)) or (header_chars and header_chars == '#' and line.startswith(header_chars) and line[1] != '#'):
                    header_found = True
                    column_names_list = line.split('\t')
                    column_names = {i:column_names_list[i] for i in range(len(column_names_list))}
                    sql_column_spec = my.get_sql_column_spec(column_names_list)
                continue
            elif skip and line.startswith(skip):
                continue
            else: #data line
                if first_data_line == None:
                    first_data_line = line_count
                fields = line.split('\t')
                indices = [i for i in range(len(fields))]
                column_widths = {i: max(len(fields[i]),column_widths.setdefault(i,0)) for i in indices}
                sql_data_dict = my.get_sql_data_dict(column_names, fields)
                my.update_sql_column_spec(sql_column_spec, sql_data_dict)
                pass
        #output.write('col_number\tcol_name\tcol_width\n')
        #output.write(''.join(['{}\t{}\t{}\n'.format(col+1, column_names[col], column_widths[col]) for col in sorted(column_widths)]))

        if not (column_names_list and sql_column_spec):
            raise SqlColumnsError('Could not identify columns. Try using --header_line if there is a header, --column_names or --no_column_names, if there is not.')

        mysql_scripts = my.get_mysql_scripts(
            table_name=table,
            index_base_name=table,
            primary_key=primary_key,
            indexes=indexes,
            columns_out=column_names_list,
            columns_out_spec=sql_column_spec,
            rows_to_delete=first_data_line-1
            )
        mysql.write(mysql_scripts[0])

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'calculate maximum column widths for tab-delimited file',
        epilog = 'pypeline.column_widths version 1.0β1 ©2011-2013 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required = True,
        help='input  file(s)[.gz]')
    parser.add_argument('--table', '-t',
        help='input  file[.gz]')
    parser.add_argument('--primary_key', '-pk',
        help='primary key in the form of comma-separated list of column names, e.g., col1,col3')
    parser.add_argument('--indexes', '-x', nargs='*',
        help='space-separated list of indexes in the form of comma-separated lists of column names, e.g., name,employee_number ssn')
    parser.add_argument('--skip',
        help='skip lines that start with character(s)')
    parser.add_argument('--header_line', default = '#',
        help='1-based row number of column names or # if the first row that starts with a single # character or string that header row starts with')
    parser.add_argument('--column_names', nargs='*',
        help='list of column names (when file does not contain column name header)')
    parser.add_argument('--no_column_names', action='store_true', default=False,
        help='file does not contain column name header, so use "col1, col2, ..."')
    args = parser.parse_args()
    
    if not my.file_exists(args.input):
        raise SqlColumnsError('The file '+args.input+' does not exist.')
    if args.table and len(args.table) > 64:
        raise SqlColumnsError('Table name '+args.table+' is too long (max 64).')
    run(input=args.input, table=args.table, primary_key=args.primary_key, indexes=args.indexes,
        skip=args.skip, header_line=args.header_line,
        column_names=args.column_names, no_column_names=args.no_column_names,)

if __name__ == "__main__": sys.exit(main())
