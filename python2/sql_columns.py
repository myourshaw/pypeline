#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import time
from warnings import warn
import copy
import my

class SqlColumnsError(Exception): pass

#--header_line 28 -i /scratch1/tmp/myourshaw/resources/rgd/20140328/rgd_homo_terms_bp.txt
#-db v -s vax -ci GeneSymbol -i /scratch1/tmp/myourshaw/resources/clinvar/gene_condition_source_id.txt
#--header_line 1 -i /scratch1/vax/75/hpa/hpa_subcellular_location.txt
#-i /scratch1/tmp/myourshaw/resources/intervals/Nimblegen/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_capture.bed -s nimblegen -pk chrom chromStart chromEnd --column_names chrom chromStart chromEnd name --ignore_lines 1
#-i /scratch1/tmp/myourshaw/rnaseq_20140227/mus/mCD24.txt --header_line 1 -x gene

def run(input, database='database', schema='schema', table=None,
        primary_key=None, clustered_index=None, unique_indexes=None, indexes=None,
        skip_string=None, ignore_lines=None, header_line='#', column_names=None, no_column_names=False,
        no_microsoft=False, no_mysql=False, engine='MyISAM'):
        
    header_line_number = int(header_line) if my.is_int(header_line) and int(header_line) > 0 else None
    header_chars = header_line if not my.is_int(header_line) else None
    
    files = my.unglob(input)
    
    if len(files) > 1 and (table or primary_key or clustered_index or indexes):
        warn("""table, primary_key, clustered_index, and indexes will be applied to all files
Is this what you want?
Consider using a for loop with a different set of these parameters for each file file. For example,
for p in \
"-i file1.txt --table1 --primary_key k1,k2 --input file1" \
"-i file2.txt --table2 --primary_key k3,k4 --input file2" \
; do
python sql_columns.py ${p};
done""")
    
    for file in files:
        print('Analyzing {}.'.format(file))
        microsoft_sql_file = file+'.microsoft.sql'
        mysql_file = file+'.mysql'
        microsoft_done_file = my.done_file(microsoft_sql_file)
        mysql_done_file = my.done_file(mysql_file)
        
        with my.open_gz_or_text(file) as input, open(microsoft_sql_file, 'w') as microsoft, open(file+'.mysql', 'w') as mysql:
            column_widths = {}
            column_names_list = None
            line_count = 0
            first_data_line = None
            header_found = False if header_line_number or header_chars else True
            if column_names:
                header_found = True
                column_names_list = column_names
                these_column_names = {i:column_names_list[i] for i in range(len(column_names_list))}
                sql_column_spec = my.get_sql_column_spec(column_names_list)
            else:
                these_column_names = {}
            for line in input:
                line_count += 1
                line = line.rstrip('\r\n')
                if not bool(line.strip()):
                    continue
                if not header_found:
                    if no_column_names:
                        header_found = True
                        column_names_list = ['col'+str(i+1) for i in range(len(line.split('\t')))]
                        these_column_names = {i:column_names_list[i] for i in range(len(column_names_list))}
                        sql_column_spec = my.get_sql_column_spec(column_names_list)
                    elif (header_line_number and header_line_number == line_count) or (header_chars and header_chars != '#' and line.startswith(header_chars)) or (header_chars and header_chars == '#' and line.startswith(header_chars) and line[1] != '#'):
                        header_found = True
                        column_names_list = line.split('\t')
                        these_column_names = {i:column_names_list[i] for i in range(len(column_names_list))}
                        sql_column_spec = my.get_sql_column_spec(column_names_list)
                    continue
                elif (skip_string and line.startswith(skip_string)) or (ignore_lines and line_count <= ignore_lines):
                    continue
                else: #data line
                    if first_data_line == None:
                        first_data_line = line_count
                    fields = line.split('\t')
                    ii = [i for i in range(len(fields))]
                    column_widths = {i: max(len(fields[i]),column_widths.setdefault(i,0)) for i in ii}
                    if len(fields) != len(these_column_names):
                        raise SqlColumnsError(
                            ('Number of fields in line {} is not the same as the number of expected columns.\n' +
                            'Columns:\t{}\n' +
                            'Fields:\t{}').format(
                                line_count, '\t'.join(these_column_names.values()), '\t'.join(fields)))
                    sql_data_dict = my.get_sql_data_dict(these_column_names, fields)
                    my.update_sql_column_spec(sql_column_spec, sql_data_dict)

            if not (column_names_list and sql_column_spec):
                raise SqlColumnsError('Could not identify columns. Try using --header_line, --column_names, or --no_column_names.')
            microsoft_sql_server_scripts = my.get_microsoft_sql_server_scripts(
                data_file=file,
                database=database,
                schema=schema,
                table_name=table if table else my.r_strip(os.path.basename(file),'.txt'),
                primary_key=primary_key,
                clustered_index=clustered_index,
                unique_indexes=unique_indexes,
                indexes=indexes,
                columns_out=column_names_list,
                columns_out_spec=sql_column_spec,
                rows_to_delete=first_data_line-2
                )
            if not no_microsoft:
                microsoft.write('\n\n'.join(microsoft_sql_server_scripts))

            mysql_scripts = my.get_mysql_scripts(
                data_file=file,
                database=database,
                schema=schema,
                table_name=table if table else os.path.splitext(os.path.basename(file))[0],
                primary_key=primary_key,
                clustered_index=clustered_index,
                unique_indexes=unique_indexes,
                indexes=indexes,
                columns_out=column_names_list,
                columns_out_spec=sql_column_spec,
                rows_to_delete=first_data_line-1,
                engine=engine,
                )
            if not no_mysql:
                mysql.write('\n\n'.join(mysql_scripts))
        with open(microsoft_done_file, 'w'), open(mysql_done_file, 'w'):
            pass


def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'generate CREATE TABLE and import statements',
        epilog = 'pypeline.column_widths version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required = True, nargs='+',
        help='input  file(s)[.gz]')
    parser.add_argument('--database', '-db', default = 'database',
        help='database (default: database)')
    parser.add_argument('--schema', '-s', default = 'schema',
        help='schema (default: schema)')
    parser.add_argument('--table', '-t',
        help='table name (default: file name minus .txt)')
    parser.add_argument('--primary_key', '-pk', nargs='*',
        help='primary key columns, space-separated')
    parser.add_argument('--clustered_index', '-ci', nargs='*',
        help='clustered_index columns space-separated')
    parser.add_argument('--unique_indexes', '-ui', nargs='*',
        help='space-separated list of unique indexes in the form of comma-separated lists of column names, e.g., name,employee_number ssn')
    parser.add_argument('--indexes', '-x', nargs='*',
        help='space-separated list of indexes in the form of comma-separated lists of column names, e.g., name,employee_number ssn')
    parser.add_argument('--skip_string',
        help='skip lines that start with this string')
    parser.add_argument('--ignore_lines', type=int,
        help='ignore the first n lines')
    parser.add_argument('--header_line', default = '#',
        help='1-based row number of column names or # if the first row that starts with a single # character')
    parser.add_argument('--column_names', nargs='*',
        help='list of column names (when file does not contain column name header)')
    parser.add_argument('--no_column_names', action='store_true', default=False,
        help='file does not contain column name header, so use "col1, col2, ..."')
    parser.add_argument('--no_microsoft', action='store_true', default=False,
        help='do no produce microsoft file')
    parser.add_argument('--no_mysql', action='store_true', default=False,
        help='do no produce mysql file')
    parser.add_argument('--engine', default = 'MyISAM',
        help='mysql database engine')
    args = parser.parse_args()

    run(input=args.input, table=args.table, database=args.database, schema=args.schema,
        primary_key=args.primary_key, clustered_index=args.clustered_index, unique_indexes=args.unique_indexes, indexes=args.indexes,
        skip_string=args.skip_string, ignore_lines=args.ignore_lines, header_line=args.header_line, column_names=args.column_names, no_column_names=args.no_column_names,
        no_microsoft=args.no_microsoft, no_mysql=args.no_mysql, engine=args.engine)

if __name__ == "__main__": sys.exit(main())
