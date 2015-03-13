#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import gzip
import warnings
import my
import sql_columns

class VaxMergePostProcessError(Exception): pass

#-i /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_downloaded_20120816/VCF/parts_00-All.vcf.gz_20120913175843_29gBAW/*.*.00-All.vcf.gz.vex.vcf.vax -o /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_downloaded_20120816/VCF/00-All.vcf.gz.vax.test

name = 'pypeline.vcf2vax'
version = 74
copyright = '©2011-2013 Michael Yourshaw all rights reserved'
description = 'clean up vax files file by replacing "-" with "" and merge files with single set of #-prefixed header lines taken from the first file'
run_time = my.localtime_stamp()

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = description,
        epilog = '{} version {} {}'.format(name,version,copyright))
    parser.add_argument('--input', '-i', nargs='+', required=True,
        help='input vax[.gz] files (first columns must conform to VCF 4.1 format and left-most non-VCF column, if present, must be named Uploaded_variation')
    parser.add_argument('--output', '-o', required=True,
        help='a single processed vax output with metadata and header from the first file and "-" replaced by "" (default:<input>.processed.vax)')
    parser.add_argument('--compress', action='store_true', default=False,
        help='gzip output (default: False)')
    args = parser.parse_args()
    
    input_files = my.unglob(args.input,sort_unique=False)
    my.makedir(os.path.dirname(args.output))
    
    mssql_create_processed_vax_tables_output = '{}.mssql_create_tables.sql'.format(args.output)
    mysql_create_processed_vax_tables_output = '{}.mysql_create_tables.sql'.format(args.output)
    mysqlimport_processed_vax_output = '{}.mysqlimport.sh'.format(args.output)

    columns_in = {}
    processed_vax_firstrow = 0 #0-based index of first data row in output

    skip_header = False
    file_in_count = 0
    line_in_count_total = 0
    record_in_count_total = 0
    line_out_count = 0
    record_out_count = 0

    with gzip.open(args.output, 'w') if args.compress else open(args.output, 'w') as o:
        for input_file in input_files:
            if not my.file_exists(input_file):
               raise VaxMergePostProcessError("vax merge failed: file {} doesn't exist".format(input_file))
            ##DEBUG
            #if file_in_count >= 10: break
            ##DEBUG
            file_in_count += 1
            line_in_count = 0
            record_in_count = 0
            with my.open_gz_or_text(input_file) as input:
                for line in input:
                    line_in_count += 1
                    ##DEBUG
                    #if record_in_count >= 100: break
                    ##DEBUG
                    if line.startswith('#'):
                        #metadata or column header
                        if skip_header:
                            continue
                        elif line.startswith('##'):
                            #metadata
                            o.write(line)
                            line_out_count += 1
                            processed_vax_firstrow += 1
                        elif line.startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                            #column header
                            o.write('## VAX version {} Copyright 2011-2012 Michael Yourshaw all rights reserved\n'.format(version))
                            line_out_count += 1
                            processed_vax_firstrow += 1
                            o.write('## processed VAX output with fields == "-" changed to "" produced at {}\n'.format(run_time))
                            line_out_count += 1
                            processed_vax_firstrow += 1
                            o.write('## VAX --output {}\n'.format(args.output))
                            line_out_count += 1
                            processed_vax_firstrow += 1
                            o.write(line)
                            line_out_count += 1
                            processed_vax_firstrow += 1
                            col_list = line.rstrip('\n').split('\t')
                            has_genotypes= line.startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t') and len(col_list) > 9 and col_list[9] != 'Uploaded_variation'
                            columns_in = {col_list[x]: x for x in range(len(col_list))}
                            #processed_vax_columns_out_spec = my.get_sql_column_spec(columns_in)
                            #processed_vax_columns_out_spec = {col_list[x]: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for x in range(len(col_list))}
                            columns_out = col_list
                        elif not columns_in:
                            raise VaxMergePostProcessError('data encountered before column names at line {}'.format(vax_line_in_count))
                    else:
                        #data record
                        record_in_count += 1
                        data_list = ['' if d == '-' else d for d in line.rstrip('\n').split('\t')]
                        #data_in = my.get_sql_data_dict(columns_in, data_list)
                        #data_in = {k: data_list[columns_in[k]] for k in columns_in.keys()}
                        o.write('\t'.join(data_list)+'\n')
                        line_out_count += 1
                        record_out_count += 1
                        #my.update_sql_column_spec(processed_vax_columns_out_spec, data_in)
                        #accumulate_sql_spec(data_in, processed_vax_columns_out_spec)
                print '{} had {} lines and {} data records'.format(input_file,line_in_count,record_in_count)
                line_in_count_total += line_in_count
                record_in_count_total += record_in_count
            #do not copy metadata and column header after first input file
            skip_header = True
            
    sql_columns.run(args.output, database='vaxdatabase', schema='vaxschema', table=os.path.basename(args.output), clustered_index=['CHROM','POS','REF','ALT','Feature'], indexes=[['Symbol'],['Consequence_rank']])
    #write_microsoft_sql_server_scripts(mssql_create_processed_vax_tables_output, '[dbo].[{}]'.format(os.path.basename(args.output).replace('.','_')), os.path.basename(args.output).replace('.','_'), ['CHROM','POS','REF','ALT','Feature'], [['HGNC'],['Consequence_rank']], col_list, processed_vax_columns_out_spec, processed_vax_firstrow-1)
    #write_mysql_scripts(mysql_create_processed_vax_tables_output, mysqlimport_processed_vax_output, '{}'.format(os.path.basename(args.output).replace('.','_')[:64]), os.path.basename(args.output).replace('.','_')[:16], [['CHROM','POS','REF','ALT','Feature'],['HGNC'],['Consequence_rank']], col_list, processed_vax_columns_out_spec, processed_vax_firstrow)
        
    #finished
    print """vax_post_process done.
    {} input vep files had a total of {} lines and {} records
    vax output file {} had {} lines and {} records""".format(
        file_in_count, line_in_count_total, record_in_count_total,
        args.output, line_out_count, record_out_count)
 
#def accumulate_sql_spec(data_dic, spec_dic):
#    for c in data_dic.keys():
#        this_value = data_dic[c]
#        spec = spec_dic.get(c)
#        if spec and this_value != '':
#            spec['size_min'] = len(this_value) if spec['size_min'] == None else min(len(this_value),spec['size_min'])
#            spec['size_max'] = len(this_value) if spec['size_max'] == None else max(len(this_value),spec['size_max'])
#            if my.is_number(this_value):
#                spec['min'] = float(this_value) if spec['min'] == None else min(float(this_value),spec['min'])
#                spec['max'] = float(this_value) if spec['max'] == None else max(float(this_value),spec['max'])
#            if not spec['type'] == 'varchar':
#                if not my.is_number(this_value):
#                    spec['type'] = 'varchar'
#                elif not spec['type'] == 'float':
#                    if not my.is_int(this_value):
#                        spec['type'] = 'float'
#                    else:
#                        spec['type'] = 'int'

##microsoft sql server
#def write_microsoft_sql_server_scripts(output_file, table_name, index_base_name, clustered_index, indices, columns_out, columns_out_spec, rows_to_delete):
#    with open(output_file, 'w') as ms:
#        ms.write("""SET ANSI_NULLS ON
#GO
#SET QUOTED_IDENTIFIER ON
#GO
#SET ANSI_PADDING ON
#GO
#--substitute your [schema].[table_name] for {0}
#IF OBJECT_ID(N'{0}', N'U') IS NOT NULL 
#DROP TABLE {0};
#GO
#CREATE TABLE {0} (
#""".format(table_name))
#        for c in columns_out:
#            spec = columns_out_spec.get(c)
#            if spec:
#                if c == columns_out[0] and c.startswith('#'):
#                    c = c[1:]
#                if spec['type'] == 'varchar':
#                    if spec['size_min'] == spec['size_max'] and spec['size_max']<=1000:
#                        ms_type = '[CHAR] ({})'.format(spec['size_max'])
#                    else:
#                        ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
#                elif spec['type'] == 'float':
#                    if spec['min'] >= -3.4e+38 and spec['max'] <= 3.4e+38:
#                        ms_type = '[REAL]'
#                    elif spec['min'] >= -1.79e+308 and spec['max'] <= 1.79e+308:
#                        ms_type = '[FLOAT]'
#                    else:
#                        ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
#                elif spec['type'] == 'int':
#                    if spec['min'] >= 0 and spec['max'] <= 1:
#                        ms_type = '[BIT]'
#                    elif spec['min'] >= 0 and spec['max'] <= 255:
#                        ms_type = '[TINYINT]'
#                    elif spec['min'] >= -32768 and spec['max'] <= 32767:
#                        ms_type = '[SMALLINT]'
#                    elif spec['min'] >= -2147483648 and spec['max'] <= 2147483647:
#                        ms_type = '[INT]'
#                    elif spec['min'] >= -9223372036854775808 and spec['max'] <= 9223372036854775807:
#                        ms_type = '[BIGINT]'
#                    else:
#                        ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
#                else:
#                    if c in ('compound_het'):
#                        ms_type = '[BIT]'
#                    else:
#                        ms_type = '[CHAR] (1)'
#                ms_spec = '\t[{}] {} NULL{}\n'.format(c, ms_type, ',' if c != columns_out[-1] else '')
#                ms.write(ms_spec)
#        ms.write(') ON [PRIMARY]\nGO\n\nSET ANSI_PADDING OFF\nGO\n')
#        ms.write("""
#
#--substitute your [schema].[table_name] for {0}
#--substitute your path to unzipped file for {1}
#--delete first {2} rows of comments before importing
#
#DECLARE @bulk_cmd varchar(1000)
#SET @bulk_cmd = 'BULK INSERT {0}
#FROM ''{1}'' 
#WITH (ROWTERMINATOR = '''+CHAR(10)+''',
#FIRSTROW=2)'
#EXEC(@bulk_cmd)
#GO
#""".format(table_name, r'C:\path\to\file\to\be\imported', rows_to_delete))
#
#        #create index scripts
#        #clustered_index is a list of columns
#        if clustered_index:
#            ix_name = 'IX_'+index_base_name+'_'+'_'.join(clustered_index)
#            ix_cols = ',\n\t'.join(['['+c+'] ASC' for c in clustered_index])
#            ms.write("""CREATE CLUSTERED INDEX [{0}] ON {1}
#(
#    {2}
#)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
#GO
#""".format(ix_name, table_name, ix_cols))
#        #indices is a list of indices containing lists of columns
#        for i in indices:
#            ix_name = 'IX_'+index_base_name+'_'+'_'.join(i)
#            ix_cols = ',\n\t'.join(['['+c+'] ASC' for c in i])
#            ms.write("""CREATE NONCLUSTERED INDEX [{0}] ON {1}
#(
#    {2}
#)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
#GO
#""".format(ix_name, table_name, ix_cols))
#
##mysql
#def write_mysql_scripts(create_tables_output_file, mysqlimport_output_file, table_name, index_base_name, indices, columns_out, columns_out_spec, rows_to_delete):
#    with open(create_tables_output_file, 'w') as mysql, open(mysqlimport_output_file, 'w') as mysqlimport:
#        mysql.write("""
#delimiter $$
#DROP TABLE IF EXISTS `{0}`$$
#CREATE TABLE `{0}` (
#""".format(table_name))  #table name length limit of 64 characters
#        for c in columns_out:
#            spec = columns_out_spec.get(c)
#            if spec:
#                if c == columns_out[0] and c.startswith('#'):
#                    c = c[1:]
#                if spec['type'] == 'varchar':
#                    if spec['size_max']<=500:  #row size limit of 65535 bytes, therefore text preferred over varchar, choosing arbitrary limit
#                        my_type = 'varchar({})'.format(spec['size_max'])
#                    else:
#                        my_type = 'text'
#                elif spec['type'] == 'float':
#                    my_type = 'float'
#                elif spec['type'] == 'int':
#                    if spec['min'] >= 0 and spec['max'] <= 1:
#                        my_type = 'tinyint(1)'
#                    elif spec['min'] >= 0 and spec['max'] <= 255:
#                        my_type = 'tinyint unsigned'
#                    elif spec['min'] >= -32768 and spec['max'] <= 32767:
#                        my_type = 'smallint'
#                    elif spec['min'] >= -8388608 and spec['max'] <= 8388607:
#                        my_type = 'mediumint'
#                    elif spec['min'] >= -2147483648 and spec['max'] <= 2147483647:
#                        my_type = 'int'
#                    elif spec['min'] >= -9223372036854775808 and spec['max'] <= 9223372036854775807:
#                        my_type = 'bigint'
#                    elif spec['size_max']<=500:
#                        my_type = 'varchar({})'.format(spec['size_max'])
#                    else:
#                        my_type = 'text'
#                else:
#                    if c in ('compound_het'):
#                        my_type = 'tinyint(1)'
#                    else:
#                        my_type = 'char(1)'
#                my_spec = '\t`{}` {} NULL,\n'.format(c, my_type)
#                mysql.write(my_spec)
#        mysql.write("""
#KEY `IX_vax_CHROM_POS_REF_ALT_Feature` (`CHROM`,`POS`,`REF`,`ALT`,`Feature`),
#KEY `IX_vax_HGNC` (`HGNC`),
#KEY `IX_vax_Consequence_rank` (`Consequence_rank`)""")
#        mysql.write("""
#) ENGINE=MyISAM DEFAULT CHARSET=latin1$$""")
#        #mysqlimport
#        mysqlimport.write("mysqlimport --password=password \\\n--columns=")
#        for c in columns_out:
#            if c == columns_out[0] and c.startswith('#'):
#                    c = c[1:]
#            mysqlimport.write('\'`{}`\'{}'.format(c, ',' if c != columns_out[-1] else ''))
#        mysqlimport.write(' --delete --local --verbose --ignore-lines={} \\\nschema_name \\\n/path/to/{}/file.{}\n'.format(rows_to_delete, index_base_name, index_base_name))
        


if __name__ == "__main__": sys.exit(main())
