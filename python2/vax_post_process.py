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

class VaxPostProcessError(Exception): pass

#-i /scratch1/tmp/myourshaw/hlee/vcfs/hlee55gmd28bipolar20ocdtic18.analysis_ready.vcf.vax
#-i /home/myourshaw/lab/pypeline/vax_test/gmd_test.vcf.vax
#-i /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/chr22.ln.vcf.vax.gz
#--flax_consequence_threshhold 100 -i /home/myourshaw/lab/pypeline/vax_test/intergenic_1kg_chr22.vcf.vax
#--flax_consequence_threshhold -1 -i /home/myourshaw/lab/pypeline/vax_test/gmd_test.vcf.vax
#-i /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_downloaded_20120816/VCF/00-All.vcf.gz.vax --processed_vax /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_downloaded_20120816/VCF/00-All.vcf.gz.processed.vax.test

def main():

    version = 2.6
    run_time = my.localtime_stamp()

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'clean up vax file file by replacing "-" with ""',
        epilog = 'pypeline.vax2flax version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='input vcf[.gz] or vax[.gz] file (first columns must conform to VCF 4.1 format and left-most non-VCF column, if present, must be named Uploaded_variation')
    parser.add_argument('--processed_vax',
        help='processed vax output with "-" replaced by "" (default:<input>.processed.vax)')
    parser.add_argument('--count_uniques', action='store_true', default=False,
        help='count unique CHROM, POS, REF, ALT (default:False)')
   #parser.add_argument('--vax_consequence_threshhold', type=int, default=None,
    #    help='only output revised vax records where Consequence_rank <= vax_consequence_threshhold; negative number for all consequences (default: None)')
    args = parser.parse_args()
    
    if not args.processed_vax:
        args.processed_vax = args.input+'.processed.vax'
    #if not args.vax_consequence_threshhold or not my.is_number(args.vax_consequence_threshhold) or args.vax_consequence_threshhold < 0:
    #    args.vax_consequence_threshhold = None
        
    mssql_create_processed_vax_tables_output = '{}.mssql_create_tables.sql'.format(args.processed_vax)
    mysql_create_processed_vax_tables_output = '{}.mysql_create_tables.sql'.format(args.processed_vax)
    mysqlimport_processed_vax_output = '{}.mysqlimport.sh'.format(args.processed_vax)

    additional_columns = ['SAMPLE', 'ALLELE1', 'PHASE', 'ALLELE2', 'ZYGOSITY', 'compound_het', 'GT', 'GT_INFO']

    columns_in = {}
    sample_column_count = 0
    processed_vax_firstrow = 0 #0-based index of first data row

    #gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    gt_sep_re = re.compile(r'([/|])')

    with my.open_gz_or_text(args.input) as i, open(args.processed_vax, 'w') as o:
        vax_line_in_count = 0
        vax_record_in_count = 0
        processed_vax_record_out_count = 0
        unique_variants_in = set()
        processed_vax_unique_variants_out = set()
        for line in i:
            ##DEBUG
            #if flax_record_out_count > 100: break
            ##DEBUG
            vax_line_in_count += 1
            if line.startswith('##'):
                o.write(line)
                processed_vax_firstrow += 1
            elif line.startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                o.write('## VAX version {} Copyright 2011-2012 Michael Yourshaw all rights reserved\n'.format(version))
                processed_vax_firstrow += 1
                o.write('## processed VAX output with fields == "-" changed to "" produced at {}\n'.format(run_time))
                processed_vax_firstrow += 1
                o.write('## VAX --input {}\n'.format(args.input))
                processed_vax_firstrow += 1
                o.write('## VAX --processed_vax {}\n'.format(args.processed_vax))
                processed_vax_firstrow += 1
                #o.write('## VAX --vax_consequence_threshhold {}\n'.format(args.vax_consequence_threshhold))
                #processed_vax_firstrow += 1
                o.write(line)
                processed_vax_firstrow += 1
                col_list = line.rstrip('\n').split('\t')
                has_genotypes= line.startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t') and len(col_list) > 9 and col_list[9] != 'Uploaded_variation'
                columns_in = {col_list[x]: x for x in range(len(col_list))}
                processed_vax_columns_out_spec = {col_list[x]: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for x in range(len(col_list))}
                columns_out = col_list
            elif not columns_in:
                raise VaxPostProcessError('data encountered before column names at line {}'.format(vax_line_in_count))
            else:
                #data record
                vax_record_in_count += 1
                data_list = ['' if d == '-' else d for d in line.rstrip('\n').split('\t')]
                data_in = {k: data_list[columns_in[k]] for k in columns_in.keys()}
                o.write('\t'.join(data_list)+'\n')
                if args.count_uniques:
                    unique_variants_in.add((data_in['#CHROM'],data_in['POS'],data_in['REF'],data_in['ALT']))
                    processed_vax_unique_variants_out.add((data_in['#CHROM'],data_in['POS'],data_in['REF'],data_in['ALT']))
                processed_vax_record_out_count += 1
                accumulate_sql_spec(data_in, processed_vax_columns_out_spec)

    write_microsoft_sql_server_scripts(mssql_create_processed_vax_tables_output, '[dbo].[{}]'.format(os.path.basename(args.processed_vax)), 'vax', (('CHROM','POS','REF','ALT','Feature'),('HGNC'),('Consequence_rank')), col_list, processed_vax_columns_out_spec, processed_vax_firstrow-1)
    write_mysql_scripts(mysql_create_processed_vax_tables_output, mysqlimport_processed_vax_output, '{}'.format(os.path.basename(args.processed_vax)[:64]), 'processed_vax', (('CHROM','POS','REF','ALT','Feature'),('HGNC'),('Consequence_rank')), col_list, processed_vax_columns_out_spec, processed_vax_firstrow)
        
    #finished
    print """vax_post_process done.
    input vax file {} had {} lines, {} records, {} unique CHROM, POS, REF, ALT
    output processed_vax file {} had {} records, {} unique CHROM, POS, REF, ALT""".format(
        args.input, vax_line_in_count, vax_record_in_count, len(unique_variants_in),
        args.processed_vax, processed_vax_record_out_count, len(processed_vax_unique_variants_out))
 
def accumulate_sql_spec(data_dic, spec_dic):
    for c in data_dic.keys():
        this_value = data_dic[c]
        spec = spec_dic.get(c)
        if spec and this_value != '':
            spec['size_min'] = len(this_value) if spec['size_min'] == None else min(len(this_value),spec['size_min'])
            spec['size_max'] = len(this_value) if spec['size_max'] == None else max(len(this_value),spec['size_max'])
            if my.is_number(this_value):
                spec['min'] = float(this_value) if spec['min'] == None else min(float(this_value),spec['min'])
                spec['max'] = float(this_value) if spec['max'] == None else max(float(this_value),spec['max'])
            if not spec['type'] == 'varchar':
                if not my.is_number(this_value):
                    spec['type'] = 'varchar'
                elif not spec['type'] == 'float':
                    if not my.is_int(this_value):
                        spec['type'] = 'float'
                    else:
                        spec['type'] = 'int'

#microsoft sql server
def write_microsoft_sql_server_scripts(output_file, table_name, index_base_name, sample_indices, columns_out, columns_out_spec, rows_to_delete):
    with open(output_file, 'w') as ms:
        ms.write("""SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
SET ANSI_PADDING ON
GO
--substitute your [schema].[table_name] for {0}
IF OBJECT_ID(N'{0}', N'U') IS NOT NULL 
DROP TABLE {0};
GO
CREATE TABLE {0} (
""".format(table_name))
        for c in columns_out:
            spec = columns_out_spec.get(c)
            if spec:
                if c == columns_out[0] and c.startswith('#'):
                    c = c[1:]
                if spec['type'] == 'varchar':
                    if spec['size_min'] == spec['size_max'] and spec['size_max']<=1000:
                        ms_type = '[CHAR] ({})'.format(spec['size_max'])
                    else:
                        ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
                elif spec['type'] == 'float':
                    if spec['min'] >= -3.4e+38 and spec['max'] <= 3.4e+38:
                        ms_type = '[REAL]'
                    elif spec['min'] >= -1.79e+308 and spec['max'] <= 1.79e+308:
                        ms_type = '[FLOAT]'
                    else:
                        ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
                elif spec['type'] == 'int':
                    if spec['min'] >= 0 and spec['max'] <= 1:
                        ms_type = '[BIT]'
                    elif spec['min'] >= 0 and spec['max'] <= 255:
                        ms_type = '[TINYINT]'
                    elif spec['min'] >= -32768 and spec['max'] <= 32767:
                        ms_type = '[SMALLINT]'
                    elif spec['min'] >= -2147483648 and spec['max'] <= 2147483647:
                        ms_type = '[INT]'
                    elif spec['min'] >= -9223372036854775808 and spec['max'] <= 9223372036854775807:
                        ms_type = '[BIGINT]'
                    else:
                        ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
                else:
                    if c in ('compound_het'):
                        ms_type = '[BIT]'
                    else:
                        ms_type = '[CHAR] (1)'
                ms_spec = '\t[{}] {} NULL{}\n'.format(c, ms_type, ',' if c != columns_out[-1] else '')
                ms.write(ms_spec)
        ms.write(') ON [PRIMARY]\nGO\n\nSET ANSI_PADDING OFF\nGO\n')
        ms.write("""

--substitute your [schema].[table_name] for {0}
--substitute your path to unzipped file for {1}
--delete first {2} rows of comments before importing

DECLARE @bulk_cmd varchar(1000)
SET @bulk_cmd = 'BULK INSERT {0}
FROM ''{1}'' 
WITH (ROWTERMINATOR = '''+CHAR(10)+''',
FIRSTROW=2)'
EXEC(@bulk_cmd)
GO

CREATE NONCLUSTERED INDEX [IX_{3}_CHROM_POS_REF_ALT_Feature] ON {0}
(
	[CHROM] ASC,
	[POS] ASC,
	[REF] ASC,
	[ALT] ASC,
	[Feature] ASC
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO

CREATE NONCLUSTERED INDEX [IX_{3}_HGNC] ON {0}
(
	[HGNC] ASC
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO

CREATE NONCLUSTERED INDEX [IX_{3}_Consequence_rank] ON {0}
(
	[Consequence_rank] ASC
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO
""".format(table_name, r'C:\path\to\file\to\be\imported', rows_to_delete, index_base_name))

        if sample_indices:
            ms.write("""
CREATE NONCLUSTERED INDEX [IX_{3}_SAMPLE_CHROM_POS_REF_ALT_Feature] ON {0}
(
	[SAMPLE] ASC,
	[CHROM] ASC,
	[POS] ASC,
	[REF] ASC,
	[ALT] ASC,
	[Feature] ASC
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO

CREATE NONCLUSTERED INDEX [IX_{3}_ZYGOSITY] ON {0}
(
	[ZYGOSITY] ASC
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO

CREATE NONCLUSTERED INDEX [IX_{3}_SAMPLE_Feature_POS] ON {0}
(
	[SAMPLE] ASC,
	[Feature] ASC,
	[POS] ASC
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO

CREATE NONCLUSTERED INDEX [IX_{3}_compound_het] ON {0}
(
	[compound_het] ASC
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO
""".format(table_name, r'C:\path\to\file\to\be\imported', rows_to_delete, index_base_name))

#mysql
def write_mysql_scripts(create_tables_output_file, mysqlimport_output_file, table_name, index_base_name, sample_indices, columns_out, columns_out_spec, rows_to_delete):
    with open(create_tables_output_file, 'w') as mysql, open(mysqlimport_output_file, 'w') as mysqlimport:
        mysql.write("""
delimiter $$
DROP TABLE IF EXISTS `{0}`$$
CREATE TABLE `{0}` (
""".format(table_name))  #table name length limit of 64 characters
        for c in columns_out:
            spec = columns_out_spec.get(c)
            if spec:
                if c == columns_out[0] and c.startswith('#'):
                    c = c[1:]
                if spec['type'] == 'varchar':
                    if spec['size_max']<=500:  #row size limit of 65535 bytes, therefore text preferred over varchar, choosing arbitrary limit
                        my_type = 'varchar({})'.format(spec['size_max'])
                    else:
                        my_type = 'text'
                elif spec['type'] == 'float':
                    my_type = 'float'
                elif spec['type'] == 'int':
                    if spec['min'] >= 0 and spec['max'] <= 1:
                        my_type = 'tinyint(1)'
                    elif spec['min'] >= 0 and spec['max'] <= 255:
                        my_type = 'tinyint unsigned'
                    elif spec['min'] >= -32768 and spec['max'] <= 32767:
                        my_type = 'smallint'
                    elif spec['min'] >= -8388608 and spec['max'] <= 8388607:
                        my_type = 'mediumint'
                    elif spec['min'] >= -2147483648 and spec['max'] <= 2147483647:
                        my_type = 'int'
                    elif spec['min'] >= -9223372036854775808 and spec['max'] <= 9223372036854775807:
                        my_type = 'bigint'
                    elif spec['size_max']<=500:
                        my_type = 'varchar({})'.format(spec['size_max'])
                    else:
                        my_type = 'text'
                else:
                    if c in ('compound_het'):
                        my_type = 'tinyint(1)'
                    else:
                        my_type = 'char(1)'
                my_spec = '\t`{}` {} NULL,\n'.format(c, my_type)
                mysql.write(my_spec)
        mysql.write("""
KEY `IX_flax` (`CHROM`,`POS`,`REF`,`ALT`,`Feature`),""")
        if sample_indices:
            mysql.write("""
KEY `IX_flax_1` (`SAMPLE`,`CHROM`,`POS`,`REF`,`ALT`,`Feature`),
KEY `IX_flax_2` (`SAMPLE`,`Feature`,`POS`),
KEY `IX_flax_3` (`ZYGOSITY`),
KEY `IX_flax_4` (`compound_het`),""")
        mysql.write("""
KEY `IX_flax_5` (`HGNC`),
KEY `IX_flax_6` (`Consequence_rank`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1$$""")
        #mysqlimport
        mysqlimport.write("mysqlimport --password=password \\\n--columns=")
        for c in columns_out:
            if c == columns_out[0] and c.startswith('#'):
                    c = c[1:]
            mysqlimport.write('\'`{}`\'{}'.format(c, ',' if c != columns_out[-1] else ''))
        mysqlimport.write(' --delete --local --verbose --ignore-lines={} \\\nschema_name \\\n/path/to/{}/file.{}\n'.format(rows_to_delete, index_base_name, index_base_name))
        


if __name__ == "__main__": sys.exit(main())
