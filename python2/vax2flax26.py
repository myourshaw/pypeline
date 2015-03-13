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

class Vax2Flax(Exception): pass

#-i /scratch1/tmp/myourshaw/hlee/vcfs/hlee55gmd28bipolar20ocdtic18.analysis_ready.vcf.vax
#-i /home/myourshaw/lab/pypeline/vax_test/gmd_test.vcf.vax
#-i /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/chr22.ln.vcf.vax.gz
#--flax_consequence_threshhold 100 -i /home/myourshaw/lab/pypeline/vax_test/intergenic_1kg_chr22.vcf.vax
#--flax_consequence_threshhold -1 -i /home/myourshaw/lab/pypeline/vax_test/gmd_test.vcf.vax

def main():

    version = 2.6
    run_time = my.localtime_stamp()

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'flatten VAX output file by samples, and filter by Consequence_rank',
        epilog = 'pypeline.vax2flax version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='input vcf[.gz] or vax[.gz] file (first columns must conform to VCF 4.1 format and left-most non-VCF column, if present, must be named Uploaded_variation')
    parser.add_argument('--output', '-o',
        help='output file (default:<input>.flax.gz)')
    parser.add_argument('--revised_vax',
        help='revised vax output with "-" replaced by "" and optionally filtered by consequence (default:<input>.revised_vax.gz)')
    parser.add_argument('--no_revised_vax', action='store_true', default=False,
        help='do not create revised vax output with "-" replaced by "" and optionally filtered by consequence (default: False)')
    parser.add_argument('--vax_consequence_threshhold', type=int, default=None,
        help='only output revised vax records where Consequence_rank <= vax_consequence_threshhold; negative number for all consequences (default: None)')
    parser.add_argument('--flax_consequence_threshhold', type=int, default=None, #14, #Incomplete terminal codon variant in Ensembl 68
        help='only output flax records where Consequence_rank <= flax_consequence_threshhold; negative number for all consequences (default: 9 = non_synonymous_codon)')
    parser.add_argument('--flax_suppress_uncalled_genotypes', action='store_true', default=False,
        help='do not output uncalled (e.g., ./.) genotypes flax records (default: False)')
    parser.add_argument('--flax_suppress_homozygous_reference_genotypes', action='store_true', default=False,
        help='do not output homozygous_reference flax records (e.g., 0/0) genotypes (default: False)')
    parser.add_argument('--flax_include_all_genotype_columns', action='store_true', default=False,
        help='include the genotype columns for all samples in flax records (default: False)')
    parser.add_argument('--revised_vax_uncompressed', action='store_true', default=False,
        help='do not gzip revised vax output (default: False)')
    parser.add_argument('--flax_uncompressed', action='store_true', default=False,
        help='do not gzip flax output (default: False)')
    args = parser.parse_args()
    
    if not args.output:
        args.output = '{}.flax{}'.format(args.input,'.gz' if not args.flax_uncompressed else '')
    if not args.vax_consequence_threshhold or not my.is_number(args.vax_consequence_threshhold) or args.vax_consequence_threshhold < 0:
        args.vax_consequence_threshhold = None
    if not args.flax_consequence_threshhold or not my.is_number(args.flax_consequence_threshhold) or args.flax_consequence_threshhold < 0:
        args.flax_consequence_threshhold = None
    if args.flax_uncompressed and args.output.endswith('.gz'):
        args.output = args.output[:-3]
    if not args.flax_uncompressed and not args.output.endswith('.gz'):
        args.output = args.output+'.gz'
    if not args.revised_vax:
        args.revised_vax = '{}.revised_vax{}'.format(args.input,'.gz' if not args.revised_vax_uncompressed else '')
    if args.revised_vax_uncompressed and args.revised_vax.endswith('.gz'):
        args.revised_vax = args.output[:-3]
    if not args.revised_vax_uncompressed and not args.revised_vax.endswith('.gz'):
        args.revised_vax = args.output+'.gz'
        
    mssql_create_flax_tables_output = '{}.flax.mssql_create_tables.sql'.format(args.input)
    mssql_create_revised_vax_tables_output = '{}.revised_vax.mssql_create_tables.sql'.format(args.input)
    mysql_create_flax_tables_output = '{}.flax.mysql_create_tables.sql'.format(args.input)
    mysql_create_revised_vax_tables_output = '{}.revised_vax.mysql_create_tables.sql'.format(args.input)
    mysqlimport_flax_output = '{}.flax.mysqlimport.sh'.format(args.input)
    mysqlimport_revised_vax_output = '{}.revised_vax.mysqlimport.sh'.format(args.input)

    additional_columns = ['SAMPLE', 'ALLELE1', 'PHASE', 'ALLELE2', 'ZYGOSITY', 'compound_het', 'GT', 'GT_INFO']

    columns_in = {}
    sample_column_count = 0
    flax_firstrow = 0 #0-based index of first data row
    revised_vax_firstrow = 0 #0-based index of first data row

    #gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    gt_sep_re = re.compile(r'([/|])')

    with my.open_gz_or_text(args.input) as i, gzip.open(args.output, 'w') if not args.flax_uncompressed else open(args.output, 'w') as o:
        vax = None if args.no_revised_vax else gzip.open(args.revised_vax, 'w') if not args.revised_vax_uncompressed else open(args.revised_vax, 'w')
        vax_line_in_count = 0
        vax_record_in_count = 0
        flax_record_out_count = 0
        revised_vax_record_out_count = 0
        unique_variants_in = set()
        flax_unique_variants_out = set()
        revised_vax_unique_variants_out = set()
        for line in i:
            ##DEBUG
            #if flax_record_out_count > 100: break
            ##DEBUG
            vax_line_in_count += 1
            if line.startswith('##'):
                o.write(line)
                if vax:
                    vax.write(line)
                flax_firstrow += 1
                revised_vax_firstrow += 1
            elif line.startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                if vax:
                    vax.write('## VAX version {} Copyright 2011-2012 Michael Yourshaw all rights reserved\n'.format(version))
                    revised_vax_firstrow += 1
                    vax.write('## revised VAX output with fields == "-" changed to "" produced at {}\n'.format(run_time))
                    revised_vax_firstrow += 1
                    vax.write('## VAX --input {}\n'.format(args.input))
                    revised_vax_firstrow += 1
                    vax.write('## VAX --revised_vax {}\n'.format(args.revised_vax))
                    revised_vax_firstrow += 1
                    vax.write('## VAX --vax_consequence_threshhold {}\n'.format(args.vax_consequence_threshhold))
                    revised_vax_firstrow += 1
                    vax.write(line)
                    revised_vax_firstrow += 1
                col_list = line.rstrip('\n').split('\t')
                has_genotypes= line.startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t') and len(col_list) > 9 and col_list[9] != 'Uploaded_variation'
                columns_in = {col_list[x]: x for x in range(len(col_list))}
                revised_vax_columns_out_spec = {col_list[x]: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for x in range(len(col_list))}
                columns_out = col_list
                o.write('## FLAX version {} Copyright 2011-2012 Michael Yourshaw all rights reserved\n'.format(version))
                flax_firstrow += 1
                o.write('## FLAX output produced at {}\n'.format(run_time))
                flax_firstrow += 1
                o.write('## FLAX --input {}\n'.format(args.input))
                flax_firstrow += 1
                o.write('## FLAX --output {}\n'.format(args.output))
                flax_firstrow += 1
                o.write('## FLAX --flax_consequence_threshhold {}\n'.format(args.flax_consequence_threshhold))
                flax_firstrow += 1
                o.write('## FLAX --flax_suppress_uncalled_genotypes {}\n'.format(args.flax_suppress_uncalled_genotypes))
                flax_firstrow += 1
                o.write('## FLAX --flax_suppress_homozygous_reference_genotypes {}\n'.format(args.flax_suppress_homozygous_reference_genotypes))
                flax_firstrow += 1
                if has_genotypes:
                    columns_out = columns_out + additional_columns
                    sample_columns = col_list[9:columns_in.get('Uploaded_variation',len(col_list))]
                    if not args.flax_include_all_genotype_columns:
                        columns_out = [c for c in columns_out if c not in sample_columns]
                    o.write('## SAMPLE       : ID of sample from VCF column header (FLAX)\n')
                    flax_firstrow += 1
                    o.write('## ALLELE1      : first allele (FLAX)\n')
                    flax_firstrow += 1
                    o.write('## PHASE        : / = genotype unphased | = genotype phased (FLAX)\n')
                    flax_firstrow += 1
                    o.write('## ALLELE2      : second allele, if present (Multiple alleles >diploid not used) (FLAX)\n')
                    flax_firstrow += 1
                    o.write('## ZYGOSITY     : 0 = homozygous or hemizygous REF, 1 = heterozygous, 2 = homozygous or hemizygous ALT (FLAX)\n')
                    flax_firstrow += 1
                    o.write('## compound_het : 1 if variant involved in compound heterozygous group (always blank in this version) (FLAX)\n')
                    flax_firstrow += 1
                    o.write('## GT           : genotype from sample field (FLAX)\n')
                    flax_firstrow += 1
                    o.write('## GT_INFO      : complete contents of sample field (FLAX)\n')
                    flax_firstrow += 1
                flax_columns_out_spec = {columns_out[x]: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for x in range(len(columns_out))}
                if len(columns_out) != len(set(columns_out)):
                    raise Vax2Flax('duplicated column names in {}'.format(','.join(columns_out)))
                o.write('\t'.join(columns_out)+'\n')
                flax_firstrow += 1
            elif not columns_in:
                raise Vax2Flax('data encountered before column names at line {}'.format(vax_line_in_count))
            else:
                #data record
                vax_record_in_count += 1
                data_list = ['' if d == '-' else d for d in line.rstrip('\n').split('\t')]
                data_in = {k: data_list[columns_in[k]] for k in columns_in.keys()}
                unique_variants_in.add((data_in['#CHROM'],data_in['POS'],data_in['REF'],data_in['ALT']))
                alt_list = data_in['ALT'].rstrip(',').split(',')
                for alt in alt_list:
                    data_in['ALT'] = alt
                    if vax and not (args.vax_consequence_threshhold and columns_in.get('Consequence_rank') and my.is_int(data_in['Consequence_rank']) and int(data_in['Consequence_rank']) > args.vax_consequence_threshhold):
                        vax.write('\t'.join(data_list)+'\n')
                        revised_vax_unique_variants_out.add((data_in['#CHROM'],data_in['POS'],data_in['REF'],data_in['ALT']))
                        revised_vax_record_out_count += 1
                        accumulate_sql_spec(data_in, revised_vax_columns_out_spec)
                if args.flax_consequence_threshhold and columns_in.get('Consequence_rank') and my.is_int(data_in['Consequence_rank']) and int(data_in['Consequence_rank']) > args.flax_consequence_threshhold:
                    continue
                data_in['compound_het'] = ''
                ref = data_in['REF']
                these_alleles = [ref] + alt_list
                for alt in alt_list:
                    data_in['ALT'] = alt
                    flax_unique_variants_out.add((data_in['#CHROM'],data_in['POS'],data_in['REF'],data_in['ALT']))
                    if has_genotypes:
                        for sample in sample_columns:
                            data_in['SAMPLE'] = sample
                            data_in['GT_INFO'] = data_in[sample]
                            if data_in['FORMAT'].split(':')[0].upper() == 'GT':
                                data_in['GT'] = data_in[sample].split(':')[0]
                                these_allele_indices = gt_sep_re.split(data_in['GT'])
                                allele1_index = these_allele_indices[0] if len(these_allele_indices) > 0 else ''
                                phase = these_allele_indices[1] if len(these_allele_indices) > 1 else ''
                                allele2_index = these_allele_indices[2] if len(these_allele_indices) > 2 else ''
                                allele1 = '.' if allele1_index == '.' else these_alleles[int(allele1_index)] if my.is_int(allele1_index) else ''
                                data_in['ALLELE1'] = allele1
                                allele2 = '.' if allele2_index == '.' else these_alleles[int(allele2_index)] if my.is_int(allele2_index) else ''
                                data_in['ALLELE2'] = allele2
                                data_in['PHASE'] = phase
                                data_in['ZYGOSITY'] = '0' if allele1 == ref and (allele2 == ref or allele2 == '') else '2' if allele1 == alt and (allele2 == alt or allele2 == '') else '' if allele1 == '.' or allele1 == '' or allele2 == '.' or allele2 == '' else '1' if allele1 != allele2 and allele1_index < allele2_index else ''
                                if args.flax_suppress_homozygous_reference_genotypes and data_in['ZYGOSITY'] == '0':
                                    continue
                                if args.flax_suppress_uncalled_genotypes and data_in['ZYGOSITY'] == '':
                                    continue
                            data_out = [data_in.get(c,'') for c in columns_out]
                            o.write('\t'.join(data_out)+'\n')
                            flax_record_out_count += 1
                            accumulate_sql_spec(data_in, flax_columns_out_spec)
                    else:
                        #no genotypes
                        accumulate_sql_spec(data_in, flax_columns_out_spec)
    if vax:
        vax.close()

    write_microsoft_sql_server_scripts(mssql_create_flax_tables_output, '[dbo].[{}]'.format(os.path.basename(args.output)), 'flax', True, columns_out, flax_columns_out_spec, flax_firstrow-1)
    write_microsoft_sql_server_scripts(mssql_create_revised_vax_tables_output, '[dbo].[{}]'.format(os.path.basename(args.revised_vax)), 'vax', False, col_list, revised_vax_columns_out_spec, revised_vax_firstrow-1)
    
    write_mysql_scripts(mysql_create_flax_tables_output, mysqlimport_flax_output, '{}'.format(os.path.basename(args.output)[:64]), 'flax', True, columns_out, flax_columns_out_spec, flax_firstrow)
    write_mysql_scripts(mysql_create_revised_vax_tables_output, mysqlimport_revised_vax_output, '{}'.format(os.path.basename(args.output)[:64]), 'revised_vax', False, col_list, revised_vax_columns_out_spec, revised_vax_firstrow)
        
    #finished
    print """vax2flax done.
    input vax file {} had {} lines, {} records, {} unique CHROM, POS, REF, ALT
    output flax file {} had {} records, {} unique CHROM, POS, REF, ALT
    output revised_vax file {} had {} records, {} unique CHROM, POS, REF, ALT""".format(
        args.input, vax_line_in_count, vax_record_in_count, len(unique_variants_in),
        args.output, flax_record_out_count, len(flax_unique_variants_out),
        args.revised_vax, revised_vax_record_out_count, len(revised_vax_unique_variants_out))
 
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
