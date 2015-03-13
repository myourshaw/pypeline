#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import gzip
import re
from string import maketrans
import my
import job
import merge_header_files
import vax_post_process

class Vcf2GtxError(Exception): pass

name = 'pypeline.vcf2vex'
version = 2.6
copyright = '©2011-2012 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()

#-i /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_downloaded_20120816/VCF/ByPopulation/YRI-1412-Y.vcf.gz -o /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_downloaded_20120816/VCF/ByPopulation/dbsnp137.gtx --minimal

def run(vcfs, gtx, minimal=False):
        
    config = my.get_config()
        
    my_name = 'vax'
        
    #retun value: a list of job ids that can be used by downstream hold_jid
    job_ids = []
    
    vcfs_glob = vcfs
    vcfs = my.unglob(vcfs)
    if not vcfs:
        raise Vcf2GtxError("no files in {}".format(vcfs_glob))
    
    my.makedir(os.path.dirname(gtx))
    mssql_create_gtx_tables_output = '{}.mssql_create_tables.sql'.format(gtx)
    mysql_create_gtx_tables_output = '{}.mysql_create_tables.sql'.format(gtx)
    mysqlimport_gtx_output = '{}.mysqlimport.sh'.format(gtx)
    
    if minimal:
        gtx_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'GT', 'SAMPLE', 'ALLELE1', 'PHASE', 'ALLELE2']
    else:
        gtx_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT', 'GT', 'SAMPLE', 'ALLELE1', 'PHASE', 'ALLELE2', 'GT_SAMPLE_COUNT', 'GT_ALLELE_COUNT', 'GT_ALLELE1_COUNT', 'GT_ALLELE2_COUNT', 'GT_ALLELE1_SAMPLE_COUNT', 'GT_ALLELE2_SAMPLE_COUNT', 'GT_GENOTYPE_SAMPLE_COUNT', 'GT_ALLELE1_FREQUENCY', 'GT_ALLELE2_FREQUENCY']
    gtx_columns_spec = {c: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for c in gtx_columns}
    gtx_firstrow = 0
    #gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    gt_sep_re = re.compile(r'([/|])')
    info_translation = maketrans(';=','__')

    with open(gtx,'w') as gtx_out:
        gtx_out.write('\t'.join(gtx_columns)+'\n')
        gtx_firstrow += 1
        for vcf in vcfs:
            #create gtx file and sql scripts
            
            print 'processing: '+vcf
            
            columns_in = {}
            sample_column_count = 0
            
            with my.open_gz_or_text(vcf) as vcf_in:
                vcf_line_count = 0
                vcf_record_count = 0
                gtx_record_count = 0
                for line in vcf_in:
                    ##DEBUG
                    #if vcf_record_count > 100: break
                    ##DEBUG
                    vcf_line_count += 1
                    if line.startswith('##'):
                        continue
                    elif line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'):
                        col_list = line.rstrip('\n').split('\t')
                        if len(col_list) > 9 and col_list[9] != 'Uploaded_variation':
                            columns_in = {col_list[x]: x for x in range(len(col_list))}
                            sample_columns = col_list[9:columns_in.get('Uploaded_variation',len(col_list))]
                            if len(sample_columns) != len(set(sample_columns)):
                                raise Vcf2Vex('duplicated sample column names in vcf file {} [{}]'.format(','.join(vcf, sample_columns)))
                    elif not columns_in:
                        raise Vcf2Vex('data encountered before column names at line {}'.format(vcf_line_in_count))
                    else:
                        #data record
                        vcf_record_count += 1
                        data_list = ['' if d == '-' else d for d in line.rstrip('\n').split('\t')]
                        data_in = {k: data_list[columns_in[k]] for k in columns_in.keys()}
                        ref = data_in['REF']
                        alt_list = data_in['ALT'].rstrip(',').split(',')
                        these_alleles = [ref] + alt_list
                                
                        #gtx file expands GT columns to one line per sample
                        if data_in['FORMAT'].split(':')[0].upper() == 'GT':
                            #calculate statistics for all samples at this locus
                            sample_count = 0
                            allele_count = 0
                            allele_counts = {}
                            allele_sample_counts = {}
                            sample_genotype_counts = {}
                            ref_allele_count = 0
                            alt_allele_count = 0
                            het_sample_count = 0
                            hom_ref_sample_count = 0
                            hom_alt_sample_count = 0
                            ref_allele_frequency = 0.0
                            alt_allele_frequency = 0.0
                            for gt in [data_in[s].split(':')[0] for s in sample_columns]:
                                these_allele_indices = gt_sep_re.split(gt)
                                allele1_index = these_allele_indices[0] if len(these_allele_indices) > 0 else ''
                                phase = these_allele_indices[1] if len(these_allele_indices) > 1 else ''
                                allele2_index = these_allele_indices[2] if len(these_allele_indices) > 2 else ''
                                allele1 = '.' if allele1_index == '.' else these_alleles[int(allele1_index)] if my.is_int(allele1_index) else ''
                                allele2 = '.' if allele2_index == '.' else these_alleles[int(allele2_index)] if my.is_int(allele2_index) else ''
                                genotype_key = tuple(sorted((allele1,allele2)))
                                if allele1 != '.' and allele2 != '.':
                                    sample_count+=1
                                    allele_count+=2
                                    allele_counts[allele1] = allele_counts.get(allele1,0) + 1
                                    allele_counts[allele2] = allele_counts.get(allele2,0) + 1
                                    allele_sample_counts[allele1] = allele_sample_counts.get(allele1,0) + 1
                                    if allele1 != allele2:
                                        allele_sample_counts[allele2] = allele_sample_counts.get(allele2,0) + 1
                                    sample_genotype_counts[genotype_key] = sample_genotype_counts.get(genotype_key,0) + 1
                                    
                            #write lines for each sample at this locus
                            for this_sample in sample_columns:
                                gtx_data_out = {c: data_in.get(c,'') for c in gtx_columns}
                                gtx_data_out['SAMPLE'] = this_sample
                                gtx_data_out['INFO'] = 'SAMPLE={};'.format(this_sample.translate(info_translation)) + data_in['INFO']
                                gt = data_in[this_sample].split(':')[0]
                                if minimal:
                                    gtx_data_out['GT'] = gt
                                else:
                                    gtx_data_out['GT'] = data_in[this_sample]
                                these_allele_indices = gt_sep_re.split(gt)
                                allele1_index = these_allele_indices[0] if len(these_allele_indices) > 0 else ''
                                phase = these_allele_indices[1] if len(these_allele_indices) > 1 else ''
                                allele2_index = these_allele_indices[2] if len(these_allele_indices) > 2 else ''
                                allele1 = '.' if allele1_index == '.' else these_alleles[int(allele1_index)] if my.is_int(allele1_index) else ''
                                gtx_data_out['ALLELE1'] = allele1
                                allele2 = '.' if allele2_index == '.' else these_alleles[int(allele2_index)] if my.is_int(allele2_index) else ''
                                gtx_data_out['ALLELE2'] = allele2
                                gtx_data_out['PHASE'] = phase
                                if allele1 != '.' and allele2 != '.':
                                    genotype_key = tuple(sorted((allele1,allele2)))
                                    gtx_data_out['GT_SAMPLE_COUNT'] = str(sample_count)
                                    gtx_data_out['GT_ALLELE_COUNT'] = str(allele_count)
                                    gtx_data_out['GT_ALLELE1_COUNT'] = str(allele_counts[allele1])
                                    gtx_data_out['GT_ALLELE2_COUNT'] = str(allele_counts[allele2])
                                    gtx_data_out['GT_ALLELE1_SAMPLE_COUNT'] = str(allele_sample_counts[allele1])
                                    gtx_data_out['GT_ALLELE2_SAMPLE_COUNT'] = str(allele_sample_counts[allele2])
                                    gtx_data_out['GT_GENOTYPE_SAMPLE_COUNT'] = str(sample_genotype_counts[genotype_key])
                                    if allele_count != 0:
                                        allele1_frequency = float(allele_counts[allele1])/float(allele_count)
                                        allele2_frequency = float(allele_counts[allele2])/float(allele_count)
                                    gtx_data_out['GT_ALLELE1_FREQUENCY'] = str(allele1_frequency)
                                    gtx_data_out['GT_ALLELE2_FREQUENCY'] = str(allele2_frequency)
                                gtx_out.write('\t'.join([gtx_data_out.get(c,'') for c in gtx_columns])+'\n')
                                gtx_record_count += 1
                                accumulate_sql_spec(gtx_data_out, gtx_columns_spec)
    
    write_microsoft_sql_server_scripts(mssql_create_gtx_tables_output, '[dbo].[{}]'.format(os.path.basename(gtx)), 'gtx', [['CHROM','POS','REF','ALT'],['SAMPLE'],['ALLELE1'],['ALLELE2']], gtx_columns, gtx_columns_spec, gtx_firstrow-1)
    
    write_mysql_scripts(mysql_create_gtx_tables_output, mysqlimport_gtx_output, '{}'.format(os.path.basename(gtx)[:64]), 'gtx', [['CHROM','POS','REF','ALT'],['SAMPLE'],['ALLELE1'],['ALLELE2']], gtx_columns, gtx_columns_spec, gtx_firstrow)

    print 'done'
    
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
def write_microsoft_sql_server_scripts(output_file, table_name, index_base_name, indices, columns_out, columns_out_spec, rows_to_delete):
    with open(output_file, 'w') as ms:
        
        #create table script
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
        
        #import data script
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
""".format(table_name, r'C:\path\to\file\to\be\imported', rows_to_delete))

        #create index scripts
        #indices is a tuple of indices containing tuples of columns
        for i in indices:
            ix_name = 'IX_'+index_base_name+'_'+'_'.join(i)
            ix_cols = ',\n\t'.join(['['+c+'] ASC' for c in i])
            ms.write("""CREATE NONCLUSTERED INDEX [{0}] ON {1}
(
    {2}
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO
""".format(ix_name, table_name, ix_cols))


#mysql
def write_mysql_scripts(create_tables_output_file, mysqlimport_output_file, table_name, index_base_name, indices, columns_out, columns_out_spec, rows_to_delete):
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
        for i in indices:
            ix_name = 'IX_'+index_base_name+'_'+'_'.join(i)
            ix_cols = '`,`'.join(i)
            mysql.write("""
KEY `{0}` ({1}),""".format(ix_name, ix_cols))
        mysql.write("""
) ENGINE=MyISAM DEFAULT CHARSET=latin1$$""")
        #mysqlimport
        mysqlimport.write("mysqlimport --password=password \\\n--columns=")
        for c in columns_out:
            if c == columns_out[0] and c.startswith('#'):
                    c = c[1:]
            mysqlimport.write('\'`{}`\'{}'.format(c, ',' if c != columns_out[-1] else ''))
        mysqlimport.write(' --delete --local --verbose --ignore-lines={} \\\nschema_name \\\n/path/to/{}/file.{}\n'.format(rows_to_delete, index_base_name, index_base_name))

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'merge vcf files and expand sample columns',
        epilog = 'pypeline.vcf2vex version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    #input VCF parameter
    parser.add_argument('--vcfs', '-v', '-i', nargs='+', required=True,
        help='input vcf file(s); output will be <vcf>.vax')
    parser.add_argument('--gtx', '-g', '-o', required=True,
        help='output gtx file')
    parser.add_argument('--minimal', action='store_true', default=False,
        help='minimal gtx columns (CHROM,POS,ID,REF,ALT,GT,SAMPLE')
    args = parser.parse_args()

    #dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    job_ids = run(vcfs=args.vcfs, gtx=args.gtx, minimal=args.minimal)


if __name__ == "__main__": sys.exit(main())
