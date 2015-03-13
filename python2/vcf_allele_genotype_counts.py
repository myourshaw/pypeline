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

class VcfAlleleGenotypeCountsError(Exception): pass

name = 'pypeline.vcf_allele_genotype_counts'
version = 2.6
copyright = '©2011-2012 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()

#-i /scratch0/tmp/myourshaw/niehs/20120903/niehs.95samples.polymorphic.filtered.indels.vcf -o /scratch0/tmp/myourshaw/niehs/20120903/niehs.95samples.polymorphic.filtered.indels.vcf
#-i /scratch0/tmp/myourshaw/1000genomes/phase1_integrated_release_version3_20120430/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -o /scratch0/tmp/myourshaw/1000genomes/phase1_integrated_release_version3_20120430/allele_genotype_counts/test
#-i /data/storage-1-01/archive/myourshaw/scratch0/1000genomes/phase1_analysis_results_integrated_call_sets_20121012/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz -o /data/storage-1-01/archive/myourshaw/scratch0/1000genomes/phase1_analysis_results_integrated_call_sets_20121012/test
#-i /data/storage-1-01/archive/myourshaw/scratch0/1000genomes/phase1_analysis_results_integrated_call_sets_20121012/ALL.chrMT.phase1_samtools_si.20101123.snps.low_coverage.genotypes.vcf.gz -o /data/storage-1-01/archive/myourshaw/scratch0/1000genomes/phase1_analysis_results_integrated_call_sets_20121012/ALL.chrMT.phase1_samtools_si.20101123.snps.low_coverage.genotypes.vcf.gz.test

def run(vcfs, output, uncompressed=False):
        
    config = my.get_config()
        
    my_name = 'vax'
        
    #retun value: a list of job ids that can be used by downstream hold_jid
    job_ids = []
    
    vcfs_glob = vcfs
    vcfs = my.unglob(vcfs)
    if not vcfs:
        raise VcfAlleleGenotypeCountsError("no files in {}".format(vcfs_glob))
    
    my.makedir(os.path.dirname(output))
    allele_file = output + '.allele_counts.txt' + ('' if uncompressed else '.gz')
    genotype_file = output+ '.genotype_counts.txt' + ('' if uncompressed else '.gz')
    mssql_create_allele_tables_output = '{}.mssql_create_tables.sql'.format(allele_file)
    mysql_create_allele_tables_output = '{}.mysql_create_tables.sql'.format(allele_file)
    mssql_create_genotype_tables_output = '{}.mssql_create_tables.sql'.format(genotype_file)
    mysql_create_genotype_tables_output = '{}.mysql_create_tables.sql'.format(genotype_file)
    
    allele_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'allele', 'allele_count']
    genotype_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'allele1', 'allele2', 'sample_count']
    #allele_columns_spec = {c: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for c in allele_columns}
    allele_columns_spec = my.get_sql_column_spec(allele_columns)
    #genotype_columns_spec = {c: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for c in genotype_columns}
    genotype_columns_spec = my.get_sql_column_spec(genotype_columns)
    allele_firstrow = 0
    genotype_firstrow = 0
    #gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    gt_sep_re = re.compile(r'([/|])')
    gt_re = re.compile(r'[/|]')
    info_translation = maketrans(';=','__')

    #with open(allele_file,'w') as allele_out, open(genotype_file,'w') as genotype_out:
    with open(allele_file,'w') if uncompressed else gzip.open(allele_file, 'w') as allele_out, open(genotype_file,'w') if uncompressed else gzip.open(genotype_file, 'w') as genotype_out:
        allele_out.write('\t'.join(allele_columns)+'\n')
        allele_firstrow += 1
        genotype_out.write('\t'.join(genotype_columns)+'\n')
        genotype_firstrow += 1
        for vcf in vcfs:
            
            print 'processing: '+vcf
            
            columns_in = {}
            sample_column_count = 0
            
            with my.open_gz_or_text(vcf) as vcf_in:
                vcf_line_count = 0
                vcf_record_count = 0
                allele_record_count = 0
                genotype_record_count = 0
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
                                raise VcfAlleleGenotypeCountsError('duplicated sample column names in vcf file {} [{}]'.format(','.join(vcf, sample_columns)))
                    elif line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                        raise VcfAlleleGenotypeCountsError('no genotypes in {}'.format(vcf))
                    elif not columns_in:
                        raise VcfAlleleGenotypeCountsError('data encountered before column names in {} at line {}'.format(vcf, vcf_line_count))
                    else:
                        #data record
                        vcf_record_count += 1
                        data_list = ['' if d == '-' else d for d in line.rstrip('\n').split('\t')]
                        data_in = {k: data_list[columns_in[k]] for k in columns_in.keys()}
                        chrom,pos,id,ref,alt,qual,filter,info,format = data_in['#CHROM'],data_in['POS'],data_in['ID'],data_in['REF'],data_in['ALT'],data_in['QUAL'],data_in['FILTER'],data_in['INFO'],data_in['FORMAT']
                        genotype_data_out = {}
                        genotype_data_out['#CHROM'],genotype_data_out['POS'],genotype_data_out['ID'],genotype_data_out['REF'],genotype_data_out['ALT'],genotype_data_out['QUAL'],genotype_data_out['FILTER'],genotype_data_out['INFO'],genotype_data_out['FORMAT'] = chrom,pos,id,ref,alt,qual,filter,info,format
                        allele_data_out = {}
                        allele_data_out['#CHROM'],allele_data_out['POS'],allele_data_out['ID'],allele_data_out['REF'],allele_data_out['ALT'],allele_data_out['QUAL'],allele_data_out['FILTER'],allele_data_out['INFO'],allele_data_out['FORMAT'] = chrom,pos,id,ref,alt,qual,filter,info,format
                        alt_list = data_in['ALT'].rstrip(',').split(',')
                        these_alleles = [ref] + alt_list
                        #gtx file expands GT columns to one line per sample
                        if data_in['FORMAT'].split(':')[0].upper() == 'GT':
                            allele_counts = {}
                            gt_counts = {}
                            these_gts = [tuple(sorted(gt_sep_re.split(data_in[s].split(':')[0])[::2])) for s in sample_columns]
                            for gt in these_gts:
                                gt_counts[gt] = gt_counts.get(gt,0)+1
                            for gt in sorted(gt_counts.keys()):
                                allele1_index = gt[0] if len(gt) > 0 else ''
                                allele2_index = gt[1] if len(gt) > 1 else ''
                                allele1 = '.' if allele1_index == '.' else these_alleles[int(allele1_index)] if my.is_int(allele1_index) else ''
                                allele2 = '.' if allele2_index == '.' else these_alleles[int(allele2_index)] if my.is_int(allele2_index) else ''
                                genotype_data_out['allele1'],genotype_data_out['allele2'] = allele1,allele2
                                genotype_data_out['sample_count'] = str(gt_counts[gt])
                                #genotype_key = tuple(sorted((allele1,allele2)))
                                if allele1 != '.' and allele2 != '.':
                                    genotype_out.write('\t'.join([genotype_data_out.get(c,'') for c in genotype_columns])+'\n')
                                    genotype_record_count += 1
                                    my.update_sql_column_spec(genotype_columns_spec, genotype_data_out)
                                if allele1 != '.' and allele1 != '':
                                    allele_counts[allele1] = allele_counts.get(allele1,0) + gt_counts[gt]
                                if allele2 != '.' and allele2 != '':
                                    allele_counts[allele2] = allele_counts.get(allele2,0) + gt_counts[gt]
                            for allele in sorted(allele_counts.keys()):
                                allele_data_out['allele'] = allele
                                allele_data_out['allele_count'] = str(allele_counts[allele])
                                allele_out.write('\t'.join([allele_data_out.get(c,'') for c in allele_columns])+'\n')
                                allele_record_count += 1
                                my.update_sql_column_spec(allele_columns_spec, allele_data_out)
    
    with open(mssql_create_allele_tables_output,'w') as o:
        #write_microsoft_sql_server_scripts(mssql_create_allele_tables_output, '[dbo].[{}]'.format(os.path.basename(allele_file)), 'allele_counts', [['CHROM','POS','REF','ALT']], allele_columns, allele_columns_spec, allele_firstrow-1)
        microsoft_sql_server_scripts = my.get_microsoft_sql_server_scripts(
            table_name='[dbo].[{}]'.format(os.path.basename(allele_file).replace('.','_')),
            index_base_name='{}_allele_counts'.format(os.path.basename(genotype_file).replace('.','_')),
            clustered_index=['CHROM','POS','REF','allele'],
            indices=[],
            columns_out=allele_columns,
            columns_out_spec=allele_columns_spec,
            rows_to_delete=allele_firstrow-1
            )
        o.write('\n\n'.join(microsoft_sql_server_scripts))
    
    with open(mssql_create_genotype_tables_output,'w') as o:
        #write_microsoft_sql_server_scripts(mssql_create_genotype_tables_output, '[dbo].[{}]'.format(os.path.basename(genotype_file)), 'genotype_counts', [['CHROM','POS','REF','allele1','allele2']], genotype_columns, genotype_columns_spec, genotype_firstrow-1)
        microsoft_sql_server_scripts = my.get_microsoft_sql_server_scripts(
            table_name='[dbo].[{}]'.format(os.path.basename(genotype_file).replace('.','_')),
            index_base_name='{}_genotype_counts'.format(os.path.basename(genotype_file).replace('.','_')),
            clustered_index=['CHROM','POS','REF','allele1','allele2'],
            indices=[],
            columns_out=genotype_columns,
            columns_out_spec=genotype_columns_spec,
            rows_to_delete=genotype_firstrow-1
            )
        o.write('\n\n'.join(microsoft_sql_server_scripts))
   
    with open(mysql_create_allele_tables_output,'w') as o:
        #write_mysql_scripts(mysql_create_allele_tables_output, mysqlimport_allele_output, '{}'.format(os.path.basename(allele_file)[:64]), 'allele_counts', [['CHROM','POS','REF','ALT']], allele_columns, allele_columns_spec, allele_firstrow)
        mysql_scripts = my.get_mysql_scripts(
            table_name=os.path.basename(allele_file).replace('.','_')[:64],
            index_base_name='allele_counts',
            indices=[['CHROM','POS','REF','allele']],
            columns_out=allele_columns,
            columns_out_spec=allele_columns_spec,
            rows_to_delete=allele_firstrow
            )
        o.write('\n\n'.join(mysql_scripts))

    with open(mysql_create_genotype_tables_output,'w') as o:
        #write_mysql_scripts(mysql_create_genotype_tables_output, mysqlimport_genotype_output, '{}'.format(os.path.basename(genotype_file)[:64]), 'genotype_counts', [['CHROM','POS','REF','allele1','allele2']], genotype_columns, genotype_columns_spec, genotype_firstrow)
        mysql_scripts = my.get_mysql_scripts(
            table_name=os.path.basename(genotype_file).replace('.','_')[:64],
            index_base_name='genotype_counts',
            indices=[['CHROM','POS','REF','allele1','allele2']],
            columns_out=genotype_columns,
            columns_out_spec=genotype_columns_spec,
            rows_to_delete=genotype_firstrow
            )
        o.write('\n\n'.join(mysql_scripts))

    print 'done'
    
def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'count alleles and genotypes in vcf files',
        epilog = '{} version {} ©2011-2012 Michael Yourshaw all rights reserved'.format(name,version))
    #input VCF parameter
    parser.add_argument('--vcfs', '-v', '-i', nargs='+', required=True,
        help='input vcf file(s); output will be <vcf>.vax')
    parser.add_argument('--output', '-o', required=True,
        help='output directory')
    parser.add_argument('--uncompressed', action='store_true', default=False,
        help='do not gzip output (default: False)')
    args = parser.parse_args()

    #dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    job_ids = run(vcfs=args.vcfs, output=args.output, uncompressed=args.uncompressed)


if __name__ == "__main__": sys.exit(main())
