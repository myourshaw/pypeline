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
import vax_merge_post_process
import math

class Vcf2VaxError(Exception): pass

name = 'pypeline.vcf2vax'
version = 71
copyright = '©2011-2013 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()


#--no_vep --no_gtx -i /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_downloaded_20120816/VCF/clinvar_00-latest.vcf.gz
#--email myourshaw@ucla.edu --vep_ini /share/apps/myourshaw/pypeline-current/vax/vep_lite.ini -v /scratch0/tmp/myourshaw/fug/fug.vcf
#--email myourshaw@ucla.edu --vep_ini /share/apps/myourshaw/pypeline-current/vax/vep.ini -v /scratch0/tmp/myourshaw/fug/fug.vcf
#-v /home/myourshaw/lab/pypeline/vax_test/kitchen_sink/kitchen_sink_GT_test.vcf
#--no_vex --no_vep -v /scratch1/tmp/myourshaw/mmjj_20130514/vcfs/vax_tmp/mmjj_20130514.hc.analysis_ready.vcf
#-v /home/myourshaw/lab/pypeline/vax_test/big-vcf/mmjj_20130514.hc.analysis_ready.vcf

def run(vcfs,
        max_part_records=20000, max_vex_parts=40, no_gtx=False, no_vex=False, no_vep=False,
        vep=None, vep_plugins=None, vep_ini=None, vep_params='',
        #vax_consequence_threshhold=None,
        deletetempfiles=False,
        python=None, perl=None,
        email=None):
        
    config = my.get_config()
        
    my_name = 'vax'
        
    #retun value: a list of job ids that can be used by downstream hold_jid
    job_ids = []
    
    vcfs_glob = vcfs
    vcfs = my.unglob(vcfs)
    if not vcfs:
        raise Vcf2VaxError("no files in {}".format(vcfs_glob))
    if not vep:
        vep = config.get('vax','vax')
    if not my.file_exists(vep):
        raise Vcf2VaxError("can't find variant effect predictor perl script {}".format(vep))
    if not vep_ini:
        vep_ini = config.get('vax','vep_ini')
    if not my.file_exists(vep_ini):
        raise Vcf2VaxError("can't find vep.ini {}".format(vep_ini))
    if not vep_plugins:
        vep_plugins = config.get('vax','vep_plugins')
    if not vep_params:
        vep_params = ''
    if not python:
        python = config.get('DEFAULT','python')
    if not perl:
        perl = config.get('DEFAULT','perl')
        
    if not no_vep and no_vex and max_part_records == None:
        raise Vcf2VaxError('Cannot run vep without creating vex and/or parts files')
        
    part_files_to_remove = []
    parts_dirs_to_delete = []
    job_dirs = []
    vex_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    gtx_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT', 'GT', 'SAMPLE', 'ALLELE1', 'PHASE', 'ALLELE2', 'AD1', 'AD2', 'DP', 'GT_SAMPLE_COUNT', 'GT_ALLELE_COUNT', 'GT_ALLELE1_COUNT', 'GT_ALLELE2_COUNT', 'GT_ALLELE1_SAMPLE_COUNT', 'GT_ALLELE2_SAMPLE_COUNT', 'GT_GENOTYPE_SAMPLE_COUNT', 'GT_ALLELE1_FREQUENCY', 'GT_ALLELE2_FREQUENCY']
    #gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    gt_sep_re = re.compile(r'([/|])')
    info_translation = maketrans(';=','__')

    print '{} {} version {}. Copyright {}.'.format(run_time, name, version, copyright)
    
    for vcf in vcfs:
        print 'reading {}'.format(vcf)
        vcf_job_ids = []
        #create vex.vcf file (CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO, with multiple ALTs expanded to separate lines)
        #and files split by chromosome with number of lines limited to max_part_records
        #for input to vax (variant effect predictor)
        #also create gtx file (SAMPLE,CHROM,POS,ID,REF,ALLELE1,PHASE,ALLELE2,QUAL,FILTER,INFO,FORMAT,GT)
        #to be joined to vax annotations in a database
        #and create sql scripts for vex and gtx tables
        
        
        #output directories and files
        output_dir = os.path.dirname(vcf)
        vcf_base = os.path.basename(vcf)
        my.makedir(output_dir)
        vex = vcf+'.vex.vcf'
        gtx = vcf+'.gtx'
        job_dir =  my.unique_dir('jobs_{}_{}'.format(vcf_base,my.localtime_squish()), output_dir)
        job_dirs.append(job_dir)
        parts_dir =  my.unique_dir('parts_{}_{}'.format(vcf_base,my.localtime_squish()), output_dir)
        if max_part_records != None:
            my.makedir(parts_dir)
            parts_dirs_to_delete.append(parts_dir)
        vax = vcf+'.vep.vax'
        processed_vax = vcf+'.vax'
            
        mssql_create_vex_tables_output = '{}.mssql_create_tables.sql'.format(vex)
        mssql_create_gtx_tables_output = '{}.mssql_create_tables.sql'.format(gtx)
        mysql_create_vex_tables_output = '{}.mysql_create_tables.sql'.format(vex)
        mysql_create_gtx_tables_output = '{}.mysql_create_tables.sql'.format(gtx)
        mysqlimport_vex_output = '{}.mysqlimport.sh'.format(vex)
        mysqlimport_gtx_output = '{}.mysqlimport.sh'.format(gtx)
        
        columns_in = {}
        sample_column_count = 0
        
        metadata_header = [] # input lines that start with #
        vex_metadata = [
        '##INFO=<ID=oALT,Number=1,Type=String,Description="Original VCF file\'s ALT field, which may contain multiple alleles (the ALT column in this file will have only one allele per record)">',
        '## VEX file created at {} by {} version {}. Copyright {}.'.format(run_time, name, version, copyright),
        '## VEX file has one line per ALT allele and omits any FORMAT and genotype columns that were in the original VCF file',
        '## VEX --vcf {}'.format(vcf),
        '## VEX --vex {}'.format(vex),
        ]
        gtx_metadata = [
        '##INFO=<ID=oALT,Number=1,Type=String,Description="Original VCF file ALT field, which may contain multiple alleles (same as ALT field in this file)">',
        '##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample ID taken from genotype column heading in VCF file">',
        '## GTX file created at {} by {} version {}. Copyright {}.'.format(run_time, name, version, copyright),
        '## GTX file has one line per sample and genotype statistics for all samples',
        '## GTX --vcf {}'.format(vcf),
        '## GTX --gtx {}'.format(gtx),
        '## AD1 ALLELE1 depth (VCF file)',
        '## AD2 ALLELE2 depth (VCF file)',
        '## DP Approximate read depth (reads with MQ=255 or with bad mates are filtered) (VCF file)',
        '## GT_SAMPLE_COUNT number of samples  in called genotypes (VCF file)',
        '## GT_ALLELE_COUNT number of alleles in called genotypes (VCF file)',
        '## GT_ALLELE1_COUNT number of allele1s in called genotypes (VCF file)',
        '## GT_ALLELE2_COUNT number of allele2s in called genotypes (VCF file)',
        '## GT_ALLELE1_SAMPLE_COUNT number of samples with allele1 in called genotypes (VCF file)',
        '## GT_ALLELE2_SAMPLE_COUNT number of samples with allele2 in called genotypes (VCF file)',
        '## GT_GENOTYPE_SAMPLE_COUNT number of samples with this genotype in called genotypes (VCF file)',
        '## GT_ALLELE1_FREQUENCY allele1 frequency in called genotypes (VCF file)',
        '## GT_ALLELE2_FREQUENCY allele2 frequency in called genotypes (VCF file)',
        ]
        this_chrom = None
        this_part_fh = None
        this_part_number = 0
        this_part_records = 0
        #list of split vcfs for input to vax
        part_files = []
        
        if max_vex_parts and not no_vep:
            vcf_file_len = sum(1 for line in my.open_gz_or_text(vcf))
            max_part_records = max(math.ceil(float(vcf_file_len)/float(max_vex_parts)), max_part_records)
        this_part_record_count = 0
        
        with my.open_gz_or_text(vcf) as vcf_in: #, open(vex, 'w') as vex_out:
            if not no_vex:
                vex_out = open(vex, 'w')
            vex_columns_spec = {c: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for c in vex_columns}
            gtx_columns_spec = {c: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for c in gtx_columns}
            vcf_line_count = 0
            vcf_record_count = 0
            vex_record_count = 0
            gtx_record_count = 0
            vex_firstrow = 0
            gtx_firstrow = 0
            for line in vcf_in:
                ##DEBUG
                #if vcf_record_count > 100: break
                ##DEBUG
                vcf_line_count += 1
                if line.startswith('##'):
                    metadata_header.append(line)
                    if not no_vex:
                        vex_out.write(line)
                        vex_firstrow += 1
                    continue
                elif line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                    col_list = line.rstrip('\n').split('\t')
                    has_genotypes= line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t') and len(col_list) > 9 and col_list[9] != 'Uploaded_variation'
                    columns_in = {col_list[x]: x for x in range(len(col_list))}
                    if not no_vex:
                        print 'writing {}'.format(vex)
                        vex_out.write('\n'.join(vex_metadata)+'\n')
                        vex_firstrow += len(vex_metadata)
                        vex_out.write('\t'.join(vex_columns)+'\n')
                        vex_firstrow += 1
                    if not no_gtx and has_genotypes:
                        print 'writing {}'.format(gtx)
                        gtx_out = open(gtx, 'w')
                        gtx_out.write(''.join(metadata_header))
                        gtx_firstrow += len(metadata_header)
                        gtx_out.write('\n'.join(gtx_metadata)+'\n')
                        gtx_firstrow += len(gtx_metadata)
                        gtx_out.write('\t'.join(gtx_columns)+'\n')
                        gtx_firstrow += 1
                        sample_columns = col_list[9:columns_in.get('Uploaded_variation',len(col_list))]
                        if len(sample_columns) != len(set(sample_columns)):
                            raise Vcf2VaxError('duplicated sample column names in vcf file {} [{}]'.format(','.join(vcf, sample_columns)))
                elif not columns_in:
                    raise Vcf2VaxError('data encountered before column names at line {}'.format(vcf_line_in_count))
                else:
                    #data record
                    vcf_record_count += 1
                    data_list = ['' if d == '-' else d for d in line.rstrip('\n').split('\t')]
                    data_in = {k: data_list[columns_in[k]] for k in columns_in.keys()}
                    ref = data_in['REF']
                    alt_list = data_in['ALT'].rstrip(',').split(',')
                    these_alleles = [ref] + alt_list
                    #vex files expand ALT to one line for each alt allele
                    for this_alt in alt_list:
                        vex_data_out = {c: data_in[c] for c in vex_columns}
                        vex_data_out['ALT'] = this_alt
                        vex_data_out['INFO'] = 'oALT={};{}'.format(data_in['ALT'], data_in['INFO'])
                        if not no_vex:
                            vex_out.write('\t'.join([vex_data_out[c] for c in vex_columns])+'\n')
                            vex_record_count += 1
                            accumulate_sql_spec(vex_data_out, vex_columns_spec)
                        if not no_vep and (max_vex_parts or max_part_records != None):
                            if (max_vex_parts and (this_part_record_count == 0 or this_part_record_count >= max_part_records)) \
                            or ( not max_vex_parts and (vex_data_out['#CHROM'] != this_chrom or (max_part_records != 0 and this_part_record_count >= max_part_records))):
                                if this_part_fh:
                                    this_part_fh.close()
                                this_chrom = vex_data_out['#CHROM']
                                this_part_file = os.path.join(parts_dir,'{:0=8}.{}.{}'.format(this_part_number, this_chrom,os.path.basename(vex)))
                                part_files.append(this_part_file)
                                part_files_to_remove.append(this_part_file)
                                this_part_fh = open(this_part_file, 'w')
                                this_part_number += 1
                                this_part_record_count = 0
                                this_part_fh.write(''.join(metadata_header))
                                this_part_fh.write('\t'.join(vex_columns)+'\n')
                            this_part_fh.write('\t'.join([vex_data_out[c] for c in vex_columns])+'\n')
                            this_part_record_count += 1
                        elif not no_vep:
                            part_files.append(vex)
                            
                    #gtx file expands GT columns to one line per sample
                    if has_genotypes and not no_gtx:
                        format_split = [f.upper() for f in data_in['FORMAT'].split(':')]
                        if format_split[0] == 'GT':
                            AD = format_split.index('AD') if 'AD' in format_split else None
                            DP = format_split.index('DP') if 'DP' in format_split else None
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
                            for gt_info in [data_in[s].split(':') for s in sample_columns]:
                                gt = gt_info[0]
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
                                #gtx_data_out['ALT'] = this_alt
                                gtx_data_out['SAMPLE'] = this_sample
                                gtx_data_out['INFO'] = 'oALT={};SAMPLE={};{}'.format(data_in['ALT'], this_sample.translate(info_translation), data_in['INFO'])
                                gtx_data_out['GT'] = data_in[this_sample]
                                gt_split = data_in[this_sample].split(':')
                                gt = gt_split[0]
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
                                    #gtx_data_out['ZYGOSITY'] = '0' if allele1 == ref and (allele2 == ref or allele2 == '') else '2' if allele1 == alt and (allele2 == alt or allele2 == '') else '' if allele1 == '.' or allele1 == '' or allele2 == '.' or allele2 == '' else '1' if allele1 != allele2 and allele1_index < allele2_index else ''
                                    ad = gt_split[AD].split(',') if AD != None else None
                                    gtx_data_out['AD1'] = ad[0] if ad and len(ad) > 0 and ad[0] != '.' else ''
                                    gtx_data_out['AD2'] = ad[1] if ad and len(ad) > 1 and ad[1] != '.' else ''
                                    gtx_data_out['DP'] = gt_split[DP] if DP != None else ''
                                    if gtx_data_out['DP'] == '' and my.is_int(gtx_data_out['AD1']) and my.is_int(gtx_data_out['AD2']):
                                        gtx_data_out['DP'] = str(int(gtx_data_out['AD1']) + int(gtx_data_out['AD2']))
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
                                else:
                                    gtx_data_out['AD1'] = ''
                                    gtx_data_out['AD2'] = ''
                                    gtx_data_out['DP'] = ''
                                    gtx_data_out['GT_SAMPLE_COUNT'] = ''
                                    gtx_data_out['GT_ALLELE_COUNT'] = ''
                                    gtx_data_out['GT_ALLELE1_COUNT'] = ''
                                    gtx_data_out['GT_ALLELE2_COUNT'] = ''
                                    gtx_data_out['GT_ALLELE1_SAMPLE_COUNT'] = ''
                                    gtx_data_out['GT_ALLELE2_SAMPLE_COUNT'] = ''
                                    gtx_data_out['GT_GENOTYPE_SAMPLE_COUNT'] = ''
                                    gtx_data_out['GT_ALLELE1_FREQUENCY'] = ''
                                    gtx_data_out['GT_ALLELE2_FREQUENCY'] = ''
                                gtx_out.write('\t'.join([gtx_data_out.get(c,'') for c in gtx_columns])+'\n')
                                gtx_record_count += 1
                                accumulate_sql_spec(gtx_data_out, gtx_columns_spec)
        if not no_vex:
            vex_out.close()
        if  not no_gtx and has_genotypes:
            gtx_out.close()
    
        if not no_vex:
            write_microsoft_sql_server_scripts(mssql_create_vex_tables_output, '[dbo].[{}]'.format(os.path.basename(vex).replace('.','_')), os.path.basename(vex).replace('.','_'), ['CHROM','POS','REF','ALT'], [['ID']], vex_columns, vex_columns_spec, vex_firstrow-1)
            write_mysql_scripts(mysql_create_vex_tables_output, mysqlimport_vex_output, '{}'.format(os.path.basename(vex).replace('.','_')[:64]), os.path.basename(vex).replace('.','_')[:16], (['CHROM','POS','REF','ALT'],['ID']), vex_columns, vex_columns_spec, vex_firstrow)
        
        if  not no_gtx and has_genotypes:
            write_microsoft_sql_server_scripts(mssql_create_gtx_tables_output, '[dbo].[{}]'.format(os.path.basename(gtx).replace('.','_')), os.path.basename(gtx).replace('.','_'), ['SAMPLE','CHROM','POS','REF','ALT','ALLELE1','ALLELE2'], [['CHROM','POS','REF','ALT','ALLELE1','ALLELE2']], gtx_columns, gtx_columns_spec, gtx_firstrow-1)
            write_mysql_scripts(mysql_create_gtx_tables_output, mysqlimport_gtx_output, '{}'.format(os.path.basename(gtx).replace('.','_')[:64]), os.path.basename(gtx).replace('.','_')[:16], [['SAMPLE','CHROM','POS','REF','ALT','ALLELE1','ALLELE2'], ['CHROM','POS','REF','ALT','ALLELE1','ALLELE2']], gtx_columns, gtx_columns_spec, gtx_firstrow)
            
        #finished this vcf
        print 'input vcf file {} had {} lines, {} records.'.format(vcf, vcf_line_count, vcf_record_count)
        if not no_vex:
            print 'output vex.vcf file {} had {} records.'.format(vex, vex_record_count)
        if not no_gtx and has_genotypes:
            print 'output gtx file {} had {} records.'.format(gtx, gtx_record_count)
        if not no_vep:
            print 'temporary files to be annotated are in {}'.format(parts_dir)
        
        #annotate each of the split vex.vcf parts with variant effect predictor
        if not no_vep:
            cmds = []
            vax_parts = []
            for vcf_part in part_files:
                vax_part = vcf_part+'.vax'
                vax_parts += [vax_part]
                cmds.append('{} {} -input_file {} -output_file {} --config {} {}'.format(perl, vep, vcf_part, vax_part, vep_ini, vep_params))
            part_files_to_remove += vax_parts
            job_name = 'vax_'+vcf_base
            job = my.run_job(cmds, job_name, job_dir)
            vax_jobid = job.jobId
            vcf_job_ids.append(job.jobId)
            print 'submitted job {} to annotate {}'.format(vax_jobid, vcf)
            
            #merge and post-process annotated vax files
            cmd= '{} {} -i {} -o {}'.format(
                python, vax_merge_post_process.__file__, ' '.join(vax_parts), processed_vax)
                #python, vax_merge_post_process.__file__, vax_parts+'*.*.*.vax', processed_vax)
            job_name = 'vax_merge_'+vcf_base
            job = my.run_job(cmd, job_name, job_dir, hold_jid=vax_jobid, email=email)
            merge_vax_jobid = job.jobId
            vcf_job_ids.append(job.jobId)
            print 'submitted job {} to merge and post-process annotations into {}'.format(merge_vax_jobid, processed_vax)
            
            job_ids += vcf_job_ids

    #remove temporary files
    if deletetempfiles:
        cmd = 'rm -rf {}'.format(' '.join(parts_dirs_to_delete))
        job_name = 'vax_rm_'+vcf_base
        job = my.run_job(cmd, job_name, job_dir, hold_jid=job_ids)
        rm_jobid = job.jobId
        vcf_job_ids.append(job.jobId)
        job_ids.append(job.jobId)
        print 'submitted job {} to remove temporary files from {}'.format(rm_jobid, ' '.join(parts_dirs_to_delete))
    else:
        print 'when all jobs are done, please manually remove temporary files from {}\nin future you can use the --deletetempfiles option to auto-delete the parts_* directories'.format(' '.join(parts_dirs_to_delete))

    print 'when all jobs are done, check STDERR output in {} to validate the run'.format(' '.join(job_dirs))
    return job_ids

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
def write_microsoft_sql_server_scripts(output_file, table_name, index_base_name, clustered_index, indices, columns_out, columns_out_spec, rows_to_delete):
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
        #clustered_index is a list of columns
        if clustered_index:
            ix_name = 'IX_'+index_base_name+'_'+'_'.join(clustered_index)
            ix_cols = ',\n\t'.join(['['+c+'] ASC' for c in clustered_index])
            ms.write("""CREATE CLUSTERED INDEX [{0}] ON {1}
(
    {2}
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO
""".format(ix_name, table_name, ix_cols))
        #indices is a list of indices containing lists of columns
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
        description = 'annotate one or more vcf files with Ensembl Variant Effect Predictor + extras',
        epilog = '{} version {} {}'.format(name, version, copyright))
    #input VCF parameter
    parser.add_argument('--vcfs', '-v', '-i', nargs='+', required=True,
        help='input vcf file(s); output will be <vcf>.vax')
    #vcf2vex parameters
    parser.add_argument('--max_part_records', type=int, default=20000,
        help='maximum number of records per vex chromosome part file; if None do not create parts; if 0 split on chromosome only (default:20000)')
    parser.add_argument('--max_vex_parts', type=int, default=40,
        help='maximum number of vex parts (used to limit number of cluster nodes); overrides max_part_records (default:40)')
    parser.add_argument('--no_gtx', action='store_true', default=False,
        help='do not create gtx files (default:False)')
    parser.add_argument('--no_vex', action='store_true', default=False,
        help='do not create vex files (default:False)')
    parser.add_argument('--no_vep', action='store_true', default=False,
        help='do not run variant effect predictor, just create gtx and or vex files files (default:False)')
    #variant effect predictor parameters
    parser.add_argument('--vep',
        help='Ensembl variant effect predictor (default: config->vax.vax)')
    parser.add_argument('--vep_plugins',
        help='path to VEP_plugins directory (default: config->vax.vep_plugins)')
    parser.add_argument('--vep_ini',
        help='path to vep.ini, which contains VEP parameters and plugins to run (default: config->vax.vep_ini)')
    parser.add_argument('--vep_params',
        help='a single quoted string with VEP parameters to override vep_ini (default: None)')
    #output parameters
    #parser.add_argument('--vax_consequence_threshhold', type=int, default=None,
    #    help='only output revised_vax records where Consequence_rank <= vax_consequence_threshhold (default: None)')
    parser.add_argument('--deletetempfiles', action='store_true', default=False,
        help='delete temporary files (default: False)')
    #executables
    parser.add_argument('--python',
        help='path to python executable (default: config->python)')
    parser.add_argument('--perl',
        help='path to perl executable (default: config->perl)')
    args = parser.parse_args()

    #dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    job_ids = run(vcfs=args.vcfs,
        max_part_records=args.max_part_records, max_vex_parts=args.max_vex_parts, no_gtx=args.no_gtx, no_vex=args.no_vex, no_vep=args.no_vep,
        vep=args.vep, vep_plugins=args.vep_plugins, vep_ini=args.vep_ini, vep_params=args.vep_params,
        #vax_consequence_threshhold=args.vax_consequence_threshhold,
        deletetempfiles=args.deletetempfiles,
        python=args.python, perl=args.perl,
        email=args.email)


if __name__ == "__main__": sys.exit(main())
