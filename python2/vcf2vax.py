#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import gzip
#from math import ceil
import re
from string import maketrans
from warnings import warn
import my
import job
import sql_columns
import vax_merge_post_process

class Vcf2VaxError(Exception): pass

name = 'pypeline.vcf2vax'
version = 75
copyright = '©2011-2014 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()


#--no_vax --no_gtx -i /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_downloaded_20120816/VCF/clinvar_00-latest.vcf.gz
#--email myourshaw@ucla.edu --vep_ini /share/apps/myourshaw/pypeline-current/vax/vep_lite.ini -v /scratch0/tmp/myourshaw/fug/fug.vcf
#--email myourshaw@ucla.edu --vep_ini /share/apps/myourshaw/pypeline-current/vax/vep.ini -v /scratch0/tmp/myourshaw/fug/fug.vcf
#-v /home/myourshaw/lab/pypeline/vax_test/kitchen_sink/kitchen_sink_GT_test.vcf
#--no_vex --no_vax -v /scratch1/tmp/myourshaw/mmjj_20130514/vcfs/vax_tmp/mmjj_20130514.hc.analysis_ready.vcf
#-v /home/myourshaw/lab/pypeline/vax_test/big-vcf/mmjj_20130514.hc.analysis_ready.vcf
#-v /share/apps/myourshaw/vax/test/t.vcf
#-v /scratch1/tmp/myourshaw/schwa/combined-06052014.schwa.hc.analysis_ready.vcf --no_vax
#--no_vax --no_vex -v /share/apps/myourshaw/vax/test/multialleic_test.vcf

def run(vcfs,
        max_part_records=1000, max_simultaneous_jobs=20, no_gtx=False, no_vex=False, no_vax=False,
        vep=None, vep_plugins=None, vep_ini=None, vep_params='',
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
        vep = config.get('vax','vep_pl')
    if not my.file_exists(vep):
        raise Vcf2VaxError("can't find variant effect predictor perl script {}".format(vep))
    if not no_vax and not vep_ini:
        vep_ini = config.get('vax','vep_no_Alignment_ini')
    if not no_vax and not my.file_exists(vep_ini):
        raise Vcf2VaxError("can't find vep.ini {}".format(vep_ini))
    if not no_vax and not vep_plugins:
        vep_plugins = config.get('vax','vep_plugins')
    if not vep_params:
        vep_params = ''
    if not python:
        python = config.get('DEFAULT','python')
    if not perl:
        perl = config.get('DEFAULT','perl')
    
    #enforce 1-40 simultaneous jobs
    max_simultaneous_jobs = int(max_simultaneous_jobs) if max_simultaneous_jobs and int(max_simultaneous_jobs) >= 1 and int(max_simultaneous_jobs) <=40 else 20
    #enforce 100-5000 records per job
    max_part_records = int(max_part_records) if max_part_records and int(max_part_records) >= 1 and int(max_part_records) <= 5000 else 1000
    
    print '{} {} version {}. Copyright {}.'.format(run_time, name, version, copyright)
    print """IMPORTANT:
This application will create hidden files with names like '.filename.done'.
After a restart, files with an associated .done file will not be recreated.
To force recreation of an existing file, delete the associated .done file.
"""
    
    parts_dirs = []
    job_dirs = []
    all_vcfs_all_vcf_parts = []
    all_vcfs_all_vax_parts = []
    all_vcfs_vax_parts_to_vep = []
    all_vcfs_all_vax = []
    all_vax_job_ids = []
    job_ids = []
    vex_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    gtx_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT', 'GT', 'SAMPLE', 'ALLELE1', 'PHASE', 'ALLELE2', 'ADREF', 'AD1', 'AD2', 'DP', 'GQ', 'PLRR', 'PLRA1', 'PLA1A1', 'PLRA2', 'PLA1A2', 'PLA2A2', 'GT_SAMPLE_COUNT', 'GT_ALLELE_COUNT', 'GT_ALLELE1_COUNT', 'GT_ALLELE2_COUNT', 'GT_ALLELE1_SAMPLE_COUNT', 'GT_ALLELE2_SAMPLE_COUNT', 'GT_GENOTYPE_SAMPLE_COUNT', 'GT_ALLELE1_FREQUENCY', 'GT_ALLELE2_FREQUENCY']
    gt_sep_re = re.compile(r'([/|])')
    info_translation = maketrans(';=','__')
    
    if len(vcfs) > 1:
        warn("""Each vcf file will be processed separately to produce unmerged vex, gx, and vax files.
This may not be what you want.""")
        
    for vcf in vcfs:
        
        #create vex.vcf file (CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO, with multiple ALTs expanded to separate lines)
        #and files split by chromosome with number of lines limited to max_part_records
        #for input to vax (variant effect predictor)
        #also create gtx file (SAMPLE,CHROM,POS,ID,REF,ALLELE1,PHASE,ALLELE2,QUAL,FILTER,INFO,FORMAT,GT)
        #to be joined to vax annotations in a database
        #and create sql scripts for vex, gtx and vax tables
        
        this_vcf_job_ids = []
        
        #output directories and files
        output_dir = os.path.dirname(vcf)
        vcf_base = os.path.basename(vcf)
        my.makedir(output_dir)
        
        tmp_dir = os.path.join(output_dir, 'tmp_vax_'+vcf_base)
        my.makedir(tmp_dir)
        
        job_dir = os.path.join(tmp_dir, 'jobs')
        my.makedir(job_dir)
        job_dirs.append(job_dir)
        
        parts_dir = os.path.join(tmp_dir, 'parts')
        my.makedir(parts_dir)
        parts_dirs.append(parts_dir)
        
        vex = my.swap_ext(vcf, '.vcf', '.vex')
        vex_done = my.done_file(vex)
        if no_vex or (my.file_exists(vex) and my.file_exists(vex_done)):
            create_vex = False
            print ('{} already created. Skipping.'.format(vex))
        else:
            create_vex = True

        gtx = my.swap_ext(vcf, '.vcf', '.gtx')
        gtx_done = my.done_file(gtx)
        if no_gtx or (my.file_exists(gtx) and my.file_exists(gtx_done)):
            create_gtx = False
            print ('{} already created. Skipping.'.format(gtx))
        else:
            create_gtx = True

        parts_done = my.done_file(os.path.join(parts_dir, 'all_parts'))
        if no_vax or(my.file_exists(parts_done)):
            create_parts = False
        else:
            create_parts = True

        vax = vcf+'.vax'
        vax_done = my.done_file(vax)
        if no_vax or (my.file_exists(vax) and my.file_exists(vax_done)):
            create_vax = False
            print ('{} already created. Skipping.'.format(vax))
        else:
            create_vax = True

        if not create_vex and not create_gtx and not create_parts and not create_vax:
            warn('Nothing to do for VCF file {}. If you are trying to re-reun vax, you may neeed to delete one or more .done files.'.format(vcf))
            next
        
        vex_microsoft_sql_file = vex+'.microsoft.sql'
        vex_mysql_file = vex+'.mysql'
        vex_microsoft_done_file = my.done_file(vex_microsoft_sql_file)
        vex_mysql_done_file = my.done_file(vex_mysql_file)

        gtx_microsoft_sql_file = gtx+'.microsoft.sql'
        gtx_mysql_file = gtx+'.mysql'
        gtx_microsoft_done_file = my.done_file(gtx_microsoft_sql_file)
        gtx_mysql_done_file = my.done_file(gtx_mysql_file)

        vax_microsoft_sql_file = vax+'.microsoft.sql'
        vax_mysql_file = vax+'.mysql'
        vax_microsoft_done_file = my.done_file(vax_microsoft_sql_file)
        vax_mysql_done_file = my.done_file(vax_mysql_file)

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
        '## ADREF REF depth (VCF file)',
        '## AD1 ALLELE1 depth (VCF file)',
        '## AD2 ALLELE2 depth (VCF file)',
        '## DP Approximate read depth (reads with MQ=255 or with bad mates are filtered) (VCF file)',
        '## GQ Genotype Quality (VCF file)',
        '## PLRR Normalized, Phred-scaled likelihood for genotype REF/REF given the set of alleles defined in the REF and ALT fields. (VCF file)',
        '## PLRA1 Normalized, Phred-scaled likelihood for genotype REF/ALLELE1 given the set of alleles defined in the REF and ALT fields. (VCF file)',
        '## PLA1A1 Normalized, Phred-scaled likelihood for genotype ALLELE1/ALLELE1 given the set of alleles defined in the REF and ALT fields. (VCF file)',
        '## PLRA2 Normalized, Phred-scaled likelihood for genotype REF/ALLELE2 given the set of alleles defined in the REF and ALT fields. (VCF file)',
        '## PLA1A2 Normalized, Phred-scaled likelihood for genotype ALLELE1/ALLELE2 given the set of alleles defined in the REF and ALT fields. (VCF file)',
        '## PLA2A2 Normalized, Phred-scaled likelihood for genotype ALLELE2/ALLELE2 given the set of alleles defined in the REF and ALT fields. (VCF file)',
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
        
        print 'Reading {}'.format(vcf)
        with my.open_gz_or_text(vcf) as vcf_in:

            if create_vex:
                vex_out = open(vex, 'w')
                
            vcf_line_count = 0
            vcf_record_count = 0
            vex_record_count = 0
            gtx_record_count = 0
            vex_firstrow = 0
            gtx_firstrow = 0
            this_chrom = None
            that_vex_chrom = None
            this_part_fh = None
            this_part_number = 0
            this_part_records = 0
            this_part_record_count = 0
            #list of split vcfs for input to vax
            this_vcf_part_files_to_vep = []
            this_vcf_all_part_files = []
        
            if create_vex or create_gtx or create_parts or create_vax:
                for line in vcf_in:
                    ##DEBUG
                    #if vcf_record_count > 100: break
                    ##DEBUG
                    vcf_line_count += 1
                    #metadata
                    if line.startswith('##'):
                        metadata_header.append(line)
                        if create_vex:
                            vex_out.write(line)
                            vex_firstrow += 1
                        continue
                    #column headers
                    elif line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                        col_list = line.rstrip('\n').split('\t')
                        has_genotypes= line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t') and len(col_list) > 9 and col_list[9] != 'Uploaded_variation'
                        if not has_genotypes:
                            create_gtx = False
                            if not(create_vex or create_gtx or create_parts or create_vax):
                                warn('Nothing to do for VCF file {}. If you are trying to re-reun vax, you may neeed to delete one or more .done files.'.format(vcf))
                                break
                        columns_in = {col_list[x]: x for x in range(len(col_list))}
                        #vex metadata and column header
                        if create_vex:
                            print 'Writing {}. This can take a long time for large VCFs.'.format(vex)
                            vex_out.write('\n'.join(vex_metadata)+'\n')
                            vex_firstrow += len(vex_metadata)
                            vex_out.write('\t'.join(vex_columns)+'\n')
                            vex_firstrow += 1
                        #gtx metadata and column header
                        if create_gtx:
                            print 'Writing {}. This can take a long time for large VCFs.'.format(gtx)
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
                        raise Vcf2VaxError('Data encountered before column names at line {}'.format(vcf_line_count))
                    else:
                        #data record
                        vcf_record_count += 1
                        data_list = ['' if d == '-' else d for d in line.rstrip('\n').split('\t')]
                        data_in = {k: data_list[columns_in[k]] for k in columns_in.keys()}
                        ref = data_in['REF']
                        alt_list = data_in['ALT'].rstrip(',').split(',')
                        these_alleles = [ref] + alt_list
                        these_alleles_indices = {i: these_alleles[i] for i in range(len(these_alleles))}
                        
                        #vex files expand ALT to one line for each alt allele
                        for this_alt in alt_list:
                            vex_data_out = {c: data_in[c] for c in vex_columns}
                            if vex_data_out['#CHROM'] != that_vex_chrom:
                                that_vex_chrom = vex_data_out['#CHROM']
                                print('contig {}'.format(that_vex_chrom))
                            vex_data_out['ALT'] = this_alt
                            vex_data_out['INFO'] = 'oALT={};{}'.format(data_in['ALT'], data_in['INFO'])
                            if create_vex:
                                vex_out.write('\t'.join([vex_data_out[c] for c in vex_columns])+'\n')
                                vex_record_count += 1
                                
                            #are we making vex part files for input to vep?
                            if create_parts:
                                #start a new part file?
                                if (
                                    this_part_record_count == 0
                                    or this_part_record_count >= max_part_records
                                    or vex_data_out['#CHROM'] != this_chrom
                                ):
                                    #close the previous part file
                                    if this_part_fh:
                                        this_part_fh.close()
                                        with open(this_part_done, 'w'):
                                            pass
                                    #create the new part file, if necessary
                                    this_chrom = vex_data_out['#CHROM']
                                    this_part_file = os.path.join(parts_dir,'{:0=8}.{}.{}'.format(this_part_number, this_chrom, os.path.basename(vex)))
                                    this_part_done = my.done_file(this_part_file)
                                    all_vcfs_all_vcf_parts.append(this_part_file)
                                    this_vcf_all_part_files.append(this_part_file)
                                    all_vcfs_all_vcf_parts.append(this_part_file)
                                    this_part_number += 1
                                    this_part_record_count = 0
                                    if my.file_exists(this_part_file) and my.file_exists(this_part_done):
                                        print ('{} already created. Skipping.'.format(this_part_file))
                                        skipping_part = True
                                    else:
                                        print(this_part_file)
                                        skipping_part = False
                                        this_vcf_part_files_to_vep.append(this_part_file)
                                        this_part_fh = open(this_part_file, 'w')
                                        this_part_fh.write(''.join(metadata_header))
                                        this_part_fh.write('\t'.join(vex_columns)+'\n')
                                #part data record
                                this_part_record_count += 1
                                if not skipping_part:
                                    this_part_fh.write('\t'.join([vex_data_out[c] for c in vex_columns])+'\n')
                
                            #gtx file expands GT columns to one line per sample
                            if create_gtx:
                                format_split = [f.upper() for f in data_in['FORMAT'].split(':')]
                                #ignore records with no GT (shouldn't happen with GATK-produced VCFs)
                                if format_split[0] == 'GT':
                                    AD = format_split.index('AD') if 'AD' in format_split else None
                                    DP = format_split.index('DP') if 'DP' in format_split else None
                                    GQ = format_split.index('GQ') if 'GQ' in format_split else None
                                    PL = format_split.index('PL') if 'PL' in format_split else None
                                    #calculate statistics for all samples at this locus
                                    sample_count = 0
                                    allele_count = 0
                                    allele_counts = {}
                                    allele_sample_counts = {}
                                    sample_genotype_counts = {}
                                    #write lines for each sample at this locus
                                    for this_sample in sample_columns:
                                        gtx_data_out = {c: data_in.get(c,'') for c in gtx_columns}
                                        #gtx_data_out['ALT'] = this_alt
                                        gtx_data_out['SAMPLE'] = this_sample
                                        gtx_data_out['INFO'] = 'oALT={};SAMPLE={};{}'.format(data_in['ALT'], this_sample.translate(info_translation), data_in['INFO'])
                                        gtx_data_out['GT'] = data_in[this_sample]
                                        this_gt_data = data_in[this_sample].split(':')
                                        gt = this_gt_data[0]
                                        this_sample_allele_indices = gt_sep_re.split(gt)
                                        allele1_index_str = this_sample_allele_indices[0] if len(this_sample_allele_indices) > 0 else ''
                                        allele1_index = int(allele1_index_str) if my.is_int(allele1_index_str) else None
                                        phase = this_sample_allele_indices[1] if len(this_sample_allele_indices) > 1 else ''
                                        allele2_index_str = this_sample_allele_indices[2] if len(this_sample_allele_indices) > 2 else ''
                                        allele2_index = int(allele2_index_str) if my.is_int(allele2_index_str) else None
                                        allele1 = '.' if allele1_index_str == '.' else these_alleles[allele1_index] if allele1_index != None else ''
                                        gtx_data_out['ALLELE1'] = allele1
                                        allele2 = '.' if allele2_index_str == '.' else these_alleles[allele2_index] if allele2_index_str != None else ''
                                        gtx_data_out['ALLELE2'] = allele2
                                        gtx_data_out['PHASE'] = phase
                                        #allele and genotype statistics for all samples in VCF
                                        genotype_key = tuple(sorted((allele1,allele2)))
                                        if allele1 != '.' and allele2 != '.':
                                            # genotype_key = tuple(sorted((allele1,allele2)))
                                            #gtx_data_out['ZYGOSITY'] = '0' if allele1 == ref and (allele2 == ref or allele2 == '') else '2' if allele1 == alt and (allele2 == alt or allele2 == '') else '' if allele1 == '.' or allele1 == '' or allele2 == '.' or allele2 == '' else '1' if allele1 != allele2 and allele1_index < allele2_index else ''
                                            sample_count+=1
                                            allele_count+=2
                                            allele_counts[allele1] = allele_counts.get(allele1,0) + 1
                                            allele_counts[allele2] = allele_counts.get(allele2,0) + 1
                                            allele_sample_counts[allele1] = allele_sample_counts.get(allele1,0) + 1
                                            if allele1 != allele2:
                                                allele_sample_counts[allele2] = allele_sample_counts.get(allele2,0) + 1
                                            sample_genotype_counts[genotype_key] = sample_genotype_counts.get(genotype_key,0) + 1
                                            #allele depths
                                            locus_ADs = this_gt_data[AD].split(',') if AD != None else None
                                            gtx_data_out['ADREF'] = locus_ADs[0] if locus_ADs and len(locus_ADs) > 0 and locus_ADs[0] != '.' else ''
                                            gtx_data_out['AD1'] = locus_ADs[allele1_index] if locus_ADs and allele1_index != None and len(locus_ADs) > allele1_index and locus_ADs[allele1_index] != '.' else ''
                                            gtx_data_out['AD2'] = locus_ADs[allele2_index] if locus_ADs and allele2_index != None and len(locus_ADs) > allele2_index and locus_ADs[allele2_index] != '.' else ''
                                            gtx_data_out['DP'] = this_gt_data[DP] if DP != None else ''
                                            gtx_data_out['GQ'] = this_gt_data[GQ] if GQ != None else ''

                                            #phred-scaled genotype likelihoods
                                            #the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j
                                            locus_PLs = map(int, this_gt_data[PL].split(',')) if PL != None else None
                                            if locus_PLs:
                                                locusPLRR = locus_PLs[0]
                                                gtx_data_out['PLRR'] = locusPLRR
                                                if len(locus_PLs) > allele1_index:
                                                    locusPLRA1 = locus_PLs[(allele1_index*(allele1_index+1)/2)+0]
                                                    locusPLA1A1 = locus_PLs[(allele1_index*(allele1_index+1)/2)+allele1_index]
                                                    gtx_data_out['PLRA1'] = locusPLRA1
                                                    gtx_data_out['PLA1A1'] = locusPLA1A1
                                                else:
                                                    gtx_data_out['PLRA1'] = ''
                                                    gtx_data_out['PLA1A1'] = ''
                                                if len(locus_PLs) > allele2_index:
                                                    locusPLRA2 = locus_PLs[(allele2_index*(allele2_index+1)/2)+0]
                                                    locusPLA2A2 = locus_PLs[(allele2_index*(allele2_index+1)/2)+allele2_index]
                                                    gtx_data_out['PLRA2'] = locusPLRA2
                                                    gtx_data_out['PLA2A2'] = locusPLA2A2
                                                else:
                                                    gtx_data_out['PLRA2'] = ''
                                                    gtx_data_out['PLA2A2'] = ''
                                                if len(locus_PLs) > max(allele1_index,allele2_index):
                                                    locusPLA1A2 = locus_PLs[(allele2_index*(allele2_index+1)/2)+allele1_index]
                                                    gtx_data_out['PLA1A2'] = locusPLA1A2
                                                else:
                                                    gtx_data_out['PLA1A2'] = ''
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
                                            gtx_data_out['ADREF'] = ''
                                            gtx_data_out['AD1'] = ''
                                            gtx_data_out['AD2'] = ''
                                            gtx_data_out['DP'] = ''
                                            gtx_data_out['GQ'] = ''
                                            gtx_data_out['PLRR'] = ''
                                            gtx_data_out['PLRA1'] = ''
                                            gtx_data_out['PLA1A1'] = ''
                                            gtx_data_out['PLRA2'] = ''
                                            gtx_data_out['PLA1A2'] = ''
                                            gtx_data_out['PLA2A2'] = ''
                                            gtx_data_out['GT_SAMPLE_COUNT'] = ''
                                            gtx_data_out['GT_ALLELE_COUNT'] = ''
                                            gtx_data_out['GT_ALLELE1_COUNT'] = ''
                                            gtx_data_out['GT_ALLELE2_COUNT'] = ''
                                            gtx_data_out['GT_ALLELE1_SAMPLE_COUNT'] = ''
                                            gtx_data_out['GT_ALLELE2_SAMPLE_COUNT'] = ''
                                            gtx_data_out['GT_GENOTYPE_SAMPLE_COUNT'] = ''
                                            gtx_data_out['GT_ALLELE1_FREQUENCY'] = ''
                                            gtx_data_out['GT_ALLELE2_FREQUENCY'] = ''
                                        gtx_out.write('\t'.join(map(str,[gtx_data_out.get(c,'') for c in gtx_columns]))+'\n')
                                        gtx_record_count += 1

        #close last part
        if this_part_fh:
            this_part_fh.close()
            with open(this_part_done, 'w'):
                pass

       #finished reading this vcf
        if create_vex or create_gtx or create_parts:
            print 'input vcf file {} had {} lines, {} records.'.format(vcf, vcf_line_count, vcf_record_count)
        if create_vex:
            vex_out.close()
            with open(vex_done, 'w'):
                pass
            print 'Output vex.vcf file {} had {} records.'.format(vex, vex_record_count)
            
        if  create_gtx:
            gtx_out.close()
            with open(gtx_done, 'w'):
                pass
            print 'Output gtx file {} had {} records.'.format(gtx, gtx_record_count)
    
        if  create_parts:
            with open(parts_done, 'w'):
                pass
            print '{} part files done.'.format(len(this_vcf_all_part_files))
    
        #sql script jobs
        vex_microsoft_sql_file = vex+'.microsoft.sql'
        vex_mysql_file = vex+'.mysql'
        vex_microsoft_done_file = my.done_file(vex_microsoft_sql_file)
        vex_mysql_done_file = my.done_file(vex_mysql_file)

        gtx_microsoft_sql_file = gtx+'.microsoft.sql'
        gtx_mysql_file = gtx+'.mysql'
        gtx_microsoft_done_file = my.done_file(gtx_microsoft_sql_file)
        gtx_mysql_done_file = my.done_file(gtx_mysql_file)

        sql_cmds = []
        sql_columns_py = os.path.join(os.path.dirname(__file__), 'sql_columns.py')
        job_name = 'vax_sql_cols'

        #commands to make sql scripts
        if (my.file_exists(vex) and my.file_exists(vex_done) and
            (not my.file_exists(vex_microsoft_sql_file)
            or not my.file_exists(vex_microsoft_done_file)
            or not my.file_exists(vex_mysql_file)
            or not my.file_exists(vex_mysql_done_file))
            ):
            table = os.path.basename(vex)
            sql_cmds.append('{} {} --input {} --database vax_data --schema vax_data --table {} --clustered_index CHROM POS'.format(
                python, sql_columns_py, vex, table,))
            job_name += '_'+table
        
        if  (my.file_exists(gtx) and my.file_exists(gtx_done) and
            (not my.file_exists(gtx_microsoft_sql_file)
            or not my.file_exists(gtx_microsoft_done_file)
            or not my.file_exists(gtx_mysql_file)
            or not my.file_exists(gtx_mysql_done_file))
            ):
            table = os.path.basename(gtx)
            sql_cmds.append('{} {} --input {} --database vax_data --schema vax_data --table {} --clustered_index CHROM POS SAMPLE --indexes SAMPLE,CHROM,POS'.format(
                python, sql_columns_py, gtx, table,))
            job_name += '_'+table
        
        #job to run sql script commands
        if sql_cmds:
            job = my.run_job(sql_cmds, job_name, job_dir)
            sql_jobid = job.jobId
            job_ids.append(job.jobId)
            print 'Submitted job {} to create SQL scripts.'.format(job.jobId,)
            
        #schedule jobs to annotate each of the split vex.vcf parts with variant effect predictor
        if create_vax:
            vex_part_files = my.unglob(os.path.join(parts_dir, '*.'+os.path.basename(vex)))
            this_vcf_part_files_to_vep = []
            this_vcf_vax_files_to_merge = []
            this_vcf_all_vax_parts = []
            for f in vex_part_files:
                if not my.file_exists(my.done_file(f)):
                    raise Vcf2VaxError('The vex part file {} may have problems; there is no corresponding .done file'.format(f))
                this_vcf_part_files_to_vep.append(f)
                    
            #TODO: capture job ids of currently running vax jobs and add to hold_jids
            cmds = []
            all_vcfs_all_vax.append(vax)
            this_vcf_vax_parts_to_vep = []
            this_vcf_job_ids = []
            this_job_number = 0
            #run VEP on vcf parts that have not already run
            #limit number of simultaneous jobs to keep the mysql server happy
            for (part_number, vcf_part) in enumerate(this_vcf_part_files_to_vep, start=1):
                vax_part = vcf_part+'.vax'
                
                #the .done file must be written by the variant effect predictor
                #this requires a patch to VEP
                ##MY:
                ##write file to show job finished in the form '.filename.done'
                #use File::Spec;
                #my($volume,$directories,$file) = File::Spec->splitpath( $config->{output_file} );
                #my $done_file = File::Spec->catpath( $volume,$directories,'.'.$file.'.done' );
                #open(DONE_FILE,'>', $done_file) or die "Can't create $done_file";
                #close(DONE_FILE);
                #debug("Finished!") unless defined $config->{quiet};
                vax_part_done = my.done_file(vax_part)
                
                this_vcf_all_vax_parts.append(vax_part)
                all_vcfs_all_vax_parts.append(vax_part)
                if my.file_exists(vax_part) and my.file_exists(vax_part_done):
                    print('{} exists. Skipping VEP.'.format(vax_part))
                else:
                    cmds.append('{} {} -input_file {} -output_file {} --config {} {}'.format(perl, vep, vcf_part, vax_part, vep_ini, vep_params))
                    this_vcf_vax_parts_to_vep.append(vax_part)
                    all_vcfs_vax_parts_to_vep.append(vax_part)
                    if len(cmds) >= max_simultaneous_jobs or part_number >= len(this_vcf_part_files_to_vep):
                        this_job_number += 1
                        job_name = 'vax_{}.{}'.format(this_job_number, vcf_base)
                        job = my.run_job(cmds, job_name, job_dir, hold_jid=all_vax_job_ids)
                        vax_jobid = job.jobId
                        this_vcf_job_ids.append(job.jobId)
                        all_vax_job_ids.append(job.jobId)
                        cmds = []
                        print 'Submitted job {}, batch {} to annotate {}.'.format(vax_jobid, this_job_number, vcf_base)
            
            #merge and post-process annotated vax files
            cmd= '{} {} -i {} -o {}'.format(
                python, vax_merge_post_process.__file__, ' '.join(this_vcf_all_vax_parts), vax)
            job_name = 'vax_merge_'+vcf_base
            job = my.run_job(cmd, job_name, job_dir, hold_jid=this_vcf_job_ids, email=email)
            merge_vax_jobid = job.jobId
            this_vcf_job_ids.append(job.jobId)
            print 'Submitted job {} to merge and post-process annotations into {}.'.format(merge_vax_jobid, vax)
            
            job_ids += this_vcf_job_ids

    #remove temporary files
    if deletetempfiles:
        cmd = 'rm -rf {}'.format(' '.join(parts_dirs))
        job_name = 'vax_rm_'+vcf_base
        job = my.run_job(cmd, job_name, job_dir, hold_jid=job_ids)
        rm_jobid = job.jobId
        this_vcf_job_ids.append(job.jobId)
        job_ids.append(job.jobId)
        print 'Submitted job {} to remove temporary files from {}.'.format(rm_jobid, ' '.join(parts_dir))
    else:
        print """
        When all jobs are done, please manually remove temporary files from {}.
        In future you can use the --deletetempfiles option to auto-delete the parts_* directories.""".format(' '.join(parts_dirs))

    print """
    When all jobs are done, check STDERR output in {} to validate the run.""".format(' '.join(job_dirs))
    return job_ids


def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'annotate one or more vcf files with Ensembl Variant Effect Predictor + extras',
        epilog = '{} version {} {}'.format(name, version, copyright))
    #input VCF parameter
    parser.add_argument('--vcfs', '-v', '-i', nargs='+', required=True,
        help='input vcf file(s); output will be <vcf>.vax')
    #vcf2vex parameters
    parser.add_argument('--max_part_records', type=int, default=1000,
        help='maximum number of records per vex chromosome part file (keep small to avoid locking out other cluster users; default:1000; range 100-5000)')
    parser.add_argument('--max_simultaneous_jobs', type=int, default=20,
        help='maximum number of vep jobs to run simultaneously (used to choke mysql server load, default:20; range 1-40)')
    parser.add_argument('--no_gtx', action='store_true', default=False,
        help='do not create gtx files (default:False)')
    parser.add_argument('--no_vex', action='store_true', default=False,
        help='do not create vex files (default:False)')
    parser.add_argument('--no_vax', action='store_true', default=False,
        help='do not run variant effect predictor, just create gtx and or vex files files (default:False)')
    #variant effect predictor parameters
    parser.add_argument('--vep',
        help='Ensembl variant effect predictor (default: config->vax.vep_ini)')
    parser.add_argument('--vep_plugins',
        help='path to VEP_plugins directory (default: config->vax.vep_plugins)')
    parser.add_argument('--vep_ini',
        help='path to vep.ini, which contains VEP parameters and plugins to run (default: config->vax.vep_ini)')
    parser.add_argument('--vep_params',
        help='a single quoted string with VEP parameters to override vep_ini (default: None)')
    #output parameters
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
        max_part_records=args.max_part_records, max_simultaneous_jobs=args.max_simultaneous_jobs, no_gtx=args.no_gtx, no_vex=args.no_vex, no_vax=args.no_vax,
        vep=args.vep, vep_plugins=args.vep_plugins, vep_ini=args.vep_ini, vep_params=args.vep_params,
        deletetempfiles=args.deletetempfiles,
        python=args.python, perl=args.perl,
        email=args.email)


if __name__ == "__main__": sys.exit(main())
