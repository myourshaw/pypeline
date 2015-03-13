#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import my
import job
import merge_header_files

class Vcf2VaxError(Exception): pass

#--keeptempfiles -d /scratch0/tmp/myourshaw/gmd/ -v /scratch0/tmp/myourshaw/gmd/vcfs/gmd32.analysis_ready.vcf
#--keeptempfiles -d /scratch0/tmp/myourshaw/gmd/ -v /scratch0/tmp/myourshaw/gmd/vcfs/gmd32.analysis_ready.vcf
#--keeptempfiles -d /scratch0/tmp/myourshaw/gmd/ -v /scratch0/tmp/myourshaw/gmd/vcfs/PCSK1_variants.vcf
#--records_per_job 25000 --keeptempfiles -d /scratch0/tmp/myourshaw/gmdmms/ -v /scratch0/tmp/myourshaw/gmdmms/vcfs/gmd35mms.analysis_ready.vcf
#python /home/myourshaw/lab/pypeline/python/vcf2vax.py -d /scratch0/tmp/myourshaw/gmd28mms -v /scratch0/tmp/myourshaw/gmd28mms/vcfs/gmd28mms.analysis_ready.vcf
#--keeptempfiles -d /scratch0/tmp/myourshaw/hlee -v /scratch0/tmp/myourshaw/hlee/UCLA_ClinSeq_11042011.snps.recalibrated.vcf

def run(top_dir, vcfs, records_per_job=1000, parent_job_dir=None, vax=None,
        python=None, perl=None,
        host=None, port=None, user=None, password=None, species='human',
        vw_database=None, vw_host=None, vw_password=None, vw_platform=None,
        vw_port=None, vw_user=None,
        email=None, keeptempfiles=False):
    
    dirs = my.create_default_directory_structure(top_dir)
    if not parent_job_dir or not os.path.isdir(parent_job_dir):
        parent_job_dir = dirs['jobs']
    config = my.get_config()

    my_name = 'vax'
    my_job_dir =  my.unique_dir('vax_'+my.localtime_squish(), dirs['jobs'])
    
    #retun value: a list of job ids that can be used by downstream hold_jid
    job_ids = []

    vcfs = my.unglob(vcfs)
    if not vax:
        vax = config.get('vax','script')
    if not python:
        python = config.get('DEFAULT','python')
    if not perl:
        perl = config.get('DEFAULT','perl')
    if not my.file_exists(vax):
        raise Vcf2VaxError("can't find vax perl script {}".format(vax))
    if not host:
        host = config.get('vax','host')
    if not port:
        port = config.get('vax','port')
    if not user:
        user = config.get('vax','user')
    if not password:
        password = config.get('vax','password')
    if not vw_database:
        vw_database = config.get('vax','vw_database')
    if not vw_host:
        vw_host = config.get('vax','vw_host')
    if not vw_password:
        vw_password = config.get('vax','vw_password')
    if not vw_platform:
        vw_platform = config.get('vax','vw_platform')
    if not vw_port:
        vw_port = config.get('vax','vw_port')
    if not vw_user:
        vw_user = config.get('vax','vw_user')

    all_lines_count = 0
    for vcf in vcfs:
        
        #split each vcf into files having no more than records_per_job records
        output_dir = os.path.dirname(vcf)
        vcf_base = os.path.basename(vcf)
        my.makedir(output_dir)
        vax_out = os.path.join(output_dir, vcf_base+'.vax')
        tempdir = my.unique_dir(os.path.join(output_dir,'tmp'))
        my.makedir(tempdir)
        part_base = '{}.part.'.format(os.path.join(tempdir, vcf_base))
        job_name = 'split_'+vcf_base
        job_dir = tempdir
        vcf_parts = my.split_file(vcf, records_per_job, 8, part_base, True)
        
        #annotate each of the split vcfs
        cmds = []
        vax_parts = []
        for vcf_part in vcf_parts:
            vax_part = vcf_part+'.vax'
            vax_parts += [vax_part]
            cmds.append('{} {} \
--vcf_out \
--buffer_size 1000 \
-input_file {} \
--format vcf \
-output_file {} \
--force_overwrite \
--terms SO \
--sift=b \
--polyphen=b \
--condel=b \
--regulatory \
--hgnc \
--hgvs \
--protein \
--gene \
--check_ref \
--host {} \
--port={} \
--user {} \
--password {} \
--species {} \
--no_progress \
--vw_database {} \
--vw_host {} \
--vw_password {} \
--vw_platform {} \
--vw_port {} \
--vw_user {} \
'.format(perl, vax, vcf_part, vax_part, host, port, user, password, species,
        vw_database, vw_host, vw_password, vw_platform, vw_port, vw_user))
        job_name = 'vax_'+vcf_base
        job = my.run_job(cmds, job_name, job_dir)
        vax_jobid = job.jobId
        job_ids += [job.jobId]

        #merge partial vax files
        cmd= '{} {} -i {} -o {}'.format(
            python, merge_header_files.__file__, ' '.join(vax_parts), vax_out)
        job_name = 'merge_'+vcf_base
        job = my.run_job(cmd, job_name, job_dir, hold_jid=vax_jobid)
        merge_jobid = job.jobId
        job_ids += [job.jobId]

        #remove temporary files
        if not keeptempfiles:
            cmd = 'rm -rf {}'.format(tempdir)
            job_name = 'rm_'+vcf_base
            job = my.run_job(cmd, job_name, job_dir, hold_jid=merge_jobid)
            rm_jobid = job.jobId
            job_ids += [job.jobId]
            
    return job_ids

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'annotate one or more vcf files with Ensembl Variant Effect Predictor + extras',
        epilog = 'pypeline.vcf2vax version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--vcfs', '-v', nargs='+',
        help='input vcf file(s); output will be <vcf>.vax')
    parser.add_argument('--records_per_job', type=int, default=1000,
        help='number of vcf records to process in each job (default: 1000)')
    parser.add_argument('--vax',
        help='path to Variant Annotator Xtra perl executable (default: config->vax.script)')
    parser.add_argument('--python',
        help='path to python executable (default: config->python)')
    parser.add_argument('--perl',
        help='path to perl executable (default: config->perl)')
    parser.add_argument('--host',
        help='vax database host (default: config->vax.host)')
    parser.add_argument('--port',
        help='vax database port (default: config->vax.port)')
    parser.add_argument('--user',
        help='vax database user (default: config->vax.user)')
    parser.add_argument('--password',
        help='vax database password (default: config->vax.password)')
    parser.add_argument('--species', default='human',
        help='vax database species (default: human)')
    parser.add_argument('--vw_database',
        help='variant warehouse database (default: config->vax.vw_database)')
    parser.add_argument('--vw_host',
        help='variant warehouse database host (default: config->vax.vw_host)')
    parser.add_argument('--vw_password',
        help='variant warehouse database password (default: config->vax.vw_password)')
    parser.add_argument('--vw_platform',
        help='variant warehouse database platform (default: config->vax.vw_platform)')
    parser.add_argument('--vw_port',
        help='variant warehouse database port (default: config->vax.vw_port)')
    parser.add_argument('--vw_user',
        help='variant warehouse database user (default: config->vax.vw_user)')
    parser.add_argument('--keeptempfiles', action='store_true', default=False,
        help='keep temporary files (default: False)')
    args = parser.parse_args()

    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    job_ids = run(top_dir=args.top_dir, vcfs=args.vcfs, records_per_job=args.records_per_job,
        parent_job_dir=args.parent_job_dir, vax=args.vax, python=args.python, perl=args.perl, host=args.host, port=args.port,
        user=args.user, password=args.password, species=args.species, vw_database=args.vw_database,
        vw_host=args.vw_host, vw_password=args.vw_password, vw_platform=args.vw_platform,
        vw_port=args.vw_port, vw_user=args.vw_user, email=args.email, keeptempfiles=args.keeptempfiles)
    
    #return job_ids


if __name__ == "__main__": sys.exit(main())
