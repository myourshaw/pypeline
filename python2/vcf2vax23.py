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
#--keeptempfiles -d /scratch0/tmp/myourshaw/gmd28mms -v /scratch0/tmp/myourshaw/gmd28mms/vcfs/gmd28mms.analysis_ready.vcf
#python /home/myourshaw/lab/pypeline/python/vcf2vax.py --email myourshaw@ucla.edu --keeptempfiles -d /scratch0/tmp/myourshaw/gmd28mms -v /scratch0/tmp/myourshaw/gmd28mms/vcfs/gmd28mms.analysis_ready.vcf
#-d /home/myourshaw/lab/pypeline/vax_test/ -v /home/myourshaw/lab/pypeline/vax_test/hgmd_test.vcf


def run(top_dir, vcfs, records_per_job=5000, parent_job_dir=None, vax=None,
        vep_plugins=None, vep_ini=None, vep_params='', python=None, perl=None,
        email=None, genotype_file_output=False, keeptempfiles=False):
    
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
        vax = config.get('vax','vax23')
    if not my.file_exists(vax):
        raise Vcf2VaxError("can't find vax perl script {}".format(vax))
    if not vep_ini:
        vep_ini = config.get('vax','vep_large_gt_file') if genotype_file_output else config.get('vax','vep_ini_23')
    if not my.file_exists(vep_ini):
        raise Vcf2VaxError("can't find vep.ini {}".format(vep_ini))
    if not vep_plugins:
        vep_plugins = config.get('vax','vep_plugins_23')
    if not vep_params:
        vep_params = ''
    if not python:
        python = config.get('DEFAULT','python')
    if not perl:
        perl = config.get('DEFAULT','perl')

    all_lines_count = 0
    for vcf in vcfs:
        
        #split each vcf into files having no more than records_per_job records
        print 'splitting '+vcf
        output_dir = os.path.dirname(vcf)
        vcf_base = os.path.basename(vcf)
        my.makedir(output_dir)
        vax_out = os.path.join(output_dir, vcf_base+'.vax')
        gt_out = vax_out+'.gt.txt'
        tempdir = my.unique_dir(os.path.join(output_dir,'tmp'))
        my.makedir(tempdir)
        part_base = '{}.part.'.format(os.path.join(tempdir, vcf_base))
        job_name = 'split_'+vcf_base
        job_dir = tempdir
        vcf_parts = my.split_file(vcf, records_per_job, 8, part_base, True)
        
        #annotate each of the split vcfs
        print 'submitting jobs to annotate '+vcf
        cmds = []
        vax_parts = []
        gt_parts = []
        for vcf_part in vcf_parts:
            vax_part = vcf_part+'.vax'
            vax_parts += [vax_part]
            gt_part = vax_part+'.gt.txt'
            gt_parts += [gt_part]
            cmds.append('{} {} \
-input_file {} \
-output_file {} \
--config {} \
{} \
'.format(perl, vax, vcf_part, vax_part, vep_ini, vep_params))
        job_name = 'vax_'+vcf_base
        job = my.run_job(cmds, job_name, job_dir)
        vax_jobid = job.jobId
        job_ids += [job.jobId]

        #merge partial vax files
        print 'submitting job to merge annotations into '+vax_out
        cmd= '{} {} -i {} -o {}'.format(
            python, merge_header_files.__file__, ' '.join(vax_parts), vax_out)
        job_name = 'merge_vax_'+vcf_base
        job = my.run_job(cmd, job_name, job_dir, hold_jid=vax_jobid, email=email)
        merge_vax_jobid = job.jobId
        job_ids += [job.jobId]

        #merge partial vax.gt.txt files
        if genotype_file_output:
            print 'submitting job to merge genotype annotations into '+gt_out
            cmd= '{} {} -i {} -o {}'.format(
                python, merge_header_files.__file__, ' '.join(gt_parts), gt_out)
            job_name = 'merge_gt_'+vcf_base
            job = my.run_job(cmd, job_name, job_dir, hold_jid=vax_jobid)
            merge_gt_jobid = job.jobId
            job_ids += [job.jobId]

        #remove temporary files
        if not keeptempfiles:
            print 'submitting job to remove temporary dir '+tempdir
            cmd = 'rm -rf {}'.format(tempdir)
            job_name = 'rm_'+vcf_base
            job = my.run_job(cmd, job_name, job_dir, hold_jid=job_ids)
            rm_jobid = job.jobId
            job_ids += [job.jobId]
        else:
            print 'when all jobs are done, please manually remove temporary dir '+tempdir

    return job_ids

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'annotate one or more vcf files with Ensembl Variant Effect Predictor + extras',
        epilog = 'pypeline.vcf2vax version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--vcfs', '-v', nargs='+', required=True,
        help='input vcf file(s); output will be <vcf>.vax')
    parser.add_argument('--records_per_job', type=int, default=5000,
        help='number of vcf records to process in each job (default: 5000)')
    parser.add_argument('--vax',
        help='path to Variant Annotator Xtra perl executable (default: config->vax.vax23)')
    parser.add_argument('--vep_plugins',
        help='path to VEP_plugins directory (default: config->vax.vep_plugins)')
    parser.add_argument('--vep_ini',
        help='path to vep.ini, which contains VEP parameters and plugins to run (default: config->vax.vep_ini or config->vax.vep_no_gt_ini)')
    parser.add_argument('--vep_params',
        help='a single quoted string with VEP parameters to override vep_ini (default: None)')
    parser.add_argument('--python',
        help='path to python executable (default: config->python)')
    parser.add_argument('--perl',
        help='path to perl executable (default: config->perl)')
    parser.add_argument('--genotype_file_output', action='store_true', default=False,
        help='create potentially large <vax file>.gt.txt output with records repeated for each sample ID (default: False)')
    parser.add_argument('--keeptempfiles', action='store_true', default=False,
        help='keep temporary files (default: False)')
    args = parser.parse_args()

    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    job_ids = run(top_dir=args.top_dir, vcfs=args.vcfs, records_per_job=args.records_per_job,
        parent_job_dir=args.parent_job_dir, vax=args.vax, vep_plugins=args.vep_plugins,
        vep_ini=args.vep_ini, vep_params=args.vep_params, python=args.python, perl=args.perl,
        email=args.email, genotype_file_output=args.genotype_file_output, keeptempfiles=args.keeptempfiles)


if __name__ == "__main__": sys.exit(main())
