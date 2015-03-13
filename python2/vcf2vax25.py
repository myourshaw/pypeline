#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import gzip
import my
import job
import merge_header_files
import vax2flax

class Vcf2VaxError(Exception): pass

#-d /home/myourshaw/lab/pypeline/vax_test -v /home/myourshaw/lab/pypeline/vax_test/gmd_test.vcf --email myourshaw@ucla.edu
#python /home/myourshaw/lab/pypeline/python2/vcf2vax.py --email myourshaw@ucla.edu -d /scratch1/tmp/myourshaw/hlee -v /scratch1/tmp/myourshaw/hlee/vcfs/hlee55gmd28bipolar20ocdtic18.analysis_ready.vcf
#-d /home/myourshaw/lab/pypeline/vax_test/foo -v /home/myourshaw/lab/pypeline/vax_test/kitchen_sink_GT_test.vcf --email myourshaw@ucla.edu

def run(top_dir, vcfs, records_per_job=5000, parent_job_dir=None, vax=None,
        vep_plugins=None, vep_ini=None, vep_params='',
        revised_vax=None, no_revised_vax=False, vax_consequence_threshhold=None,
        flax_consequence_threshhold=None, flax_suppress_uncalled_genotypes=False,
        flax_suppress_homozygous_reference_genotypes=False,
        flax_include_all_genotype_columns=False, vax_uncompressed=False,
        revised_vax_uncompressed=False, flax_uncompressed=False,
        no_flax=False,
        python=None, perl=None,
        email=None, deletetempfiles=False):
        
    dirs = my.create_default_directory_structure(top_dir)
    if not parent_job_dir or not os.path.isdir(parent_job_dir):
        parent_job_dir = dirs['jobs']
    config = my.get_config()
        
    my_name = 'vax'
    my_job_dir =  my.unique_dir('vax_'+my.localtime_squish(), dirs['jobs'])
        
    #retun value: a list of job ids that can be used by downstream hold_jid
    job_ids = []
    
    vcfs_glob = vcfs
    vcfs = my.unglob(vcfs)
    if not vcfs:
        raise Vcf2VaxError("no files in {}".format(vcfs_glob))
    if not vax:
        vax = config.get('vax','vax')
    if not my.file_exists(vax):
        raise Vcf2VaxError("can't find vax perl script {}".format(vax))
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
        
    tmp_files_to_remove = []
    for vcf in vcfs:
        
        vcf_job_ids = []
        #split each vcf into files having no more than records_per_job records
        print 'splitting '+vcf
        output_dir = os.path.dirname(vcf)
        vcf_base = os.path.basename(vcf)
        my.makedir(output_dir)
        vax_out = os.path.join(output_dir, vcf_base+'.vax') + ('.gz' if not vax_uncompressed else '')
        tempdir = os.path.join(output_dir,'tmp_'+vcf_base)
        my.makedir(tempdir)
        part_base = '{}.part.'.format(os.path.join(tempdir, vcf_base))
        job_name = 'split_'+vcf_base
        job_dir = tempdir
        vcf_parts = my.split_file(vcf, records_per_job, suffix_length=8, prefix=part_base, repeat_header=True, compress=False)
        tmp_files_to_remove += vcf_parts
        
        #annotate each of the split vcfs
        cmds = []
        vax_parts = []
        for vcf_part in vcf_parts:
            vax_part = vcf_part+'.vax'
            vax_parts += [vax_part]
            cmds.append('{} {} -input_file {} -output_file {} --config {} {}'.format(perl, vax, vcf_part, vax_part, vep_ini, vep_params))
        tmp_files_to_remove += vax_parts
        job_name = 'vax_'+vcf_base
        job = my.run_job(cmds, job_name, job_dir)
        vax_jobid = job.jobId
        vcf_job_ids += [job.jobId]
        print 'submitted job {} to annotate {}'.format(vax_jobid, vcf)
        
        #merge partial vax files
        cmd= '{} {} -i {} -o {} --must_have_data --all_data_records_must_have_all_columns --unique_key_data_columns 0 1 3 4'.format(
            python, merge_header_files.__file__, ' '.join(vax_parts), vax_out, )
        if not vax_uncompressed:
            cmd += ' --compress'
        job_name = 'vax_merge_'+vcf_base
        job = my.run_job(cmd, job_name, job_dir, hold_jid=vax_jobid, email=email)
        merge_vax_jobid = job.jobId
        vcf_job_ids += [job.jobId]
        print 'submitted job {} to merge annotations into {}'.format(merge_vax_jobid, vax_out)
        
        #vax2flax
        if not no_flax:
            flax_params = []
            if revised_vax:
                flax_params.append('--revised_vax {}'.format(revised_vax))
            if no_revised_vax:
                flax_params.append('--no_revised_vax')
            if vax_consequence_threshhold:
                flax_params.append('--vax_consequence_threshhold {}'.format(vax_consequence_threshhold))
            if flax_consequence_threshhold:
                flax_params.append('--flax_consequence_threshhold {}'.format(flax_consequence_threshhold))
            if flax_suppress_uncalled_genotypes:
                flax_params.append('--flax_suppress_uncalled_genotypes')
            if flax_suppress_homozygous_reference_genotypes:
                flax_params.append('--flax_suppress_homozygous_reference_genotypes')
            if flax_include_all_genotype_columns:
                flax_params.append('--flax_include_all_genotype_columns')
            if revised_vax_uncompressed:
                flax_params.append('--revised_vax_uncompressed')
            if flax_uncompressed:
                flax_params.append('--flax_uncompressed')
            cmd= '{} {} -i {} {}'.format(python, vax2flax.__file__, vax_out, ' '.join(flax_params))
            job_name = 'vax2flax_'+vcf_base
            job = my.run_job(cmd, job_name, job_dir, hold_jid=merge_vax_jobid, email=email)
            vax2flax_jobid = job.jobId
            vcf_job_ids += [job.jobId]
            print 'submitted job {} to create flax and revised vax files'.format(vax2flax_jobid)

        #remove temporary files
        if deletetempfiles:
            cmd = 'rm -rf {}'.format(' '.join(vax_parts))
            job_name = 'vax_rm_'+vcf_base
            job = my.run_job(cmd, job_name, job_dir, hold_jid=vcf_job_ids)
            rm_jobid = job.jobId
            vcf_job_ids += [job.jobId]
            print 'submitted job {} to remove temporary files from dir {}'.format(rm_jobid, tempdir)
        else:
            print 'when all jobs are done, please manually remove temporary files from {}\nin future use the --deletetempfiles option to auto-delete the tmp directory'.format(tempdir)
        job_ids += vcf_job_ids

    print 'when all jobs are done, check STDERR output in the tmp_* directories to validate the run' 
    return job_ids

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'annotate one or more vcf files with Ensembl Variant Effect Predictor + extras',
        epilog = 'pypeline.vcf2vax version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--vcfs', '-v', '-i', nargs='+', required=True,
        help='input vcf file(s); output will be <vcf>.vax')
    parser.add_argument('--records_per_job', type=int, default=5000,
        help='number of vcf records to process in each job (default: 5000)')
    parser.add_argument('--vax',
        help='path to Variant Annotator Xtra perl executable (default: config->vax.vax)')
    parser.add_argument('--vep_plugins',
        help='path to VEP_plugins directory (default: config->vax.vep_plugins)')
    parser.add_argument('--vep_ini',
        help='path to vep.ini, which contains VEP parameters and plugins to run (default: config->vax.vep_ini or config->vax.vep_no_gt_ini)')
    parser.add_argument('--vep_params',
        help='a single quoted string with VEP parameters to override vep_ini (default: None)')
    parser.add_argument('--revised_vax',
        help='revised vax output with "-" replaced by "" and optionally filtered by consequence (default:<input>.revised_vax.gz)')
    parser.add_argument('--no_revised_vax', action='store_true', default=False,
        help='do not create revised vax output with "-" replaced by "" and optionally filtered by consequence (default: False)')
    parser.add_argument('--vax_consequence_threshhold', type=int, default=None,
        help='only output revised vax records where Consequence_rank <= vax_consequence_threshhold (default: None)')
    parser.add_argument('--flax_consequence_threshhold', type=int, default=None, #14, #Incomplete terminal codon variant in Ensembl 68
        help='only output flax records where Consequence_rank <= flax_consequence_threshhold (default: 9 = non_synonymous_codon)')
    parser.add_argument('--flax_suppress_uncalled_genotypes', action='store_true', default=False,
        help='do not output uncalled (e.g., ./.) genotypes flax records (default: False)')
    parser.add_argument('--flax_suppress_homozygous_reference_genotypes', action='store_true', default=False,
        help='do not output homozygous_reference flax records (e.g., 0/0) genotypes (default: False)')
    parser.add_argument('--flax_include_all_genotype_columns', action='store_true', default=False,
        help='include the genotype columns for all samples in flax records (default: False)')
    parser.add_argument('--vax_uncompressed', action='store_true', default=False,
        help='do not gzip vax output (default: False)')
    parser.add_argument('--revised_vax_uncompressed', action='store_true', default=False,
        help='do not gzip revised vax output (default: False)')
    parser.add_argument('--flax_uncompressed', action='store_true', default=False,
        help='do not gzip flax output (default: False)')
    parser.add_argument('--no_flax', action='store_true', default=False,
        help='do not run vax2flax (default: False)')
    parser.add_argument('--python',
        help='path to python executable (default: config->python)')
    parser.add_argument('--perl',
        help='path to perl executable (default: config->perl)')
    parser.add_argument('--genotype_file_output', action='store_true', default=False,
        help='create potentially large <vax file>.gt.txt output with records repeated for each sample ID (default: False)')
    parser.add_argument('--deletetempfiles', action='store_true', default=False,
        help='delete temporary files (default: False)')
    args = parser.parse_args()

    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    job_ids = run(top_dir=args.top_dir, vcfs=args.vcfs, records_per_job=args.records_per_job,
        parent_job_dir=args.parent_job_dir, vax=args.vax, vep_plugins=args.vep_plugins,
        vep_ini=args.vep_ini, vep_params=args.vep_params,
        revised_vax=args.revised_vax, no_revised_vax=args.no_revised_vax, vax_consequence_threshhold=args.vax_consequence_threshhold,
        flax_consequence_threshhold=args.flax_consequence_threshhold, flax_suppress_uncalled_genotypes=args.flax_suppress_uncalled_genotypes,
        flax_suppress_homozygous_reference_genotypes=args.flax_suppress_homozygous_reference_genotypes,
        flax_include_all_genotype_columns=args.flax_include_all_genotype_columns, vax_uncompressed=args.vax_uncompressed,
        revised_vax_uncompressed=args.revised_vax_uncompressed, flax_uncompressed=args.flax_uncompressed,
        no_flax=args.no_flax,
        python=args.python, perl=args.perl,
        email=args.email, deletetempfiles=args.deletetempfiles)


if __name__ == "__main__": sys.exit(main())
