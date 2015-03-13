#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import ConfigParser #configparser in python 3
import tempfile
from time import localtime, strftime
import my
import job
import validate_brlmm_file
import verify_OK_file
import validate_vcf

class Bam2HomozygosityBlocksError(Exception): pass

name = 'pypeline.bam2homozygosity_blocks'
version = 1.0
copyright = '©2011-2013 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()

#-d /scratch1/tmp/myourshaw/mmjj_20130514 -i /scratch1/tmp/myourshaw/mmjj_20130514/bams/sample_bams/*.sample.bam

def run(input, linkdatagen_dir, plink_dir, job_dir,
        my_name, hold_jids, 
        python, perl, reference, samtools,
        bcftools, annotHapMap2L_b37, annotHapMap2,
        vcf2linkdatagen, linkdatagen_mps, plink,
        overwrite_all_output_files, overwrite_pileups):

    print 'Submitting jobs to calculate homozygosity blocks with linkdatagen and plink ...'
    
    hold_linkdatagen_jids = []
    annotHapMap2_dir = os.path.dirname(annotHapMap2)
    annotHapMap2_file = os.path.basename(annotHapMap2)
    
    input_files = my.unglob(input)
    bams = [i for i in input_files if i]
    if not bams:
        raise Bam2HomozygosityBlocksError('no bam files to process')
    #file to tell linkdatagen to use first (only) genotype in bcf_vcf
    linkdatagen_whichsamplesfile = os.path.join(linkdatagen_dir, '1.ws')
    with open(linkdatagen_whichsamplesfile, 'w') as ws:
        ws.write('1')
    for bam in bams:
        bcf = os.path.join(linkdatagen_dir, os.path.basename(bam)+'.bcf')
        bcf_vcf = bcf+'.vcf'
        bcf_vcf_validate = bcf_vcf+'.validate'
        brlmm = bcf_vcf+'.brlmm'
        brlmm_validate = brlmm+'.validate'
        sample_prefix = my.r_strip(my.bam_strip(bam), '.sample')
        sample_name = os.path.basename(sample_prefix)
        linkdatagen_mps_stdout = os.path.join(linkdatagen_dir, sample_name+'.Ped_HapMap2_pl.out')
        sample_plink_dir = os.path.join(plink_dir, sample_name)
        plink_makebed_stdout = os.path.join(plink_dir, sample_name+'.plink_makebed.out')
        plink_homozyg_stdout = os.path.join(plink_dir, sample_name+'.plink_homozyg.out')
        #fake ped file for linkdatagen (all female to get X)
        linkdatagen_pedfile = os.path.join(linkdatagen_dir, sample_name+'.ped')
        with open(linkdatagen_pedfile,'w') as ped:
            ped.write('0001\t001\t0\t0\t2\t2')
        if overwrite_all_output_files or overwrite_pileups or not my.check_files(files_exist=[bcf,bcf_vcf,brlmm], files_contain=[(bcf_vcf_validate,'OK'),(brlmm_validate,'OK')], files_not_before=[(bcf,bam),(bcf_vcf,bcf),(brlmm,bcf_vcf)]):
            
            #mpileup
            job_name = my_name+'_mpileup_'+os.path.basename(bcf)
            cmd = '{} mpileup -C50 -d10000 -q13 -gf {} -l {} {} > {};'.format(samtools, reference, annotHapMap2L_b37, bam, bcf)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_jids)
            mpileup_jobid = job.jobId
            hold_linkdatagen_jids += [job.jobId]
            
            #bcftools
            job_name = my_name+'_bcftools_'+os.path.basename(bcf_vcf)
            cmd = '{} view -cg -t0.5 {} > {};'.format(bcftools, bcf, bcf_vcf)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=mpileup_jobid)
            bcftools_jobid = job.jobId
            hold_linkdatagen_jids += [job.jobId]
            
            #validate bcf_vcf
            job_name = my_name+'_ValidateBcfVcfFile_'+os.path.basename(bcf_vcf)
            cmd = '{} {} -i {} -o {}'.format(python, validate_vcf.__file__, bcf_vcf, bcf_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=bcftools_jobid)
            validatesamplebcfvcffile_jobid = job.jobId
            hold_linkdatagen_jids += [job.jobId]
            
            #verify bcf_vcf validation
            job_name = my_name+'_VerifyValidateSampleBcfVcfFile_'+os.path.basename(bcf_vcf)
            cmd = '{} {} -v {}'.format(python, verify_OK_file.__file__, brlmm_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=validatesamplebcfvcffile_jobid)
            verifyvalidatesamplebcfvcffile_jobid = job.jobId
            hold_linkdatagen_jids += [job.jobId]
            
            #vcf2linkdatagen
            job_name = my_name+'_vcf2linkdatagen_'+os.path.basename(brlmm)
            cmd = '{} {} -annotfile {} -mindepth 5 -missingness 0 {} > {};'.format(perl, vcf2linkdatagen, annotHapMap2, bcf_vcf, brlmm)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=bcftools_jobid)
            vcf2linkdatagen_jobid = job.jobId
            hold_linkdatagen_jids += [job.jobId]
            
            #validate brlmm
            job_name = my_name+'_ValidateBrlmmFile_'+os.path.basename(brlmm)
            cmd = '{} {} -i {} -o {}'.format(python, validate_brlmm_file.__file__, brlmm, brlmm_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=vcf2linkdatagen_jobid)
            validatesamplebrlmmfile_jobid = job.jobId
            hold_linkdatagen_jids += [job.jobId]
            
            #verify brlmm validation
            job_name = my_name+'_VerifyValidateSampleBrlmmFile_'+os.path.basename(brlmm)
            cmd = '{} {} -v {}'.format(python, verify_OK_file.__file__, brlmm_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=validatesamplebrlmmfile_jobid)
            verifyvalidatesamplebrlmmfile_jobid = job.jobId
            hold_linkdatagen_jids += [job.jobId]
            
            #linkdatagen_mps
            job_name = my_name+'_linkdatagen_mps_'+os.path.basename(brlmm)
            cmd = '{} {} -pedfile {} -whichsamplesfile {} -callfile {} -annot_dir {} -annotfile {} -binsize 0 -merr -prog pl -uninf -outputdir {} > {};'.format(
                perl, linkdatagen_mps, linkdatagen_pedfile, linkdatagen_whichsamplesfile, brlmm, annotHapMap2_dir, annotHapMap2_file, sample_plink_dir, linkdatagen_mps_stdout)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=vcf2linkdatagen_jobid)
            linkdatagen_mps_jobid = job.jobId
            hold_linkdatagen_jids += [job.jobId]
            
            #plink makebed
            job_name = my_name+'_plink_makebed_'+os.path.basename(sample_name)
            cmd = 'cd {}_plink/; {} --file plink --make-bed --out {} > {};'.format(sample_plink_dir, plink, sample_name, plink_makebed_stdout)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=linkdatagen_mps_jobid)
            plink_makebed_jobid = job.jobId
            hold_linkdatagen_jids += [job.jobId]
            
            #plink homozyg
            job_name = my_name+'_plink_makebed_'+os.path.basename(sample_name)
            cmd = 'cd {}_plink/; {} --bfile {} --homozyg > {};'.format(sample_plink_dir, plink, sample_name, plink_homozyg_stdout)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=plink_makebed_jobid)
            plink_homozyg_jobid = job.jobId
            hold_linkdatagen_jids += [job.jobId]

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'Find blocks of homozygosity from bam files using linkdatagen and plink',
        epilog = '{} version {} {}'.format(name, version, copyright))
    parser.add_argument('--input', '-i', metavar='BAM', nargs='+', required=True,
        help='input list of bam files')
    parser.add_argument('--linkdatagen_dir',
        help='linkdatagen output directory (default: dirs["linkdatagen"]')
    parser.add_argument('--plink_dir',
        help='plink output directory (default: dirs["plink"]')
    parser.add_argument('--job_dir', '-j',
        help='parent job directory (default: dirs[jobs]')
    parser.add_argument('--my_name', default='bam2homozygosity_blocks',
        help='job name prefix')
    parser.add_argument('--hold_jids', nargs='*', default=None,
        help='job ids to wait for before running')
    parser.add_argument('--python',
        help='python executable (default: from config->DEFAULT.python)')
    parser.add_argument('--perl',
        help='perl executable (default: from config->DEFAULT.perl)')
    parser.add_argument('--reference', '-r',
        help='path to reference sequence fasta file (default: from config->reference.default_reference)')
    parser.add_argument('--samtools',
        help='path to samtools executable (default: from config->samtools.samtools)')
    parser.add_argument('--bcftools',
        help='path to bcftools executable (default: from config->samtools.bcftools)')
    parser.add_argument('--annotHapMap2L_b37',
        help='path to linkdatagen annotHapMap2L_b37 file (default: from config->linkdatagen.annotHapMap2L_b37)')
    parser.add_argument('--annotHapMap2',
        help='path to linkdatagen annotHapMap2 file (default: from config->linkdatagen.annotHapMap2)')
    parser.add_argument('--vcf2linkdatagen',
        help='path to vcf2linkdatagen.pl executable (default: from config->linkdatagen.vcf2linkdatagen)')
    parser.add_argument('--linkdatagen_mps',
        help='path to linkdatagen_mps.pl executable (default: from config->linkdatagen.linkdatagen_mps)')
    parser.add_argument('--plink',
        help='path to plink executable (default: from config->plink.plink)')
    parser.add_argument('--overwrite_all_output_files', action='store_true', default=False,
        help='overwriting any existing valid files')
    parser.add_argument('--overwrite_pileups', action='store_true', default=False,
        help='redo sample pileups, overwriting existing valid sample pileup files')

    args = parser.parse_args()

    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    args.python = args.python if args.python else config.get('DEFAULT','python')
    if not my.is_exe(args.python):
        raise Bam2HomozygosityBlocksError('{} is not the python executable'.format(args.python))
    args.perl = args.perl if args.perl else config.get('DEFAULT','perl')
    if not my.is_exe(args.perl):
        raise Bam2HomozygosityBlocksError('{} is not the perl executable'.format(args.perl))
    args.reference = args.reference if args.reference else config.get('reference','default_decoy_reference')
    if not my.file_exists(args.reference):
        raise Bam2HomozygosityBlocksError('cannot find reference sequence fasta file {}'.format(args.reference))
    args.samtools = args.samtools if args.samtools else config.get('samtools','samtools')
    if not my.is_exe(args.samtools):
        raise Bam2HomozygosityBlocksError('{} is not a samtools executable'.format(args.samtools))
    args.bcftools = args.bcftools if args.bcftools else config.get('samtools','bcftools')
    if not my.is_exe(args.bcftools):
        raise Bam2HomozygosityBlocksError('{} is not a bcftools executable'.format(args.bcftools))
    args.annotHapMap2L_b37 = args.annotHapMap2L_b37 if args.annotHapMap2L_b37 else config.get('linkdatagen','annotHapMap2L_b37')
    if not my.file_exists(args.annotHapMap2L_b37):
        raise Bam2HomozygosityBlocksError('{} is not the annotHapMap2L_b37 file'.format(args.annotHapMap2L_b37))
    args.annotHapMap2 = args.annotHapMap2 if args.annotHapMap2 else config.get('linkdatagen','annotHapMap2')
    if not my.file_exists(args.annotHapMap2):
        raise Bam2HomozygosityBlocksError('{} is not the annotHapMap2 file'.format(args.annotHapMap2))
    annotHapMap2_dir = os.path.dirname(args.annotHapMap2)
    annotHapMap2_file = os.path.basename(args.annotHapMap2)
    args.vcf2linkdatagen = args.vcf2linkdatagen if args.vcf2linkdatagen else config.get('linkdatagen','vcf2linkdatagen')
    if not my.file_exists(args.vcf2linkdatagen):
        raise Bam2HomozygosityBlocksError('{} is not the vcf2linkdatagen perl script'.format(args.vcf2linkdatagen))
    args.linkdatagen_mps = args.linkdatagen_mps if args.linkdatagen_mps else config.get('linkdatagen','linkdatagen_mps')
    if not my.file_exists(args.linkdatagen_mps):
        raise Bam2HomozygosityBlocksError('{} is not the linkdatagen_mps perl script'.format(args.linkdatagen_mps))
    args.plink = args.plink if args.plink else config.get('plink','plink')
    if not my.file_exists(args.plink):
        raise Bam2HomozygosityBlocksError('{} is not the plink executable'.format(args.plink))

    args.job_dir = args.job_dir if args.job_dir else my.unique_dir(args.my_name, dirs['jobs'])
    args.plink_dir = args.plink_dir if args.plink_dir else dirs['plink']
    args.linkdatagen_dir = args.linkdatagen_dir if args.linkdatagen_dir else dirs['linkdatagen']

    run(input=args.input, linkdatagen_dir=args.linkdatagen_dir, plink_dir=args.plink_dir, job_dir=args.job_dir,
        my_name=args.my_name, hold_jids=args.hold_jids, 
        python=args.python, perl=args.perl, reference=args.reference, samtools=args.samtools,
        bcftools=args.bcftools, annotHapMap2L_b37=args.annotHapMap2L_b37, annotHapMap2=args.annotHapMap2,
        vcf2linkdatagen=args.vcf2linkdatagen, linkdatagen_mps=args.linkdatagen_mps, plink=args.plink,
        overwrite_all_output_files=args.overwrite_all_output_files, overwrite_pileups=args.overwrite_pileups)


if __name__ == "__main__": sys.exit(main())
