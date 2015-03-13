#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import stat
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
from time import localtime, strftime
import my
import job
import picard

class Qseq2BamError(Exception): pass


#--NoFixMateInformation --NoNovoalign --novoalign /share/apps/myourshaw/novocraft-current/novoalign --index /share/apps/myourshaw/resources/ref.b37.k15s1.nix --reference /share/apps/myourshaw/resources/human_g1k_v37.fasta --samtools /share/apps/myourshaw/samtools-current/samtools --parent_job_dir /scratch0/tmp/myourshaw/gmd/jobs/qseqs2bams_20110916180057_3HY7eC --qseqs /scratch0/tmp/myourshaw/gmd/qseqs/HWI-ST430/243/6/HWI-ST430.243.6.1.ATCACGA_qseq.txt.gz /scratch0/tmp/myourshaw/gmd/qseqs/HWI-ST430/243/6/HWI-ST430.243.6.3.ATCACGA_qseq.txt.gz --readgroup_bam /scratch0/tmp/myourshaw/gmd/bams/readgroup_bams/GMD1A.GMD1A_XT_1.HWI-ST430.243.6.ATCACGA.novoalign.bam --RG "@RG\tID:HWI-ST430.243.6.ATCACGA\tCN:UCLA-Nelson\tDT:2011-05-11\tLB:GMD1A_XT_1\tPI:174\tPL:ILLUMINA\tPU:HWI-ST430.243.B817FLABXX.6.ATCACGA\tSM:GMD1A\tDS:sample=GMD1A;library=GMD1A_XT_1;library_protocol=SureSelect_XT_Illumina;sequencing_library=201105061_XT_mix_2+5%PhiX;sequencing_center=UCLA-Nelson;platform=ILLUMINA;instrument_model=Illumina HiSeq 2000;platform_unit=HWI-ST430.243.B817FLABXX.6.ATCACGA;run_date=2011-05-11;run_folder=110511_SN430_0243_B817FLABXX;machine=HWI-ST430;run=243;flowcell=B817FLABXX;lane=6;barcode=ATCACGA;read_length=50;predicted_median_insert_size=174;adapter_1=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;adapter_2=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" --read_length 50 --adapters AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC


#one or two qseq files containing a single readgroup -> readgroup bam file ready to be merged into other bamfiles of same library
#intended to be run as a job array, input parameters pre-validated by calling program
#submits each step as a job job array with one command and waits for completion

def main():

    #command line arguments
    parser = argparse.ArgumentParser(#parents=[my.default_parser()],
        description = 'Align qseq readgroup files with novoalign, fix mate information, \
            validate sam file, verify validation, calculate pre-rmdup readgroup bam metrics',
        epilog = 'pypeline.qseq2bam version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--novoalign',
        help='path to novoalign executable (default: from config->novocraft.novoalign)')
    parser.add_argument('--index',
        help='path to novoindex-produced index (default: from config->novocraft.default_index)')
    parser.add_argument('--reference', '-r',
        help='path to reference sequence fasta file (default: from config->reference.default_reference)')
    parser.add_argument('--samtools',
        help='path to samtools executable (default: from config->samtools.samtools)')
    parser.add_argument('--parent_job_dir', required=True,
        help='job directory of calling script')
    parser.add_argument('--qseqs', '-q', nargs='+',
        help='one or two qseq[.gz] files')
    parser.add_argument('--readgroup_bam', required=True,
        help='readgroup output bam file')
    parser.add_argument('--RG', required=True,
        help='@RG string for novoalign')
    parser.add_argument('--read_length', required=True, type=int,
        help='number of bases in qseq read sequence')
    parser.add_argument('--adapters', nargs=2, required=True,
        help='sequencing adapters')
    parser.add_argument('--NoNovoalign', action='store_true', default=False,
        help='do not run novoalign')
    parser.add_argument('--NoFixMateInformation', action='store_true', default=False,
        help='do not run FixMateInformation')
    parser.add_argument('--NoValidateSamFile', action='store_true', default=False,
        help='do not run ValidateSamFile')

    args = parser.parse_args()
    config = my.get_config(args)
    
    args.novoalign = args.novoalign if args.novoalign else config.get('novocraft','novoalign')
    if not my.is_exe(args.novoalign):
        raise Qseq2BamError('{} is not a novoalign executable'.format(args.novoalign))
    args.index = args.index if args.index else config.get('novocraft','default_index')
    if not my.file_exists(args.index):
        raise Qseq2BamError('cannot find novolaign genome index file {}'.format(args.index))
    himem = bool(os.path.getsize(args.index) > 6500000000)
    args.reference = args.reference if args.reference else config.get('reference','default_reference')
    if not my.file_exists(args.reference):
        raise Qseq2BamError('cannot find reference sequence fasta file {}'.format(args.reference))
    args.samtools = args.samtools if args.samtools else config.get('samtools','samtools')
    if not my.is_exe(args.samtools):
        raise Qseq2BamError('{} is not a samtools executable'.format(args.samtools))
    
    my_name = '{}_{}'.format('qseq2bam', my.localtime_squish())
    job_dir = my.unique_dir(my_name, args.parent_job_dir)
    
    #align qseq file(s) with novoalign to the genome specified by the --index argument
    if not args.NoNovoalign:
        min_quality_bases = int(round(args.read_length/2.0)) #novoalign default is ~20; this yields 25 for 50 base reads and 50 for 100 base reads
        cmd = '{} -k -o SAM "{}" -d {} -a {} -F QSEQ -l {} -H --hdrhd 1 -f {} | {} view -S1  -  > {};'\
            .format(args.novoalign, args.RG, args.index, ' '.join(args.adapters), min_quality_bases, ' '.join(args.qseqs), args.samtools, args.readgroup_bam)
        job_name = my_name+'_novoalign_'+os.path.basename(args.readgroup_bam)
        print 'Running job {}. Job_info, STDOUT, and STDERR will be in {}. ...'.format(job_name, job_dir)    
        job = my.run_job(cmd, job_name, job_dir, synchronous=True,
            processors = 8, memory = '24G' if himem else '7G')
        print 'Job {}, jobId {}, done. Job_info, STDOUT, and STDERR are in {}.'.format(job_name, job.jobId, job_dir)
    
    #fix mate information with picard, overwriting the bam file created by novoalign
    #this step sorts the bam file, creates an index,
    #and ensures that all mate-pair information is in sync between each read and its mate pair.
    if not args.NoFixMateInformation:
        job_name = my_name+'_FixMateInformation_'+os.path.basename(args.readgroup_bam)
        print 'Running job {}. Job_info, STDOUT, and STDERR will be in {}. ...'.format(job_name, job_dir)    
        picard_args = {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate', 'INPUT': args.readgroup_bam}
        job = picard.run(tool='FixMateInformation', picard_args=picard_args, job_dir=job_dir, synchronous=True)
        print 'Job {}, jobId {}, done. Job_info, STDOUT, and STDERR are in {}.'.format(job_name, job.jobId, job_dir)
    
    #validate bam file with picard
    if not args.NoValidateSamFile:
        job_name = my_name+'_ValidateSamFile_'+os.path.basename(args.readgroup_bam)
        print 'Running job {}. Job_info, STDOUT, and STDERR will be in {}. ...'.format(job_name, job_dir)    
        readgroup_bam_validate = args.readgroup_bam + '.validate'
        picard_args = {'INPUT': args.readgroup_bam, 'OUTPUT': readgroup_bam_validate}
        job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=args.reference, synchronous=True)
        rg['ValidateSamFile_jobid'] = job.jobId
        print 'Job {}, jobId {}, done. Job_info, STDOUT, and STDERR are in {}.'.format(job_name, job.jobId, job_dir)
        
        #verify validation
        with open(readgroup_bam_validate) as v:
            val_msg = v.readlines().rstrip('\n').lower()
            if val_msg != 'no errors found':
                raise Qseq2BamError('{} had validation errrors\n{}'.format(args.readgroup_bam, val_msg))


if __name__ == "__main__": sys.exit(main())
