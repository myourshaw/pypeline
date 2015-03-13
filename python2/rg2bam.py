#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import stat
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import glob
import tempfile
from time import localtime, strftime
import re
import my
import job

class Rg2BamError(Exception): pass


#--top_dir /scratch0/tmp/myourshaw/gmd --parent_jobName rgs2bams_metadata.txt_07182011142549 --parent_job_dir /scratch0/tmp/myourshaw/gmd/job_info/rgs2bams_metadata.txt_07182011142549_WVx4cb --qseqs /scratch0/tmp/myourshaw/gmd/reads/readgroup_qseqs/110511_SN430_0243_B817FLABXX.1.1.TTAGGCA_qseq.txt.gz /scratch0/tmp/myourshaw/gmd/reads/readgroup_qseqs/110511_SN430_0243_B817FLABXX.1.3.TTAGGCA_qseq.txt.gz --novoalign_output_sam /scratch0/tmp/myourshaw/gmd/reads/readgroup_bams/GMD154A.GMD154A_XT_27.110511_SN430_0243_B817FLABXX.1.TTAGGCA.novoalign.sam --readgroup_output_bam /scratch0/tmp/myourshaw/gmd/reads/readgroup_bams/GMD154A.GMD154A_XT_27.110511_SN430_0243_B817FLABXX.1.TTAGGCA.novoalign.bam --RG "@RG\tID:110511_SN430_0243_B817FLABXX.1.TTAGGCA\tCN:UCLA-Nelson\tDT:2011-05-11\tLB:GMD154A_XT_27\tPL:ILLUMINA\tPU:HWI-ST430.243.1.TTAGGCA\tSM:GMD154A\tDS:sample=GMD154A;library=GMD154A_XT_27;sequencing_center=UCLA-Nelson;platform=ILLUMINA;platform_unit=HWI-ST430.243.1.TTAGGCA;run_date=2011-05-11;experiment_run_id=110511_SN430_0243_B817FLABXX;lane=1;barcode=TTAGGCA;predicted_median_insert_size=184" --adapters AGATCGGAAGAGCACACGTCT AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA


#one or two qseq files containing a single readgroup -> readgroup bam file ready to be merged into other bamfiles of same library
#indended to be run as a job array, input parameters pre-validated by calling program
#submits each step as a job job array with one command and waits for completion

def main():

	#command line arguments
	parser = argparse.ArgumentParser(parents=[my.get_directory_parser()],
		description = 'Align qseq readgroup files with novoalign, fix mate information, validate sam file, verify validation, calculate pre-rmdup readgroup bam metrics, delete intermediate sam file produced by novoalign',
		epilog = 'pypeline.rg2bam version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--parent_jobName', required=True,
											help='job name of calling script')
	parser.add_argument('--parent_job_dir', required=True,
											help='job directory of calling script')
	parser.add_argument('--qseqs', nargs='*',
											help='one or two qseq[.gz] files')
	parser.add_argument('--novoalign_output_sam', required=True,
											help='novoalign output sam file')
	parser.add_argument('--readgroup_output_bam', required=True,
											help='readgroup output sam file')
	parser.add_argument('--novoalign',
											help='path to novoalign executable')
	parser.add_argument('--index',
											help='path to novoindex-produced index')
	parser.add_argument('--RG', required=True,
											help='@RG string for novoalign')
	parser.add_argument('--adapters', nargs='*', required=True,
											help='sequencing adapters')
	parser.add_argument('--picard_FixMateInformation',
											help='path to picard FixMateInformation jar')
	parser.add_argument('--picard_ValidateSamFile',
											help='path to picard ValidateSamFile jar')
	parser.add_argument('--ref',
											help='path to reference genome fasta file')
	parser.add_argument('--NoNovoalign', action='store_true', default=False,
											help='do not run novoalign')
	parser.add_argument('--NoFixMateInformation', action='store_true', default=False,
											help='do not run FixMateInformation')
	parser.add_argument('--NoValidateSamFile', action='store_true', default=False,
											help='do not run ValidateSamFile')
	args = parser.parse_args()
	my.set_directory_defaults(args)
	config = my.get_config(args)

	readgroup_id = [r for r in args.RG.split('\\t') if r.startswith('ID:')][0].split(':')[1]

	this_job_name = 'rg2bam_{}_{}'.format(readgroup_id, my.localtime_squish())
	this_job_dir = my.unique_dir(this_job_name, args.parent_job_dir)

	novoalign_jobName = 'novoalign_{}'.format(readgroup_id)
	novoalign_job_dir = my.unique_dir(novoalign_jobName, this_job_dir)
	picard_FixMateInformation_jobName = 'FixMateInformation_{}'.format(readgroup_id)
	picard_FixMateInformation_job_dir = my.unique_dir(picard_FixMateInformation_jobName, this_job_dir)
	picard_ValidateSamFile_jobName = 'ValidateSamFile_{}'.format(readgroup_id)
	picard_ValidateSamFile_job_dir = my.unique_dir(picard_ValidateSamFile_jobName, this_job_dir)

	
	#novoalign
	if not args.NoNovoalign:
		novoalign = args.novoalign if args.novoalign else config.get('novocraft','novoalign')
		if not my.is_exe(novoalign):
			raise Rg2BamError('{} is not a novoalign executable'.format(novoalign))
		index = args.index if args.index else config.get('novocraft','default_index')
		if not my.file_exists(index):
			raise Rg2BamError('cannot find index file {}'.format(index))
		himem = bool(os.path.getsize(index) > 6500000000)
		qseq_info = my.qseq_peek(args.qseqs[0])
		read_length = qseq_info['read_length']
		read_length = 50 if read_length <= 50 else 75 if read_length <= 75 else 100 if read_length <= 100 else read_length
		min_quality_bases = int(round(read_length/2.0)) #novoalign default is ~20; this yields 25 for 50 base reads and 50 for 100 base reads
		
		cmd = '{} -k -o SAM "{}" -d {} -a {} -F QSEQ -l {} -H --hdrhd 1 -f {} > {};'\
				 .format(novoalign, args.RG, index, ' '.join(args.adapters), min_quality_bases, ' '.join(args.qseqs), args.novoalign_output_sam)
		
		jobName = novoalign_jobName
		job_dir = novoalign_job_dir
		my.print_log('job {} started; command:\n{}'.format(jobName, cmd))
		j = job.Job(jobName = jobName,
		outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = job_dir,
		email = args.email, blockEmail = not args.sendemailforintermediatejobs,
		processors = 8, memory = '24G' if himem else '7G')
		try:
			j.executeCommand(cmd, synchronous=True)
		except Exception as e:
			my.print_err('There was an error running job {}\n{}\nintermediate files are in {}\n'.format(jobName, e, job_dir))
		else:
			with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
				s.write(format('\n'.join(j.completionMsg)))

		java_str = 'java -Xmx5g -Djava.io.tmpdir={} -jar'.format(args.tmp_dir)
		picard_std_options = 'VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=1000000'
	
	#picard FixMateInformation
	if not args.NoFixMateInformation:
		picard_jar = args.picard_FixMateInformation if args.picard_FixMateInformation else config.get('picard','FixMateInformation')
		if not my.file_exists(picard_jar):
			raise Rg2BamError('cannot find picard FixMateInformation jar')
		cmd =  '{} {} {} CREATE_INDEX=true CREATE_MD5_FILE=true SORT_ORDER=coordinate INPUT={} OUTPUT={};'\
					.format(java_str, picard_jar, picard_std_options, args.novoalign_output_sam, args.readgroup_output_bam)
		
		jobName = picard_FixMateInformation_jobName
		job_dir = picard_FixMateInformation_job_dir
		my.print_log('job {} started; command:\n{}'.format(jobName, cmd))
		j = job.Job(jobName = jobName,
		outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = job_dir,
		email = args.email, blockEmail = not args.sendemailforintermediatejobs,
		processors = 8, queuingModeBunchJobsOnNodes=False)
		try:
			j.executeCommand(cmd, synchronous=True)
		except Exception as e:
			my.print_err('There was an error running job {}\n{}\nintermediate files are in {}\n'.format(jobName, e, job_dir))
		else:
			with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
				s.write(format('\n'.join(j.completionMsg)))

	#picard ValidateSamFile
	if not args.NoFixMateInformation and not args.NoValidateSamFile:
		readgroup_output_bam_validate = '{}.validate'.format(args.readgroup_output_bam)
		picard_jar = args.picard_ValidateSamFile if args.picard_ValidateSamFile else config.get('picard','ValidateSamFile')
		if not my.file_exists(picard_jar):
			raise Rg2BamError('cannot find picard ValidateSamFile jar')
		ref = args.ref if args.ref else config.get('reference','b37_fasta')
		if not my.file_exists(ref):
			raise Rg2BamError('cannot find reference genome for picard ValidateSamFile')
		cmd =  '{} {} {} MODE=VERBOSE REFERENCE_SEQUENCE={} VALIDATE_INDEX=true INPUT={} OUTPUT={};'\
					.format(java_str, picard_jar, picard_std_options, ref, args.readgroup_output_bam, readgroup_output_bam_validate)
		
		jobName = picard_ValidateSamFile_jobName
		job_dir = picard_ValidateSamFile_job_dir
		j = job.Job(jobName = jobName,
		outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = job_dir,
		email = args.email, blockEmail = not args.sendemailforintermediatejobs)
		try:
			j.executeCommand(cmd, synchronous=True)
		except Exception as e:
			my.print_err('There was an error running job {}\n{}\nintermediate files are in {}\n'.format(jobName, e, job_dir))
		else:
			with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
				s.write(format('\n'.join(j.completionMsg)))

			#verify validation
			with open(readgroup_output_bam_validate) as v:
				val_msg = v.readlines().rstrip('\n').lower()
				if val_msg != 'no errors found':
					raise Rg2BamError('readgroup {} picard.ValidateSamFile validation errrors\n{}\ninermediate files in {}'.format(readgroup_id, val_msg, args.parent_job_dir))
				elif args.deletetempfiles:
					os.remove(args.novoalign_output_sam)


if __name__ == "__main__": sys.exit(main())
