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
from warnings import warn
import shutil
import my
import job

class RgBams2LibraryBamError(Exception): pass

#1 bam
#--top_dir /data/storage-1-02/archive/myourshaw/gmd --parent_jobName rgbamss2librarybams_metadata.txt_07202011163310 --parent_job_dir /data/storage-1-02/archive/myourshaw/gmd/job_info/rgbamss2librarybams_metadata.txt_foo --bams /data/storage-1-02/archive/myourshaw/gmd/reads/readgroup_bams/GMD151.GMD151_XT_16.110511_SN430_0243_B817FLABXX.7.ATCACGA.novoalign.bam --library_bam /scratch0/tmp/myourshaw/gmd/reads/library_bams/GMD151.GMD151_XT_16.bam --library_markdup_bam /scratch0/tmp/myourshaw/gmd/reads/library_bams/GMD151.GMD151_XT_16.markdup.bam

#6 bams
#--top_dir /data/storage-1-02/archive/myourshaw/gmd --parent_jobName rgbamss2librarybams_metadata.txt_07202011163310 --parent_job_dir /data/storage-1-02/archive/myourshaw/gmd/job_info/rgbamss2librarybams_metadata.txt_foo --bams /data/storage-1-02/archive/myourshaw/gmd/reads/readgroup_bams/GMD154A.GMD154A_XT_27.110511_SN430_0243_B817FLABXX.1.TTAGGCA.novoalign.bam /data/storage-1-02/archive/myourshaw/gmd/reads/readgroup_bams/GMD154A.GMD154A_XT_27.110511_SN430_0243_B817FLABXX.4.TTAGGCA.novoalign.bam /data/storage-1-02/archive/myourshaw/gmd/reads/readgroup_bams/GMD154A.GMD154A_XT_27.110511_SN430_0243_B817FLABXX.5.TTAGGCA.novoalign.bam /data/storage-1-02/archive/myourshaw/gmd/reads/readgroup_bams/GMD154A.GMD154A_XT_27.110511_SN430_0243_B817FLABXX.2.TTAGGCA.novoalign.bam /data/storage-1-02/archive/myourshaw/gmd/reads/readgroup_bams/GMD154A.GMD154A_XT_27.110511_SN430_0243_B817FLABXX.3.TTAGGCA.novoalign.bam /data/storage-1-02/archive/myourshaw/gmd/reads/readgroup_bams/GMD154A.GMD154A_XT_27.110623_SN860_0067_2011-100R_A81MVKABXX.8.TTAGGCA.novoalign.bam --library_bam /scratch0/tmp/myourshaw/gmd/reads/library_bams/GMD154A.GMD154A_XT_27.bam --library_markdup_bam /scratch0/tmp/myourshaw/gmd/reads/library_bams/GMD154A.GMD154A_XT_27.markdup.bam

def main():

	#command line arguments
	parser = argparse.ArgumentParser(parents=[my.get_directory_parser()],
		description = 'From a metadata file, merge multiple readgroup bam files by library and remove duplicates, calculate metrics, validate bam files',
		epilog = 'pypeline.rgs2bams version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--parent_jobName', required=True,
											help='job name of calling script')
	parser.add_argument('--parent_job_dir', required=True,
											help='job directory of calling script')
	parser.add_argument('--bams', nargs='+',
											help='input readgroup bam file(s)',)
	parser.add_argument('--library_bam', required=True,
											help='output library bam file',)
	parser.add_argument('--library_markdup_bam', required=True,
											help='output library markdup bam file',)
	parser.add_argument('--picard_MergeSamFiles',
											help='path to picard picard_MergeSamFiles jar')
	parser.add_argument('--picard_MarkDuplicates',
											help='path to picard MarkDuplicates jar')
	parser.add_argument('--picard_ValidateSamFile',
											help='path to picard ValidateSamFile jar')
	parser.add_argument('--ref',
											help='path to reference genome fasta file')
	parser.add_argument('--read_name_regex', default='\"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*\"',
											help='read name regex')
	parser.add_argument('--optical_duplicate_pixel_distance', default=100,
											help='optical duplicate pixel distance')
	args = parser.parse_args()
	my.set_directory_defaults(args)
	config = my.get_config(args)
		
	this_job_name = 'rgbams2librarybam_{}_{}'.format(os.path.basename(args.library_markdup_bam), my.localtime_squish())
	this_job_dir = my.unique_dir(this_job_name, args.parent_job_dir)

	if not args.picard_MergeSamFiles: args.picard_MergeSamFiles = config.get('picard','MergeSamFiles')
	if not my.file_exists(args.picard_MergeSamFiles):
		raise Rg2BamError('cannot find picard MergeSamFiles jar')
	if not args.picard_MarkDuplicates: args.picard_MarkDuplicates = config.get('picard','MarkDuplicates')
	if not my.file_exists(args.picard_MarkDuplicates):
		raise Rg2BamError('cannot find picard MarkDuplicates jar')
	if not args.picard_ValidateSamFile: args.picard_ValidateSamFile = config.get('picard','ValidateSamFile')
	if not my.file_exists(args.picard_ValidateSamFile):
		raise Rg2BamError('cannot find picard ValidateSamFile jar')
	if not args.ref: args.ref = config.get('reference','b37_fasta')
	if not my.file_exists(args.ref):
		raise Rg2BamError('cannot find reference genome for picard ValidateSamFile')
	args.bams = sorted(args.bams)
	if not args.read_name_regex : args.read_name_regex = "\"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*\""
	if not args.optical_duplicate_pixel_distance : args.optical_duplicate_pixel_distance = 100
	java_str = 'java -Xmx5g -Djava.io.tmpdir={} -jar'.format(args.tmp_dir)
	picard_std_options = 'VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=1000000'
	picard_MergeSamFiles_jobName = 'MergeSamFiles_{}'.format(os.path.basename(args.library_bam))
	picard_MergeSamFiles_job_dir = my.unique_dir(picard_MergeSamFiles_jobName, this_job_dir)
	picard_MarkDuplicates_jobName = 'MarkDuplicates_{}'.format(os.path.basename(args.library_markdup_bam))
	picard_MarkDuplicates_job_dir = my.unique_dir(picard_MarkDuplicates_jobName, this_job_dir)
	picard_ValidateSamFile_jobName = 'ValidateSamFile_{}'.format(os.path.basename(args.library_markdup_bam))
	picard_ValidateSamFile_job_dir = my.unique_dir(picard_ValidateSamFile_jobName, this_job_dir)

	#picard MergeSamFiles
	if len(args.bams) == 1:
		shutil.copy(args.bams[0], args.library_bam)
	else:
		picard_jar = args.picard_MergeSamFiles
		cmd =  '{} {} {} CREATE_INDEX=true CREATE_MD5_FILE=true SORT_ORDER=coordinate ASSUME_SORTED=false MERGE_SEQUENCE_DICTIONARIES=false USE_THREADING=true INPUT={} OUTPUT={};'\
					.format(java_str, picard_jar, picard_std_options, ' INPUT='.join(args.bams), args.library_bam)
		
		jobName = picard_MergeSamFiles_jobName
		job_dir = picard_MergeSamFiles_job_dir
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

	#picard MarkDuplicates
	picard_jar = args.picard_MarkDuplicates
	metrics_file = os.path.join(args.bam_metrics_dir, '{}.markdup.metrics'.format(os.path.basename(args.library_markdup_bam)))
	cmd =  '{} {} {} CREATE_INDEX=true CREATE_MD5_FILE=true REMOVE_DUPLICATES=false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 READ_NAME_REGEX={} OPTICAL_DUPLICATE_PIXEL_DISTANCE={} INPUT={} OUTPUT={} METRICS_FILE={};'\
				.format(java_str, picard_jar, picard_std_options, args.read_name_regex, args.optical_duplicate_pixel_distance, args.library_bam, args.library_markdup_bam, metrics_file)
	
	jobName = picard_MarkDuplicates_jobName
	job_dir = picard_MarkDuplicates_job_dir
	my.print_log('job {} started; command:\n{}'.format(jobName, cmd))
	j = job.Job(jobName = jobName,
	outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = job_dir,
	email = args.email, blockEmail = not args.sendemailforintermediatejobs,
	processors = 8, memory = '24G', queuingModeBunchJobsOnNodes=False)
	try:
		j.executeCommand(cmd, synchronous=True)
	except Exception as e:
		my.print_err('There was an error running job {}\n{}\nintermediate files are in {}\n'.format(jobName, e, job_dir))
	else:
		with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
			s.write(format('\n'.join(j.completionMsg)))

	#picard ValidateSamFile
	picard_jar = args.picard_ValidateSamFile
	this_bam = args.library_markdup_bam
	this_validate = '{}.validate'.format(args.library_markdup_bam)
	cmd =  '{} {} {} MODE=VERBOSE REFERENCE_SEQUENCE={} VALIDATE_INDEX=true INPUT={} OUTPUT={};'\
				.format(java_str, picard_jar, picard_std_options, args.ref, this_bam, this_validate)

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
		with open(this_validate) as v:
			val_msg = v.readline().rstrip('\n').lower()
		if val_msg != 'no errors found':
			raise RgBams2LibraryBamError('{} picard.ValidateSamFile validation errrors\n{}\ninermediate files in {}'.format(this_bam, val_msg, args.parent_job_dir))
		#elif args.deletetempfiles:
		#	os.remove(args.library_bam)

			
if __name__ == "__main__": sys.exit(main())
