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
import csv
import my
import job
import rgbams2librarybam

class RgBams2LibraryBamsError(Exception): pass

#-m /home/myourshaw/lab/git-myourshaw/python/metadata.txt -d /data/storage-1-02/archive/myourshaw/gmd

def main():

	#command line arguments
	parser = argparse.ArgumentParser(parents=[my.get_directory_parser()],
		description = 'From a metadata file, merge multiple readgroup bam files by library and remove duplicates, calculate metrics, validate bam files',
		epilog = 'pypeline.rgs2bams version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--metadata', '-m', required=True,
											help='metadata file',)
	parser.add_argument('--readgroups', nargs='*',
											help='readgroup_ids from metadata file to process, "" selects blank readgroup_ids (default: all)',)
	parser.add_argument('--parent_jobName',
											help='job name of calling script')
	parser.add_argument('--parent_job_dir',
											help='job directory of calling script')
	parser.add_argument('--picard_MergeSamFiles',
											help='path to picard picard_MergeSamFiles jar')
	parser.add_argument('--picard_MarkDuplicates',
											help='path to picard MarkDuplicates jar')
	parser.add_argument('--picard_ValidateSamFile',
											help='path to picard ValidateSamFile jar')
	args = parser.parse_args()
	my.set_directory_defaults(args)
	config = my.get_config(args)
		
	if not args.parent_job_dir: args.parent_job_dir = args.job_info_dir
	this_jobName = 'rgbams2librarybams_{}_{}'.format(os.path.basename(args.metadata), my.localtime_squish())
	this_job_dir = my.unique_dir(this_jobName, args.parent_job_dir)

	cmds = []
	
	if not my.file_exists(args.metadata):
		raise Rg2BamError('metadata file {} does not exist'.format(args.metadata))
	dictReader = csv.DictReader(open(args.metadata,'rb'), dialect=csv.excel_tab)
	metadata_list = [m for m in csv.DictReader(open(args.metadata,'rb'), dialect=csv.excel_tab) if not args.readgroups or m['readgroup_id'] in args.readgroups]
	if not metadata_list:
		raise Rg2BamError('no records selected from metadata file {}'.format(args.metadata))
	samples_libraries = list({(m['sample'], m['library']) for m in metadata_list})
	for sample_library in samples_libraries:
		sample_library_readgroup_bams = list({m['readgroup_bam'] for m in metadata_list if (m['sample'], m['library']) == sample_library})
		library_bam = os.path.join(args.library_bams_dir, '{}.bam'.format('.'.join(sample_library)))
		library_markdup_bam = os.path.join(args.library_bams_dir, '{}.markdup.bam'.format('.'.join(sample_library)))
		cmds += ['python {} --top_dir {} --parent_jobName {} --parent_job_dir {} --bams {} --library_bam {} --library_markdup_bam {}'
						 .format(rgbams2librarybam.__file__, args.top_dir, this_jobName, this_job_dir, ' '.join(sample_library_readgroup_bams), library_bam, library_markdup_bam)]

	jobName = this_jobName
	job_dir = this_job_dir
	my.print_log('job {} started; commands:\n{}'.format(jobName, '\n'.join(cmds)))
	j = job.Job(jobName = jobName,
	outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = job_dir,
	email = args.email, blockEmail = not args.sendemailforintermediatejobs)
	try:
		j.executeCommands(cmds, synchronous=True)
	except Exception as e:
		my.print_err('There was an error running the job {}:\n{}\nintermediate files are in {}'.format(jobName, e, job_dir))
	else:
		my.print_log('job {} finished; status:\n{}'.format(jobName, '\n'.join(j.completionMsg)))
		with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
			s.write(format('\n'.join(j.completionMsg)))
			
if __name__ == "__main__": sys.exit(main())
