#!/usr/bin/env python

import os
import sys
import stat
import argparse
import ConfigParser #configparser in python 3
import glob
import tempfile
from time import localtime, strftime
import re
import drmaa
from my import util
from getpass import getuser

#data
#-b /share/apps/myourshaw/resources/truseq_xt_barcode03.tag_list --NC_OFF -d /scratch1/tmp/myourshaw/gmd/new/qseq/unbarcode /data/storage-1-02/solexa/110511_SN430_0243_B817FLABXX/Data/Intensities/BaseCalls/s_[1-5]_1_*_qseq.txt
#-d /scratch1/tmp/myourshaw/gmd/new/qseq/unbarcode /data/storage-1-02/solexa/110511_SN430_0243_B817FLABXX/Data/Intensities/BaseCalls/s_[6-7]_1_*_qseq.txt /data/storage-1-04/archive/illumina/110616_SN860_0065_2011-101_B817EAABXX/Data/Intensities/BaseCalls/QSEQ/s_[1-7]_1_*_qseq.txt
#test
#-b /share/apps/myourshaw/resources/truseq_xt_barcode03.tag_list --NC_OFF -d /scratch1/tmp/myourshaw/gmd/new/qseq/test-unbarcode /data/storage-1-02/solexa/110511_SN430_0243_B817FLABXX/Data/Intensities/BaseCalls/s_1_1_0001_qseq.txt
#-d /scratch1/tmp/myourshaw/gmd/new/qseq/test-unbarcode /data/storage-1-02/solexa/110511_SN430_0243_B817FLABXX/Data/Intensities/BaseCalls/s_6_1_0001_qseq.txt /data/storage-1-04/archive/illumina/110616_SN860_0065_2011-101_B817EAABXX/Data/Intensities/BaseCalls/QSEQ/s_1_1_1101_qseq.txt
#--NoUnbarcodeDemux -d /scratch1/tmp/myourshaw/gmd/new/qseq/unbarcode

	
def main():
	
	USER = getuser()
	HERE = os.path.dirname(sys.argv[0])
	PROGRAM = 'qseq2bam'
	print_log('{} \xa9 2011 Michael Yourshaw all rights reserved'.format(PROGRAM))
	print_log('user: {}'.format(USER))
	print_log('command: {}'.format(' '.join(sys.argv)))
	
	#command line arguments
	parser = argparse.ArgumentParser(
		description = 'Use novobarcode to segregate qseq files by barcode and merge tiles', # main description for help
		epilog = '\xc2 2011 Michael Yourshaw all rights reserved') # displayed after help
	parser.add_argument('--ini',nargs='?', default=HERE+'/b37.ini',
											help='.ini configuration file (default: '+HERE+'/b37.ini')
	parser.add_argument('--email','-m', metavar='EMAIL', nargs='+', default=USER+'@ucla.edu',
											help='email addresses (default: '+USER+'@ucla.edu)')
	parser.add_argument('--NoEmail', action='store_true', default=False, help='do not send email notifications')
	parser.add_argument('--novobarcode','-n', nargs='?',
											help='novobarcode executable (default: b37.ini->novobarcode)')
	parser.add_argument('--barcodes','-b', nargs='?',
											help='barcode tag list (default: b37.ini->truseq_xt_barcodes)')
	parser.add_argument('--NC_OFF', action='store_const', const='--NC_OFF', default='',
											help='Turns off creation of NC folder and the writing of the unclassified reads (default: on)')
	parser.add_argument('--flowcell','-f', nargs='?',
											help='flowcell (default: from qseq directory)')
	parser.add_argument('--jobname','-N', nargs='?',default='unbarcode',
											help='job name (default: unbarcode)')
	parser.add_argument('--output_dir', '-d', nargs='?', default=os.getcwd()+'/unbarcode',
											help='output directory (default: '+os.getcwd()+'/unbarcode)')
	parser.add_argument('qseqs', metavar='QSEQ', nargs='*',
											help='input list of paired end 1 files (s_[1-8]_1_[1-3]_qseq.txt)')
	parser.add_argument('--NoUnbarcodeDemux', action='store_true', default=False, help='do not demultiplex qseqs by barcode')
	parser.add_argument('--NoUnbarcodeMerge', action='store_true', default=False, help='do not merge demuxed qseqs')
	parser.add_argument('--unbarcode_merged_dirname', nargs='?', default='unbarcode_merged_qseq',
											help='name of directory for merged qseqs (default: unbarcode_merged_qseq)')
	parser.add_argument('--unbarcode_unmerged_dirname', nargs='?', default='unbarcode_unmerged_qseq',
											help='name of directory for unmerged qseqs (default: unbarcode_unmerged_qseq)')

	args = parser.parse_args()
	
	#get defaults from configuration ini file
	config = ConfigParser.SafeConfigParser()
	config.read(args.ini)
	if args.novobarcode == None: args.novobarcode = config.get('novocraft','novobarcode')
	if args.barcodes == None: args.barcodes = config.get('novocraft','truseq_xt_barcodes')

	allq = config.get('qsub','allq')
	himemq = config.get('qsub','himemq')
	tempdir = os.path.join(config.get('DEFAULT','tmpscratch'),USER)

	FILE_MISSING_MSG = 'These required files are missing. You may need to edit {0}.\n\t{1}'
	RUN_JOB_MSG = """Job: {0} [{1}]
Command: {2}
STDOUT and STDERR in {3}.
Notification to {4}."""
	EXCEPTION_MSG = 'ERROR: {0}; {1}'

	#output directories
	if args.output_dir == None: args.output_dir = os.path.dirname(os.path.abspath(args.qseqs[0]))
	makedir(args.output_dir)
	unbarcode_tmp_dir = unique_dir(args.output_dir, 'unbarcode_tmp')
	unbarcode_merged_dir = os.path.join(args.output_dir,'unbarcode_merged')\
		if args.unbarcode_merged_dirname == None else os.path.join(args.output_dir,args.unbarcode_merged_dirname)
	makedir(unbarcode_merged_dir)
	unbarcode_unmerged_dir = os.path.join(args.output_dir,'unbarcode_unmerged')\
		if args.unbarcode_unmerged_dirname == None else os.path.join(args.output_dir,args.unbarcode_unmerged_dirname)
	makedir(unbarcode_unmerged_dir)
	qout_dir = os.path.join(unbarcode_tmp_dir,'qout')
	makedir(qout_dir)
	
	#temp file for running job array commands
	drmaa_job_array_runner_file, drmaa_job_array_runner = unique_file(os.path.join(unbarcode_tmp_dir,'drmaa_job_array_runner.sh'))
	drmaa_job_array_runner_file.write("""#!/bin/bash

#SGE job array runner
#input: $1 is a job array file with the commands to be run, one per line
#qsub -t or drmaa runBulkJobs calls this script once for each job
#with ${SGE_TASK_ID} = 1-based line number in the job array file

eval "`head -n ${SGE_TASK_ID} ${1} | tail -n 1`";
""")
	drmaa_job_array_runner_file.close()
	os.chmod(drmaa_job_array_runner, stat.S_IRWXU)	
	

	#novobarcode(multiplexed barcoded qseq files from illumina sequencer) -> unbarcode_tmp_dir/<experiment run>/<lane>/<barcode>/<qseq files> 
	if not args.NoAlign:
