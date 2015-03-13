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
from getpass import getuser

#data
#-b /share/apps/myourshaw/resources/truseq_xt_barcode03.tag_list --NC_OFF -d /scratch1/tmp/myourshaw/gmd/new/qseq/unbarcode /data/storage-1-02/solexa/110511_SN430_0243_B817FLABXX/Data/Intensities/BaseCalls/s_[1-5]_1_*_qseq.txt
#-d /scratch1/tmp/myourshaw/gmd/new/qseq/unbarcode /data/storage-1-02/solexa/110511_SN430_0243_B817FLABXX/Data/Intensities/BaseCalls/s_[6-7]_1_*_qseq.txt /data/storage-1-04/archive/illumina/110616_SN860_0065_2011-101_B817EAABXX/Data/Intensities/BaseCalls/QSEQ/s_[1-7]_1_*_qseq.txt
#test
#-b /share/apps/myourshaw/resources/truseq_xt_barcode03.tag_list --NC_OFF -d /scratch1/tmp/myourshaw/gmd/new/qseq/test-unbarcode /data/storage-1-02/solexa/110511_SN430_0243_B817FLABXX/Data/Intensities/BaseCalls/s_1_1_0001_qseq.txt
#-d /scratch1/tmp/myourshaw/gmd/new/qseq/test-unbarcode /data/storage-1-02/solexa/110511_SN430_0243_B817FLABXX/Data/Intensities/BaseCalls/s_6_1_0001_qseq.txt /data/storage-1-04/archive/illumina/110616_SN860_0065_2011-101_B817EAABXX/Data/Intensities/BaseCalls/QSEQ/s_1_1_1101_qseq.txt
#--NoUnbarcodeDemux -d /scratch1/tmp/myourshaw/gmd/new/qseq/unbarcode

class MyError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

def makedir(dir):
	try:
	  os.makedirs(dir)
	except OSError:
		if os.path.isdir(dir):
			pass
		else:
			raise

#http://stackoverflow.com/questions/183480/is-this-the-best-way-to-get-unique-version-of-filename-w-python
def unique_file(file_name):
	dirname, filename = os.path.split(file_name)
	prefix, suffix = os.path.splitext(filename)
	fd, filename = tempfile.mkstemp(suffix, prefix+"_", dirname)
	return os.fdopen(fd,'w'), filename

def list_files(pathnames):
	#list -> list
	result = []
	for pathname in pathnames:
		for path in glob.glob(os.path.expanduser(os.path.expandvars(pathname))):
			result.append(path)
	return result

def get_files(names):
	return (file for file in names if os.path.isfile(file))

def is_number(s):
	try:
		n = float(s)
		if (math.isinf(n) or math.isnan(n)): return False
		return True
	except ValueError:
		return False

#http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
def flatten(x):
	"""flatten(sequence) -> list

	Returns a single, flat list which contains all elements retrieved
	from the sequence and all recursively contained sub-sequences
	(iterables)."""

	result = []
	for el in x:
		#if isinstance(el, (list, tuple)):
		if hasattr(el, "__iter__") and not isinstance(el, basestring):
			result.extend(flatten(el))
		else:
			result.append(el)
	return result

def file_exists(f):
	if os.path.isfile(f):
		try:
			open(f)
			return True
		except IOError as e: return False
	else: return False

def list_existing_files(files):
	"""list_existing_files(sequence of paths) -> list of paths

	Returns a single, flat list of all paths that openable files
	or None if all files are missing"""
	
	result = []
	for f in flatten(files):
		if file_exists(f): result.append(f)
	return None if len(result) == 0 else result

def list_missing_files(files):
	"""list_missing_files(sequence of paths) -> list of paths

	Returns a single, flat list of all paths that are not files or cannot be opened
	or None if all files are present and accessible"""
	
	result = []
	for f in flatten(files):
		if not file_exists(f): result.append(f)
	return None if len(result) == 0 else result
	
def print_log(s):
	print '{} {}'.format(strftime("%m-%d-%Y %H:%M:%S", localtime()),s)
	
def get_reads(read):
	#qseq.txt read -> (paired_end_1_read,paired_end_2_read,barcode_read)
	dirname,basename,flowcell,lane,read,tile = split_qseq_path(read)
	return (join_qseq_path((dirname,lane,'1',tile)),
					join_qseq_path((dirname,lane,'3',tile)),
					join_qseq_path((dirname,lane,'2',tile)))

flowcell_re = re.compile(r"/(?P<flowcell>\d{6}_[^_/]+_\d+(?:_[^_/]+)??_[^_/]+)/")
qseq_re = re.compile(r"s_(?P<lane>\d)_(?P<read>\d+)_(?P<tile>\d+)_qseq\.txt$")
def split_qseq_path(path):
	#/data/storage-1-04/archive/illumina/110616_SN860_0065_2011-101_B817EAABXX/Data/Intensities/BaseCalls/QSEQ/s_7_1_2208_qseq.txt
	#/scratch1/tmp/myourshaw/gmd/new/qseq/links/s_1_1_0002_qseq.txt
	#<dirname>/s_<lane>_<read>_<tile>_qseq.txt-> (dirname,basename,flowcell,lane,read,tile)
	dirname,basename = os.path.split(path)
	match = flowcell_re.search(dirname)
	flowcell = dirname.strip('/').replace('/','.') if match == None else match.group('flowcell')
	match = qseq_re.search(basename)
	if match == None:
		print '{0} is not a qseq file'.format(path)
		sys.exit(1)
	return (dirname,basename,flowcell,match.group('lane'),match.group('read'),match.group('tile'))

def join_qseq_path(p):
	#(dirname,lane,read,tile) -> <dirname>/s_<lane>_<read>_<tile>_qseq.txt
	foo = p[0]
	bar = p[1:]
	return os.path.join(p[0],'s_{0}_qseq.txt'.format('_'.join(p[1:])))
	
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
	unbarcode_tmp_dir = unique_dir('unbarcode_tmp', args.output_dir)
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
	if not args.NoUnbarcodeDemux:
		#required files
		qseqs = list_files(args.qseqs)
		required_files = [args.novobarcode , args.barcodes]
		missing_files = list_missing_files(required_files)
		if len(qseqs) == 0 or missing_files != None:
			if len(qseqs) == 0:
				print 'no qseq files'
			else:
				print FILE_MISSING_MSG.format(args.ini, '\n\t'.join(missing_files))
			sys.exit(1)
	
		#validate existence of reads for both ends and barcode
		for qseq in qseqs:
			dirname,basename,flowcell,lane,read,tile = split_qseq_path(qseq)
			#read 1 must be in list, ignore reads 2 and 3 in the list
			if read != '1':
				continue
			paired_end_1_read,paired_end_2_read,barcode_read = get_reads(qseq)
			if not file_exists(paired_end_2_read):
				print 'paired end 2 read file {0} not found'.format(paired_end_2_read)
				sys.exit(1)
			if not file_exists(barcode_read):
				print 'barcode read file {0} not found'.format(barcode_read)
				sys.exit(1)
	
		#print run info
		print_log("""Barcode file: {1}
	Merged demux'd qseqs will be saved in: {2}/<flowcell>/<lane>/<barcode>
	STDOUT and STDERR will be saved in: {3}"""\
	.format(args.barcodes, unbarcode_merged_dir, qout_dir))
		print_log('Email notifications will {0}.'.format(\
			'not be sent' if args.NoEmail else 'be sent to: {0}'.format(args.email)))
	
		#create job array
		job_array_commands = []
		for qseq in qseqs:
			paired_end_1_read,paired_end_2_read,barcode_read = get_reads(qseq)
			dirname,basename,flowcell,lane,read,tile = split_qseq_path(paired_end_1_read)
			if read != '1': continue
			if args.flowcell != None: flowcell = args.flowcell
			novobarcode_tmp_dir = os.path.join(unbarcode_tmp_dir,flowcell,lane)
			makedir(novobarcode_tmp_dir)
			cmd='\
	{0} \
	-b {1} \
	-d {2} \
	-F QSEQ \
	-f {3} {4} \
	-i {5} \
	--ILQ_SKIP \
	--QSEQ_OUT \
	{6} \
	;'.format(args.novobarcode,args.barcodes,novobarcode_tmp_dir,paired_end_1_read,paired_end_2_read,barcode_read,args.NC_OFF)
			job_array_commands.append(cmd)
		
		#submit cluster job array for all qseq reads
		#command file for job array
		command_file, command_filename = unique_file(os.path.join(unbarcode_tmp_dir,'unbarcode.commands'))
		command_file.writelines("%s\n" % item for item in job_array_commands)
		command_file.close()
		
		normal_exit = True
		jobname = args.jobname+'.demux'
		try:
			s = drmaa.Session()
			s.initialize()
			jt = s.createJobTemplate()
			jt.workingDirectory = unbarcode_tmp_dir
			jt.remoteCommand = drmaa_job_array_runner
			jt.email = [args.email]
			jt.blockEmail = args.NoEmail
			jt.outputPath = ':'+qout_dir
			jt.errorPath = jt.outputPath
			jt.nativeSpecification = '{0} -V'.format(allq)
			jt.args = [command_filename]
			jt.jobName = jobname
			joblist = s.runBulkJobs(jt,1,len(job_array_commands),1)
			print_log('{} job array submitted; commands are in: {}\n\
	unbarcode job ids: {}'.format(args.jobname,command_filename,joblist))
			#wait for all jobs to complete, keeping accounting info
			s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
			#after completion, get status and release accounting info
			print_log('jobId\twasAborted\thasExited\texitStatus\thasSignal\tterminatedSignal\thasCoreDump\tresourceUsage') #\tjobName\tuser\tqueue\thost\tstartTime\tendTime\tuserTime\tsystemTime\twallclockTime\tcpu\tmaxVmem'
			for curjob in joblist:
				status = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
				normal_exit = normal_exit and not status.wasAborted and status.hasExited and not status.exitStatus
				print_log('{jobId}\t{wasAborted}\t{hasExited}\t{exitStatus}\t{hasSignal}\t{terminatedSignal}\t{hasCoreDump}\t{resourceUsage}'\
				.format(**status._asdict()))
		except:
			print_log( "Unexpected error in : {}".format(jobname,sys.exc_info()[0]))
			raise
		else:
			print_log('{} to {} exited {}'.format(jobname,unbarcode_merged_dir,'normally' if normal_exit else 'with fatal errors'))
			if not normal_exit:
				raise MyError('a job in {} had a fatal error'.format(jobname) )
		finally:
			s.deleteJobTemplate(jt)
			s.exit()
	#end of demux
	
	#unbarcode_tmp_dir/<experiment run>/<lane>/<barcode>/<qseq files> -> unbarcode_merged_dir/<<experiment run>.<lane>.<barcode>.pe[12]_qseq.txt
	if not args.NoUnbarcodeMerge:
		print_log('merging unbarcoded qseq files by run.lane.barcode into {}'.format(unbarcode_merged_dir))
		#gather all tile qseq files that were placed into <unbarcode temp dir>/<experiment_run>/<lane>/<barcode>
		#and merge by <experiment_run>.<lane>.<barcode>.<paired end>
		qseq_pe1_re = re.compile(r"^s_[1-8]_[1]_\d+_qseq\.txt$")
		qseq_pe1_az_re = re.compile(r"(?P<a>.*s_[1-8]_)[1](?P<z>_\d+_qseq\.txt)$")
		qseq_pe2_re = re.compile(r"^(?P<a>s_[1-8]_)[3](?P<z>_\d+_qseq\.txt)$")
		job_array_commands = []
		for root, dirs, files in os.walk(unbarcode_unmerged_dir):
			pe1 = []
			pe2 = []
			r,barcode = os.path.split(root)
			r,lane = os.path.split(r)
			r,exptdir = os.path.split(r)
			[pe1.append(os.path.join(root,f)) for f in files if qseq_pe1_re.match(f)]
			if pe1:
				pe1 = sorted(pe1)
				[pe2.append(os.path.join(root,qseq_pe1_az_re.match(p).group('a')+'3'+qseq_pe1_az_re.match(p).group('z'))) for p in pe1]
				job_array_commands.append('cat {} > {}.pe1_qseq.txt'.format(' '.join(pe1), os.path.join(unbarcode_merged_dir,'.'.join((exptdir,lane,barcode)))))
				job_array_commands.append('cat {} > {}.pe2_qseq.txt'.format(' '.join(pe2), os.path.join(unbarcode_merged_dir,'.'.join((exptdir,lane,barcode)))))
		
		command_file, command_filename = unique_file(os.path.join(unbarcode_tmp_dir,'merge.commands'))
		command_file.writelines("%s\n" % item for item in job_array_commands)
		command_file.close()
		
		#run cat commands on cluster
		normal_exit = True
		jobname = args.jobname+'.merge'
		try:
			s = drmaa.Session()
			s.initialize()
			jt = s.createJobTemplate()
			jt.workingDirectory = unbarcode_tmp_dir
			jt.remoteCommand = drmaa_job_array_runner
			jt.email = [args.email]
			jt.blockEmail = args.NoEmail
			jt.outputPath = ':'+qout_dir
			jt.errorPath = jt.outputPath
			jt.nativeSpecification = '{0} -V'.format(allq)
			jt.args = [command_filename]
			jt.jobName = jobname
			joblist = s.runBulkJobs(jt,1,len(job_array_commands),1)
			print_log('{} job array submitted; commands are in: {}\n\
		job ids: {}'.format(jobname,command_filename,joblist))
			#wait for all jobs to complete, keeping accounting info
			s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
			#after completion, get status and release accounting info
			print_log('jobId\twasAborted\thasExited\texitStatus\thasSignal\tterminatedSignal\thasCoreDump\tresourceUsage') #\tjobName\tuser\tqueue\thost\tstartTime\tendTime\tuserTime\tsystemTime\twallclockTime\tcpu\tmaxVmem'
			for j in joblist:
				status = s.wait(j, drmaa.Session.TIMEOUT_WAIT_FOREVER)
				normal_exit = normal_exit and not status.wasAborted and status.hasExited and not status.exitStatus
				print_log('{jobId}\t{wasAborted}\t{hasExited}\t{exitStatus}\t{hasSignal}\t{terminatedSignal}\t{hasCoreDump}\t{resourceUsage}'\
				.format(**status._asdict()))
		except:
			print_log( "Unexpected error in : {}".format(jobname,sys.exc_info()[0]))
			raise
		else:
			print_log('{} to {} exited {}'.format(jobname,unbarcode_merged_dir,'normally' if normal_exit else 'with fatal errors'))
			if not normal_exit:
				raise MyError('a job in {} had a fatal error'.format(jobname) )
		finally:
			s.deleteJobTemplate(jt)
			s.exit()

if __name__ == "__main__": sys.exit(main())
