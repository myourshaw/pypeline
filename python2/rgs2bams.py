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
import my
import job
import rg2bam

class Rg2BamError(Exception): pass

#1/12 lane
#-m /home/myourshaw/lab/git-myourshaw/python/metadata.txt -d /scratch0/tmp/myourshaw/gmd/bam_test -i /scratch0/tmp/myourshaw/gmd/runs_tmp/110511_SN430_0243_B817FLABXX.1.1.TTAGGCA_qseq.txt.gz /scratch0/tmp/myourshaw/gmd/runs_tmp/110511_SN430_0243_B817FLABXX.1.3.TTAGGCA_qseq.txt.gz

#lane
#-m /home/myourshaw/lab/git-myourshaw/python/metadata.txt -d /scratch0/tmp/myourshaw/gmd/bam_test -i /data/storage-1-04/archive/illumina/110616_SN860_0065_2011-101_B817EAABXX/Data/Intensities/BaseCalls/ReadGroups/110616_SN860_0065_2011-101_B817EAABXX.1.1.GATCAGA_qseq.txt.gz	/data/storage-1-04/archive/illumina/110616_SN860_0065_2011-101_B817EAABXX/Data/Intensities/BaseCalls/ReadGroups/110616_SN860_0065_2011-101_B817EAABXX.1.3.GATCAGA_qseq.txt.gz

#metadata
#-m /home/myourshaw/lab/git-myourshaw/python/metadata.txt -d /scratch0/tmp/myourshaw/gmd --email myourshaw@ucla.edu --sendemailforintermediatejobs

#metadata, limited readgroups
#--readgroups 110616_SN860_0065_2011-101_B817EAABXX.6.CAGATCA 110616_SN860_0065_2011-101_B817EAABXX.7.ACTTGAA 110616_SN860_0065_2011-101_B817EAABXX.5.GCCAATA 110616_SN860_0065_2011-101_B817EAABXX.3.GGCTACA -m /home/myourshaw/lab/git-myourshaw/python/metadata.txt -d /data/storage-1-02/archive/myourshaw/gmd --email myourshaw@ucla.edu --sendemailforintermediatejobs

def main():

	#command line arguments
	parser = argparse.ArgumentParser(parents=[my.get_directory_parser()],
		description = 'From a metadata file, align multiple readgroups of qseq readgroup files with novoalign, calculate metrics, remove duplicates with picard, re-calculate metrics, validate bam file',
		epilog = 'pypeline.rgs2bams version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--metadata', '-m', required=True,
											help='metadata file',)
	parser.add_argument('--readgroups', nargs='*',
											help='readgroup_ids from metadata file to process, "" selects blank readgroup_ids (default: all)',)
	parser.add_argument('--novoalign','-n',
											help='path to novoalign executable (default: from config->novocraft.novoalign)')
	parser.add_argument('--index','-x',
											help='path to novoindex-produced index (default: from config->novocraft.default_index)')
	parser.add_argument('--adapters', '-a',
											help='sequencing adapters (default: from config->novocraft.default_adapter_1,default_adapter_2 or novoalign adapter 1 and 2 defaults)')
	args = parser.parse_args()
	my.set_directory_defaults(args)
	config = my.get_config(args)
	
	if not args.parent_job_dir: args.parent_job_dir = args.job_info_dir
	this_jobName = 'rgs2bams_{}_{}'.format(os.path.basename(args.metadata), my.localtime_squish())
	this_job_dir = my.unique_dir(this_jobName, args.parent_job_dir)
	
	cmds = []
	
	if not my.file_exists(args.metadata):
		raise Rg2BamError('metadata file {} does not exist'.format(args.metadata))
	metadata_list = [m for m in my.get_metadata(args.metadata) if not args.readgroups or m.get('readgroup_id','').strip() in args.readgroups]
	if not metadata_list:
		raise Rg2BamError('no records selected from metadata file {}'.format(args.metadata))
	for metadata in metadata_list:
		file1 = metadata.get('readgroup_read_1_qseq', None)
		if not my.file_exists(file1):
			warn('metadata line {} [{}] readgroup_read_1_qseq missing. ignoring this record\n'.format(metadata.get('line_number','?'), metadata))
			continue
		file2 = metadata.get('readgroup_read_2_qseq', '')
		if not my.file_exists(file2):
			warn('metadata line {} [{}] readgroup_read_2_qseq missing; aligning as single-end reads\n'.format(metadata.get('line_number','?'), metadata))
		if file2 and not my.qseq_read1_matches(file1, file2):
			warn('metadata line {} [{}] first reads in paired end qseq input files do not match. ignoring this record\n:\n{}\n{}\n'.format(metadata.get('line_number','?'), metadata, my.qseq_peek(file1), my.qseq_peek(file2)))
			continue
		qseq_info = my.qseq_peek(file1)
		sample = metadata.get('sample', None)
		if not sample:
			raise Rg2BamError('no sample specified in metadata file')
		library = metadata.get('library', sample)
		experiment_run_id = metadata.get('experiment_run_id', None)
		if not experiment_run_id:
			experiment_run_id = '{}_{}'.format(qseq_info['machine'], qseq_info['run'])
		machine_input = metadata.get('machine', None)
		machine = str(qseq_info['machine'])
		run_number_input = metadata.get('run_number', None)
		run_number = str(qseq_info['run'])
		lane_input = metadata.get('lane', None)
		lane = str(qseq_info['lane'])
		if (machine_input and machine_input != machine) or (run_number_input and run_number_input != run_number) or(lane_input and lane_input != lane):
			raise Rg2BamError('metadata line {}, specified machine: {}, run: {}, lane: {} does not match qseq file machine: {}, run: {}, lane: {}'.format(metadata.get('line_number','?'), machine_input, run_number_input, lane_input, machine, run_number, lane))
		barcode = metadata.get('barcode', '')
		if not barcode and re.search(r"[_.][ACGT]+[_\.]", file1, re.I):
			warn('metadata line {}, no barcode specified but filename looks like it contains a barcode'.format(metadata.get('line_number','?'), ))
		sequencing_center = metadata.get('sequencing_center', 'UCLA')
		run_date = metadata.get('run_date', None)
		if not run_date:
			run_date = strftime('%Y-%m-%d', localtime(os.path.getctime()))
		predicted_median_insert_size = metadata.get('predicted_median_insert_size', None)
		platform = metadata.get('platform', None)
		platform_unit = metadata.get('platform_unit', None)
		if not platform_unit:
			platform_unit = '{}.{}.{}{}'.format(qseq_info['machine'], qseq_info['run'], qseq_info['lane'], '.' + barcode if barcode else '')
		adapters = []
		if not adapters:
			adapters = [metadata.get('adapter_1', ''), metadata.get('adapter_2', '')]
		if not adapters[0]:
			adapters[0] = config.get('novocraft','default_adapter_1')
		if not adapters[1]:
			config.get('novocraft','default_adapter_2')
		readgroup_id = metadata.get('readgroup_id', None)
		if not readgroup_id:
			readgroup_id = '{}.{}{}'.format(experiment_run_id, qseq_info['lane'], '.'+barcode if barcode else '')
		readgroup_description = metadata.get('readgroup_description', None)
		if not readgroup_description:
			 readgroup_description = 'sample={};library={};sequencing_center={};platform={};platform_unit={};run_date={};experiment_run_id={};lane={};barcode={};predicted_median_insert_size={}'\
			.format(sample, library, sequencing_center, platform, platform_unit, run_date, experiment_run_id, lane, barcode, predicted_median_insert_size)
		RG = '"@RG\\tID:{}\\tCN:{}\\tDT:{}\\tLB:{}\\tPI:{}\\tPL:{}\\tPU:{}\\tSM:{}\\tDS:{}"'.\
			format(readgroup_id, sequencing_center, run_date, library, predicted_median_insert_size, platform, platform_unit, sample, readgroup_description)
		
		output_name = '{}.{}.{}'.format(sample, library, readgroup_id)
		novoalign_output_sam = os.path.join(args.readgroup_bams_dir, '{}.novoalign.sam'.format(output_name))
		readgroup_output_bam = os.path.join(args.readgroup_bams_dir, '{}.novoalign.bam'.format(output_name))
		
		cmds += ['python {} --top_dir {} --parent_jobName {} --parent_job_dir {} --qseqs {} {} --novoalign_output_sam {} --readgroup_output_bam {} --RG {} --adapters {}'
						 .format(rg2bam.__file__, args.top_dir, this_jobName, this_job_dir, file1, file2 if file2 else '', novoalign_output_sam, readgroup_output_bam, RG, ' '.join(adapters))]

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
