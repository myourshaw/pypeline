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
import bam2vcf

class Bams2VcfsError(Exception): pass

#-d /scratch0/tmp/myourshaw/gmd --bams /scratch0/tmp/myourshaw/gmd/reads/library_bams/GMD159.GMD159_XT_19.markdup.bam --cohort gmd --email myourshaw@ucla.edu --sendemailforintermediatejobs

def main():

	#command line arguments
	parser = argparse.ArgumentParser(parents=[my.get_directory_parser()],
		description = 'Starting with markdup bam files (that may be specified in metadata), realign indels, recalibrate base quality scores, genotype, and filter/recalibrate variants; optionally merge cohort(s)',
		epilog = 'pypeline.bams2vcfs version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--metadata', '-m',
											help='metadata file (ignored if --bams present)',)
	parser.add_argument('--readgroups', nargs='*',
											help='readgroup_ids from metadata file to process, "" selects blank readgroup_ids (default: all) (ignored if --bams present)',)
	parser.add_argument('--cohorts', nargs='*',
											help='list of cohorts to process, "" selects blank cohort (default: all) (ignored if --bams present)',)
	parser.add_argument('--bams', nargs='+',
											help='bam file list',)
	parser.add_argument('--cohort', default='',
											help='cohort name for bams (default: cohort)',)
	parser.add_argument('--merge_cohorts', action='store_true', default=False,
											help='merge bam files by cohort into combined vcf')
	parser.add_argument('--targets',
											help='targeted intervals file')
	parser.add_argument('--gatk',
											help='path to GenomeAnalysisTK jar')
	parser.add_argument('--analyzecovariates',
											help='path to AnalyzeCovariates jar')
	parser.add_argument('--ref',
											help='path to reference genome fasta file')
	parser.add_argument('--dbsnp',
											help='path to dbsnp vcf')
	parser.add_argument('--hapmap',
											help='path to hapmap vcf')
	parser.add_argument('--onekg_indels',
											help='path to 1000 genomes indels vcf')
	parser.add_argument('--onekg_omni',
											help='path to 1000 genomes omni vcf')
	parser.add_argument('--rscript',
											help='path to Rscript executable')
	parser.add_argument('--gatkrscripts',
											help='path to gatk R scripts resources directory')
	parser.add_argument('--picard_ValidateSamFile',
											help='path to picard ValidateSamFile jar')
	args = parser.parse_args()
	my.set_directory_defaults(args)
	config = my.get_config(args)
		
	if not args.parent_job_dir: args.parent_job_dir = args.job_info_dir
	this_jobName = 'bams2vcfs_{}'.format(my.localtime_squish())
	this_job_dir = my.unique_dir(this_jobName, args.parent_job_dir)

	bams = set()
	if args.bams:
		for bam in args.bams:
			[bams.add(b) for b in glob.glob(bam)]
	elif args.metadata:
		[bams.add(m['library_markdup_bam']) for m in my.get_metadata(args.metadata) if (not args.readgroups or m['readgroup_id'] in args.readgroups) and (not args.cohorts or m['cohort'] in args.cohorts)]
	if not bams:
		raise Bams2VcfsError('no bams found')
	bams = list(bams)

	cmds = []
	
	if args.merge_cohorts:
		bams = set()
		if args.bams:
			for bam in args.bams:
				[bams.add(b) for b in glob.glob(bam)]
			bam_list = ' '.join(bams)
			cohort_bam_list = (args.cohort, bam_list)
		elif args.metadata:
			cohort_bam = [(m['cohort'],m['library_markdup_bam']) for m in my.get_metadata(args.metadata) if (not args.readgroups or m['readgroup_id'] in args.readgroups) and (not args.cohorts or m['cohort'] in args.cohorts)]
			bam_list = [m for m in my.get_metadata(args.metadata)]
			
	for bam in bams:
		cmds += ['python {} --top_dir {} --parent_jobName {} --parent_job_dir {} --bam {} --email {} {} --cohort {}'
						 .format(bam2vcf.__file__, args.top_dir, this_jobName, this_job_dir, bam, ' '.join(args.email), '--sendemailforintermediatejobs' if args.sendemailforintermediatejobs else '', args.cohort)]

	print '\n'.join(cmds)
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
		raise e
	else:
		my.print_log('job {} finished; status:\n{}'.format(jobName, '\n'.join(j.completionMsg)))
		with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
			s.write(format('\n'.join(j.completionMsg)))
			
if __name__ == "__main__": sys.exit(main())
