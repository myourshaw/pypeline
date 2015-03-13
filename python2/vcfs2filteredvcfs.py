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
import vcf2filteredvcf

class Vcfs2FilteredVcfsError(Exception): pass

#-d /scratch0/tmp/myourshaw/gmd --vcfs /scratch0/tmp/myourshaw/gmd/analysis/vcf/GMD108.GMD108_XT_8-2.gatk.realigned.recalibrated.vcf --indels_vcfs /scratch0/tmp/myourshaw/gmd/analysis/vcf/GMD108.GMD108_XT_8-2.gatk.realigned.recalibrated_indels.vcf

def main():

	#command line arguments
	parser = argparse.ArgumentParser(parents=[my.get_directory_parser()],
		description = 'raw vcfs from Unified Genotyper -> VariantFiltration -> VariantRecalibrator -> ApplyRecalibration',
		epilog = 'pypeline.vcfs2filteredvcfs version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--vcfs', nargs='+',
											help='input vcf files')
	parser.add_argument('--indels_vcfs', nargs='+',
											help='input raw indels vcf mask files from UnifiedGenotyper, in same order as --vcfs')
	args = parser.parse_args()
	my.set_directory_defaults(args)
	config = my.get_config(args)
	
	if len(args.vcfs) != len(args.indels_vcfs):
		raise Vcfs2FilteredVcfsError('number of vcf files ({}) is not equal to number of indel mask vcf files ({})'.format(len(args.vcfs), len(args.indels_vcfs)))

	this_jobName = 'vcfs2filteredvcfs_{}'.format(my.localtime_squish())
	this_job_dir = my.unique_dir(this_jobName, args.parent_job_dir)
	
	cmds = []
	for i in range(len(args.vcfs)):
		cmds += ['python {} --top_dir {} --parent_jobName {} --parent_job_dir {} --vcf {} --indels_vcf {} {} {}'
						 .format(vcf2filteredvcf.__file__, args.top_dir, this_jobName, this_job_dir, args.vcfs[i], args.indels_vcfs[i], '--email {}'.format(' '.join(args.email) if args.email else ''), '--sendemailforintermediatejobs' if args.sendemailforintermediatejobs else '')]
	
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
