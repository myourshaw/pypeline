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

class Vcf2FilteredVcfError(Exception): pass

#--top_dir /scratch0/tmp/myourshaw/gmd --parent_jobName vcfs2filteredvcfs_20110725202759 --parent_job_dir /scratch0/tmp/myourshaw/gmd/job_info/tmp_jyb2h7/vcfs2filteredvcfs_20110725202759_Zwp_2H --vcf /scratch0/tmp/myourshaw/gmd/analysis/vcf/GMD108.GMD108_XT_8-2.gatk.realigned.recalibrated.vcf --indels_vcf /scratch0/tmp/myourshaw/gmd/analysis/vcf/GMD108.GMD108_XT_8-2.gatk.realigned.recalibrated_indels.vcf --email myourshaw@ucla.edu

def main():

	#command line arguments
	parser = argparse.ArgumentParser(parents=[my.get_directory_parser()],
		description = 'GATK ',
		epilog = 'pypeline.vcf2filteredvcf version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--vcf', required=True,
											help='input vcf file')
	parser.add_argument('--indels_vcf', required=True,
											help='input raw indels vcf mask file from UnifiedGenotyper')
	parser.add_argument('--gatk',
											help='path to GenomeAnalysisTK jar')
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
	args = parser.parse_args()
	my.set_directory_defaults(args)
	config = my.get_config(args)
	
	vcfdir = os.path.dirname(args.vcf)
	vcfbase = os.path.basename(args.vcf).rstrip('.vcf')
	vcf_tmp_dir = my.unique_dir('tmp_'+vcfbase,vcfdir)

	if not args.gatk: args.gatk = config.get('gatk','GenomeAnalysisTK')
	if not my.file_exists(args.gatk):
		raise Vcf2FilteredVcfError('cannot find GenomeAnalysisTK jar')
	if not args.ref: args.ref = config.get('gatk','human_g1k_v37')
	if not my.file_exists(args.ref):
		raise Vcf2FilteredVcfError('cannot find reference genome')
	if not args.dbsnp: args.dbsnp = config.get('gatk','dbsnp')
	if not my.file_exists(args.dbsnp):
		raise Vcf2FilteredVcfError('cannot find dbsnp')
	if not args.hapmap: args.hapmap = config.get('gatk','hapmap')
	if not my.file_exists(args.hapmap):
		raise Vcf2FilteredVcfError('cannot find hapmap')
	if not args.onekg_indels: args.onekg_indels = config.get('gatk','onekg_indels')
	if not my.file_exists(args.onekg_indels):
		raise Vcf2FilteredVcfError('cannot find 1000 genomes indels')
	if not args.onekg_omni: args.onekg_omni = config.get('gatk','onekg_omni')
	if not my.file_exists(args.onekg_omni):
		raise Vcf2FilteredVcfError('cannot find 1000 genomes omni')
	if not args.rscript: args.rscript = config.get('DEFAULT','Rscript')
	if not my.file_exists(args.rscript):
		raise Vcf2FilteredVcfError('cannot find Rscript')
	if not args.gatkrscripts: args.gatkrscripts = config.get('gatk','gatk_R_scripts')
	if not os.path.isdir(args.gatkrscripts):
		raise Vcf2FilteredVcfError('cannot find gatk R script resources directory')

	snpFiltered_indelFiltered_vcf = os.path.join(args.vcf_dir, vcfbase+'.snpFiltered.indelFiltered.vcf')
	recal = os.path.join(vcf_tmp_dir, vcfbase+'.VariantRecalibrator.recal')
	tranches = os.path.join(vcf_tmp_dir, vcfbase+'.VariantRecalibrator.tranches')
	plots = os.path.join(vcf_tmp_dir, vcfbase+'.VariantRecalibrator.plots')
	recalibrated_vcf = os.path.join(args.vcf_dir, snpFiltered_indelFiltered_vcf.rstrip('.vcf')+'.recalibrated.vcf')

	this_job_name = 'vcf2filteredvcf_{}_{}'.format(os.path.basename(args.vcf), my.localtime_squish())
	this_job_dir = my.unique_dir(this_job_name, args.parent_job_dir)
	
	VariantFiltration_jobName = '{}_variantfiltration'.format(this_job_name)
	VariantRecalibrator_jobName = '{}_variantrecalibrator'.format(this_job_name)
	ApplyRecalibration_jobName = '{}_applyrecalibration'.format(this_job_name)
	
	
	#VariantFiltration
	cmd="""java -Xmx5g -Djava.io.tmpdir="{}" \
-jar "{}" \
-T VariantFiltration \
-l INFO \
-R "{}" \
-B:variant,VCF "{}" \
-o "{}" \
--clusterWindowSize 10 \
--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
--filterName "HARD_TO_VALIDATE" \
--filterExpression "SB >= -1.0" \
--filterName "StrandBiasFilter" \
--filterExpression "QUAL < 10" \
--filterName "QualFilter" \
--filterExpression "QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10" \
--filterName GATKStandard \
-B:mask,VCF "{}" \
--maskName InDel \
;""".format(args.tmp_dir, args.gatk, args.ref, args.vcf, snpFiltered_indelFiltered_vcf, args.indels_vcf)
	
	jobName = VariantFiltration_jobName
	job_dir = this_job_dir
	my.print_log('job {} started; command:\n{}'.format(jobName, cmd))
	j = job.Job(jobName = jobName,
	outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = job_dir,
	email = args.email, blockEmail = not args.sendemailforintermediatejobs,
	processors = 8, queuingModeBunchJobsOnNodes=False)
	try:
		j.executeCommand(cmd, synchronous=True)
	except Exception as e:
		my.print_err('There was an error running job {}\n{}\nintermediate files are in {}\n'.format(jobName, e, job_dir))
		raise e
	else:
		with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
			s.write(format('\n'.join(j.completionMsg)))

	#VariantRecalibrator
	cmd="""java -Xmx5g -Djava.io.tmpdir="{}" \
-jar "{}" \
-T VariantRecalibrator \
-l INFO \
-R "{}" \
-B:input,VCF "{}" \
-B:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 "{}" \
-B:omni,VCF,known=false,training=true,truth=false,prior=12.0 "{}" \
-B:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 "{}" \
-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
-recalFile "{}" \
-tranchesFile "{}" \
-rscriptFile "{}" \
-Rscript "{}" \
-resources "{}" \
--target_titv 3.2 \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
--ignore_filter HARD_TO_VALIDATE \
--ignore_filter LowQual \
--percentBadVariants 0.05 \
--maxGaussians 4 \
;""".format(args.tmp_dir, args.gatk, args.ref, snpFiltered_indelFiltered_vcf, args.hapmap, args.onekg_omni, args.dbsnp, recal, tranches, plots, args.rscript, args.gatkrscripts)
	
	jobName = VariantRecalibrator_jobName
	job_dir = this_job_dir
	my.print_log('job {} started; command:\n{}'.format(jobName, cmd))
	j = job.Job(jobName = jobName,
	outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = job_dir,
	email = args.email, blockEmail = not args.sendemailforintermediatejobs,
	processors = 8, queuingModeBunchJobsOnNodes=False)
	try:
		j.executeCommand(cmd, synchronous=True)
	except Exception as e:
		my.print_err('There was an error running job {}\n{}\nintermediate files are in {}\n'.format(jobName, e, job_dir))
		raise e
	else:
		with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
			s.write(format('\n'.join(j.completionMsg)))

	#ApplyRecalibration
	cmd="""java -Xmx5g -Djava.io.tmpdir="{}" \
-jar "{}" \
-T ApplyRecalibration \
-l INFO \
-R "{}" \
-B:input,VCF "{}" \
--ts_filter_level 99.0 \
-tranchesFile "{}" \
-recalFile "{}" \
-o "{}" \
;""".format(args.tmp_dir, args.gatk, args.ref, snpFiltered_indelFiltered_vcf, tranches, recal, recalibrated_vcf)
	
	jobName = ApplyRecalibration_jobName
	job_dir = this_job_dir
	my.print_log('job {} started; command:\n{}'.format(jobName, cmd))
	j = job.Job(jobName = jobName,
	outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = job_dir,
	email = args.email, blockEmail = not args.sendemailforintermediatejobs,
	processors = 8, queuingModeBunchJobsOnNodes=False)
	try:
		j.executeCommand(cmd, synchronous=True)
	except Exception as e:
		my.print_err('There was an error running job {}\n{}\nintermediate files are in {}\n'.format(jobName, e, job_dir))
		raise e
	else:
		with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
			s.write(format('\n'.join(j.completionMsg)))

if __name__ == "__main__": sys.exit(main())
