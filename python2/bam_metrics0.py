#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import ConfigParser #configparser in python 3
import glob
import tempfile
from time import localtime, strftime
import my
import job

class BamMetricsError(Exception): pass

#--NoDepthOfCoverage --output_prefix rg_pre-rmdup_ -d /data/storage-1-02/archive/myourshaw/gmd -i /data/storage-1-02/archive/myourshaw/gmd/reads/readgroup_bams/*.bam
# --NoDepthOfCoverage --output_prefix 

def main():
	
	#command line arguments
	parser = argparse.ArgumentParser(parents=[my.get_directory_parser()],
		description = 'Calculate bam file metrics using picard and gatk DepthOfCoverage',
		epilog = 'pypeline.bam_metrics version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--input', '-i', metavar='BAM', nargs='+',
											help='input list of bam files')
	parser.add_argument('--output_prefix',
											help='prefix for output file names')
	parser.add_argument('--ref','-R', nargs='?',
											help='reference genome fasta (default: config->b37_fasta)')
	parser.add_argument('--genelist','-G', nargs='?',
											help='GATK DepthOfCoverage geneList (default: config->refgene_b37_sorted)')
	parser.add_argument('--targetlimit','-L', nargs='?',
											help='GATK DepthOfCoverage target limit (default: bait interval list)')
	parser.add_argument('--bait','-B', metavar='BAIT_INTERVAL_LIST', nargs='?',
											help='bait interval list (default: config->sureselect_50mb_interval_list)')
	parser.add_argument('--targets','-T', metavar='TARGET_INTERVAL_LIST', nargs='+',
											help='list of target interval_list files (default: <bait interval list> config->ensembl_all_genes config->ensembl_protein_coding_genes)')
	parser.add_argument('--regex','-X', nargs='?',default=r'[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*',
											help='read name regex (default for illumina/novoalign: "[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*")')
	parser.add_argument('--pixels','-D', nargs='?', type=int, default=100,
											help='optical duplicate pixel distance (default for illumina: 100)')
	parser.add_argument('--consolidate', metavar='METRIC', nargs='*',
											help='additional files to include in consolidation (e.g., *hahahaaha.markdup.metrics')
	parser.add_argument('--NoPicardMetrics', action='store_true', default=False, help='do not calculate picard metrics')
	parser.add_argument('--Noidxstats', action='store_true', default=False, help='do not calculate idxstats')
	parser.add_argument('--NoCalculateHsMetrics', action='store_true', default=False, help='do not calculate HsMetrics')
	parser.add_argument('--NoBamIndexStats', action='store_true', default=False, help='do not calculate BamIndexStats')
	parser.add_argument('--NoCollectAlignmentSummaryMetrics', action='store_true', default=False, help='do not calculate AlignmentSummaryMetrics')
	parser.add_argument('--NoCollectInsertSizeMetrics', action='store_true', default=False, help='do not calculate InsertSizeMetrics')
	parser.add_argument('--NoQualityScoreDistribution', action='store_true', default=False, help='do not calculate QualityScoreDistribution')
	parser.add_argument('--NoMeanQualityByCycle', action='store_true', default=False, help='do not calculate MeanQualityByCycle')
	parser.add_argument('--NoCollectGcBiasMetrics', action='store_true', default=False, help='do not calculate GcBiasMetrics')
	parser.add_argument('--NoEstimateLibraryComplexity', action='store_true', default=False, help='do not calculate LibraryComplexity')
	parser.add_argument('--NoDepthOfCoverage', action='store_true', default=False, help='do not calculate DepthOfCoverage')
	args = parser.parse_args()
	my.set_directory_defaults(args)
	config = my.get_config(args)

	if not args.parent_job_dir: args.parent_job_dir = args.job_info_dir
	this_jobName = 'bammetrics_{}'.format(my.localtime_squish())
	this_job_dir = my.unique_dir(this_jobName, args.parent_job_dir)
	
	picard_jobName = 'picardMetrics_{}'.format(my.localtime_squish())
	picard_job_dir = my.unique_dir(picard_jobName, this_job_dir)
	doc_jobName = 'DepthOfCoverage_{}'.format(my.localtime_squish())
	doc_job_dir = my.unique_dir(doc_jobName, this_job_dir)

	#get defaults from configuration file
	if args.ref == None: args.ref = config.get('reference','b37_fasta')
	if args.genelist == None: args.genelist = config.get('gene','refgene_b37_sorted')
	if args.bait == None: args.bait = config.get('interval','sureselect_50mb_interval_list')
	if args.targetlimit == None: args.targetlimit = args.bait
	if args.targets == None: args.targets = [args.bait, config.get('interval','ensembl_all_genes'),config.get('interval','ensembl_protein_coding_genes')]
	if args.output_prefix == None: args.output_prefix = ''
	
	GenomeAnalysisTK = config.get('gatk','GenomeAnalysisTK')
	CalculateHsMetrics = config.get('picard','CalculateHsMetrics')
	BamIndexStats = config.get('picard','BamIndexStats')
	CollectMultipleMetrics = config.get('picard','CollectMultipleMetrics')
	CollectGcBiasMetrics = config.get('picard','CollectGcBiasMetrics')
	EstimateLibraryComplexity = config.get('picard','EstimateLibraryComplexity')
	MetricsConsolidate = config.get('picard','MetricsConsolidate')
	samtools = config.get('samtools','samtools')
	
	#constants
	CONSOLIDATED_DIR = 'consolidated_metrics'
	COVERAGE_DIR = 'DepthOfCoverage'
	
	FILE_MISSING_MSG = 'These required files are missing. You may need to edit {0}.\n\t{1}'
	RUN_JOB_MSG = """Job: {0} [{1}]
Command: {2}
STDOUT and STDERR in {3}.
Notification to {4}."""
	EXCEPTION_MSG = 'ERROR: {0}; {1}'
	
	PICARD_CMD = '$(which java) -Xmx5g -Djava.io.tmpdir='+args.tmp_dir+ \
' -jar {0} \
TMP_DIR='+args.tmp_dir+ \
' VERBOSITY=INFO \
QUIET=false \
VALIDATION_STRINGENCY=SILENT \
COMPRESSION_LEVEL=5 \
MAX_RECORDS_IN_RAM=1000000 \
INPUT={1}'
	
	PICARD_CMD_BAMINDEXSTATS = PICARD_CMD + """ > {2}"""
	PICARD_CMD_OUTPUT = PICARD_CMD + """ OUTPUT={2}"""
	PICARD_CMD_HSMETRICS = PICARD_CMD_OUTPUT + """ BAIT_INTERVALS={3} TARGET_INTERVALS={4}"""
	PICARD_CMD_MULTIPLEMETRICS = PICARD_CMD_OUTPUT + """ REFERENCE_SEQUENCE={3} {4} {5} {6} {7}"""
	PICARD_CMD_GCBIAS = PICARD_CMD_OUTPUT + """ REFERENCE_SEQUENCE={3} CHART_OUTPUT={4} SUMMARY_OUTPUT={5}"""
	PICARD_CMD_LIBRARYCOMPLEXITY = PICARD_CMD_OUTPUT + """ READ_NAME_REGEX={3} OPTICAL_DUPLICATE_PIXEL_DISTANCE={4}"""

	#required files
	bams = my.list_files(args.input)
	required_files = [args.ref, args.genelist, args.targetlimit, args.bait, args.targets, bams,
		GenomeAnalysisTK, CalculateHsMetrics, BamIndexStats, CollectMultipleMetrics,
		CollectGcBiasMetrics, EstimateLibraryComplexity, MetricsConsolidate, samtools]
	missing_files = my.list_missing_files(required_files)
	if len(bams) == 0 or missing_files != None:
		raise BamMetricsError('These required files are missing\n\t{}'.format('\n\t'.join(missing_files)))
	
	#output directories
	consolidated_dir = os.path.join(args.bam_metrics_dir, args.output_prefix+CONSOLIDATED_DIR)
	my.makedir(consolidated_dir)
	coverage_dir = os.path.join(args.bam_metrics_dir, args.output_prefix+COVERAGE_DIR)
	my.makedir(coverage_dir)

	#print run info
	print """Bam Metrics run {0}
Reference genome: {1}
Depth of coverage gene list: {2}
Depth of coverage limit: {3}
Bait: {4}
Target(s): {5}
Bam file(s): {6}
Metrics will be saved in: {7}
Consolidated metrics will be saved in: {8}
STDOUT and STDERR will be saved in: {9}"""\
.format(strftime("%a, %d %b %Y %H:%M:%S +0000", localtime()), args.ref, args.genelist, args.targetlimit, args.bait, ';'.join(args.targets), ';'.join(bams), args.bam_metrics_dir, consolidated_dir, args.qout_dir)

	cmds = []
	doc_cmds = []
	files_to_consolidate = []
	doc_files_to_consolidate = []

	#get metrics for each bam file
	#try:
	for bam in bams:
		#output files
		metricsbase = os.path.join(args.bam_metrics_dir, args.output_prefix+os.path.basename(bam)+'.metrics')
		validate = metricsbase + '.validate'
		idxstats = metricsbase + '.idxstats'
		indexstats = metricsbase + '.indexstats'
		gcbiasmetrics = metricsbase + '.gcbias.table'
		gcbiasmetricschart = metricsbase + '.gcbias.pdf'
		gcbiasmetricssummary = metricsbase + '.gcbias.summary'
		librarycomplexity = metricsbase + '.librarycomplexity'
		coverage_output_prefix = os.path.join(coverage_dir, args.output_prefix+os.path.basename(bam))
		

		if not args.NoPicardMetrics:
			#idxstats
			if not args.Noidxstats:
				cmds.append('{0} idxstats {1} > {2};'.format (samtools, bam, idxstats))
			
			#BamIndexStats
			if not args.NoBamIndexStats:
				cmds.append(PICARD_CMD_BAMINDEXSTATS.format(BamIndexStats, bam, indexstats))
	
			#CalculateHsMetrics
			if not args.NoCalculateHsMetrics:
				for target in args.targets:
					hsoutput = '{0}.target_{1}.hsmetrics'.format(metricsbase,os.path.basename(target))
					cmds.append(PICARD_CMD_HSMETRICS.format(CalculateHsMetrics, bam, hsoutput, args.bait, target))
					files_to_consolidate.append(hsoutput)
		
			#CollectMultipleMetrics
			if not args.NoCollectAlignmentSummaryMetrics or not args.NoCollectInsertSizeMetrics \
			or not args.NoQualityScoreDistribution or not args.NoMeanQualityByCycle:
				cmds.append(PICARD_CMD_MULTIPLEMETRICS.format(CollectMultipleMetrics, bam, metricsbase, args.ref, \
					'' if args.NoCollectAlignmentSummaryMetrics else 'PROGRAM=CollectAlignmentSummaryMetrics', \
					'' if args.NoCollectInsertSizeMetrics else 'PROGRAM=CollectInsertSizeMetrics', \
					'' if args.NoQualityScoreDistribution else 'PROGRAM=QualityScoreDistribution', \
					'' if args.NoMeanQualityByCycle else 'PROGRAM=MeanQualityByCycle'))
				files_to_consolidate.append(metricsbase+'.*_metrics')
	
			#CollectGcBiasMetrics
			if not args.NoCollectGcBiasMetrics:
				cmds.append(PICARD_CMD_GCBIAS.format(CollectGcBiasMetrics, bam, gcbiasmetrics, args.ref, gcbiasmetricschart, gcbiasmetricssummary))
				files_to_consolidate.append(gcbiasmetrics)
				files_to_consolidate.append(gcbiasmetricssummary)
	
			#EstimateLibraryComplexity
			if not args.NoEstimateLibraryComplexity:
				cmds.append(PICARD_CMD_LIBRARYCOMPLEXITY.format(EstimateLibraryComplexity, bam, librarycomplexity, args.regex, args.pixels))
				files_to_consolidate.append(librarycomplexity)

		#DepthOfCoverage
		if not args.NoDepthOfCoverage:
			doc_cmd='\
$(which java) -Xmx16g -Djava.io.tmpdir={0} \
-jar {1} \
-T DepthOfCoverage \
-l INFO \
-R {2} \
-I {3} \
-o {4} \
-L {5} \
-geneList {6} \
--omitDepthOutputAtEachBase \
-dcov 1000 \
-pt sample \
-ct 1 -ct 2 -ct 3 -ct 4 -ct 5 -ct 6 -ct 7 -ct 8 -ct 9 -ct 10 -ct 20 -ct 30 -ct 40 -ct 50 -ct 60 -ct 70 -ct 80 -ct 90 -ct 100 -ct 110 -ct 120 -ct 130 -ct 140 -ct 150 -ct 160 -ct 170 -ct 180 -ct 190 -ct 200 -ct 300 \
;'.format(args.tmp_dir, GenomeAnalysisTK, args.ref ,bam, coverage_output_prefix, args.targetlimit, args.genelist)
			doc_cmds.append(doc_cmd)
			doc_files_to_consolidate.append(coverage_output_prefix+'*')
	
	#submit cluster job array for all bam file metrics
	if not args.NoPicardMetrics:
		jobName = picard_jobName
		job_dir = picard_job_dir
		my.print_log('job {} started; commands:\n{}'.format(jobName, '\n'.join(cmds)))
		j = job.Job(jobName = jobName,
		outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = job_dir,
		email = args.email, blockEmail = not args.sendemailforintermediatejobs)
		try:
			j.executeCommands(cmds, synchronous=True)
		except Exception as e:
			my.print_err('There was an error running job {}\n{}\nintermediate files are in {}\n'.format(jobName, e, job_dir))
		else:
			with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
				s.write(format('\n'.join(j.completionMsg)))

	#submit cluster job array for all bam file depth of coverage
	if not args.NoDepthOfCoverage:
		jobName = doc_jobName
		job_dir = doc_job_dir
		cmds = doc_cmds
		my.print_log('job {} started; commands:\n{}'.format(jobName, '\n'.join(cmds)))
		j = job.Job(jobName = jobName,
		outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = job_dir,
		email = args.email, blockEmail = not args.sendemailforintermediatejobs,
		processors = 8, memory = '24G', queuingModeBunchJobsOnNodes=False)
		try:
			j.executeCommands(cmds, synchronous=True)
		except Exception as e:
			my.print_err('There was an error running job {}\n{}\nintermediate files are in {}\n'.format(jobName, e, job_dir))
		else:
			with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
				s.write(format('\n'.join(j.completionMsg)))
	
if __name__ == "__main__": sys.exit(main())



