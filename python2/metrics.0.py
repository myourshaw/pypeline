#!/usr/bin/env python

import os
import sys
import argparse
import ConfigParser #configparser in python 3
import glob
import tempfile
from time import localtime, strftime
import drmaa
from getpass import getuser

#/data/storage-1-02/archive/myourshaw/gmd/reads/readgroup_bams/*.bam

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
		for path in glob.glob(pathname):
			result.append(os.path.expanduser(os.path.expandvars(path)))
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

def main():
	
	USER = getuser()
	HERE = os.path.dirname(sys.argv[0])
	
	#command line arguments
	parser = argparse.ArgumentParser(
		description = 'Calculate bam file metrics using picard and gatk DepthOfCoverage', # main description for help
		epilog = '\xc2 2011 Michael Yourshaw all rights reserved') # displayed after help
	parser.add_argument('--ini',nargs='?', default=HERE+'/b37.ini',
											help='.ini configuration file (default: '+HERE+'/b37.ini')
	parser.add_argument('--email','-m', metavar='EMAIL', nargs='+', default=USER+'@ucla.edu',
											help='email addresses (default: '+USER+'@ucla.edu)')
	parser.add_argument('--NoEmail', action='store_true', default=False, help='do not send email notifications')
	parser.add_argument('--ref','-R', nargs='?',
											help='reference genome fasta (default: b37.ini->b37_fasta)')
	parser.add_argument('--genelist','-G', nargs='?',
											help='GATK DepthOfCoverage geneList (default: b37.ini->refgene_b37_sorted)')
	parser.add_argument('--targetlimit','-L', nargs='?',
											help=' (default: bait interval list)')
	parser.add_argument('--bait','-B', metavar='BAIT_INTERVAL_LIST', nargs='?',
											help='bait interval list (default: b37.ini->sureselect_50mb_interval_list)')
	parser.add_argument('--targets','-T', metavar='TARGET_INTERVAL_LIST', nargs='+',
											help='list of target interval_list files (default: <bait interval list> b37.ini->ensembl_all_genes b37.ini->ensembl_protein_coding_genes)')
	parser.add_argument('--regex','-X', nargs='?',default=r'[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*',
											help='read name regex (default for illumina/novoalign: "[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*")')
	parser.add_argument('--pixels','-D', nargs='?', type=int, default=100,
											help='optical duplicate pixel distance (default for illumina: 100)')
	parser.add_argument('--output_dir', '-O', nargs='?', default=os.getcwd()+'/metrics',
											help='output directory (default: '+os.getcwd()+'/metrics)')
	parser.add_argument('--consolidate', metavar='METRIC', nargs='+',
											help='additional files to include in consolidation (e.g., *hahahaaha.markdup.metrics')
	parser.add_argument('bams', metavar='BAM', nargs='+',
											help='input list of bam files')
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
	
	#get defaults from configuration ini file
	config = ConfigParser.SafeConfigParser()
	config.read(args.ini)
	if args.ref == None: args.ref = config.get('reference','b37_fasta')
	if args.genelist == None: args.genelist = config.get('gene','refgene_b37_sorted')
	if args.bait == None: args.bait = config.get('interval','sureselect_50mb_interval_list')
	if args.targetlimit == None: args.targetlimit = args.bait
	if args.targets == None: args.targets = [args.bait, config.get('interval','ensembl_all_genes'),config.get('interval','ensembl_protein_coding_genes')]
	
	GenomeAnalysisTK = config.get('gatk','GenomeAnalysisTK')
	CalculateHsMetrics = config.get('picard','CalculateHsMetrics')
	BamIndexStats = config.get('picard','BamIndexStats')
	CollectMultipleMetrics = config.get('picard','CollectMultipleMetrics')
	CollectGcBiasMetrics = config.get('picard','CollectGcBiasMetrics')
	EstimateLibraryComplexity = config.get('picard','EstimateLibraryComplexity')
	MetricsConsolidate = config.get('picard','MetricsConsolidate')
	samtools = config.get('samtools','samtools')
	allq = config.get('qsub','allq')
	himemq = config.get('qsub','himemq')
	tempdir = os.path.join(config.get('DEFAULT','tmpscratch'),USER)
	
	#constants
	CONSOLIDATED_DIR = 'consolidated_metrics'
	COVERAGE_DIR = 'DepthOfCoverage'
	QOUT_DIR = 'qout'
	RUNNER_FILE = os.path.join(os.getcwd(), 'drmaa_job_array_runner.sh')
	
	FILE_MISSING_MSG = 'These required files are missing. You may need to edit {0}.\n\t{1}'
	RUN_JOB_MSG = """Job: {0} [{1}]
Command: {2}
STDOUT and STDERR in {3}.
Notification to {4}."""
	EXCEPTION_MSG = 'ERROR: {0}; {1}'
	
	PICARD_CMD = '$(which java) -Xmx5g -Djava.io.tmpdir='+tempdir+ \
' -jar {0} \
TMP_DIR='+tempdir+ \
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
	bams = list_files(args.bams)
	required_files = [args.ref, args.genelist, args.targetlimit, args.bait, args.targets, bams,
		GenomeAnalysisTK, CalculateHsMetrics, BamIndexStats, CollectMultipleMetrics,
		CollectGcBiasMetrics, EstimateLibraryComplexity, MetricsConsolidate, samtools, RUNNER_FILE]
	missing_files = list_missing_files(required_files)
	if len(bams) == 0 or missing_files != None:
		print FILE_MISSING_MSG.format(args.ini, '\n\t'.join(missing_files))
		sys.exit(1)
	
	#output directories
	if args.output_dir == None: args.output_dir = os.path.dirname(os.path.abspath(args.bams[0]))
	output_dir = args.output_dir
	makedir(output_dir)
	qout_dir = os.path.join(output_dir,QOUT_DIR)
	makedir(qout_dir)
	consolidated_dir = os.path.join(output_dir,CONSOLIDATED_DIR)
	makedir(consolidated_dir)
	coverage_dir = os.path.join(output_dir,COVERAGE_DIR)
	makedir(coverage_dir)
	
	#print run info
	print """Bam Metrics run {0}
Reference: {1}
Depth of coverage gene list: {2}
Depth of coverage limit: {3}
Bait: {4}
Target(s): {5}
Bam file(s): {6}
Metrics will be saved in: {7}
Consolidated metrics will be saved in: {8}
STDOUT and STDERR will be saved in: {9}"""\
.format(strftime("%a, %d %b %Y %H:%M:%S +0000", localtime()),args.ref, args.genelist, args.targetlimit, args.bait, args.targets,bams,output_dir,consolidated_dir,qout_dir)
	print 'Email notifications will {0}.'.format(\
		'not be sent' if args.NoEmail else 'be sent to: {0}'.format(args.email))

	job_array_commands = []
	doc_job_array_commands = []
	files_to_consolidate = []
	doc_files_to_consolidate = []

	#get metrics for each bam file
	#try:
	for bam in bams:
		#output files
		metricsbase = os.path.join(output_dir, os.path.basename(bam)+'.metrics')
		validate = metricsbase + '.validate'
		idxstats = metricsbase + '.idxstats'
		indexstats = metricsbase + '.indexstats'
		gcbiasmetrics = metricsbase + '.gcbias.table'
		gcbiasmetricschart = metricsbase + '.gcbias.pdf'
		gcbiasmetricssummary = metricsbase + '.gcbias.summary'
		librarycomplexity = metricsbase + '.librarycomplexity'
		coverage_output_prefix = os.path.join(coverage_dir, os.path.basename(bam))
		

		if not args.NoPicardMetrics:
			#idxstats
			if not args.Noidxstats:
				job_array_commands.append('{0} idxstats {1} > {2};'.format (samtools, bam, idxstats))
			
			#BamIndexStats
			if not args.NoBamIndexStats:
				job_array_commands.append(PICARD_CMD_BAMINDEXSTATS.format(BamIndexStats, bam, indexstats))
	
			#CalculateHsMetrics
			if not args.NoCalculateHsMetrics:
				for target in args.targets:
					hsoutput = '{0}.target_{1}.hsmetrics'.format(metricsbase,os.path.basename(target))
					job_array_commands.append(PICARD_CMD_HSMETRICS.format(CalculateHsMetrics, bam, hsoutput, args.bait, target))
					files_to_consolidate.append(hsoutput)
		
			#CollectMultipleMetrics
			if not args.NoCollectAlignmentSummaryMetrics or not args.NoCollectInsertSizeMetrics \
			or not args.NoQualityScoreDistribution or not args.NoMeanQualityByCycle:
				job_array_commands.append(PICARD_CMD_MULTIPLEMETRICS.format(CollectMultipleMetrics, bam, metricsbase, args.ref, \
					'' if args.NoCollectAlignmentSummaryMetrics else 'PROGRAM=CollectAlignmentSummaryMetrics', \
					'' if args.NoCollectInsertSizeMetrics else 'PROGRAM=CollectInsertSizeMetrics', \
					'' if args.NoQualityScoreDistribution else 'PROGRAM=QualityScoreDistribution', \
					'' if args.NoMeanQualityByCycle else 'PROGRAM=MeanQualityByCycle'))
				files_to_consolidate.append(metricsbase+'.*_metrics')
	
			#CollectGcBiasMetrics
			if not args.NoCollectGcBiasMetrics:
				job_array_commands.append(PICARD_CMD_GCBIAS.format(CollectGcBiasMetrics, bam, gcbiasmetrics, args.ref, gcbiasmetricschart, gcbiasmetricssummary))
				files_to_consolidate.append(gcbiasmetrics)
				files_to_consolidate.append(gcbiasmetricssummary)
	
			#EstimateLibraryComplexity
			if not args.NoEstimateLibraryComplexity:
				job_array_commands.append(PICARD_CMD_LIBRARYCOMPLEXITY.format(EstimateLibraryComplexity, bam, librarycomplexity, args.regex, args.pixels))
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
;'.format(tempdir,GenomeAnalysisTK,args.ref,bam,coverage_output_prefix,args.targetlimit,args.genelist)
			doc_job_array_commands.append(doc_cmd)
			doc_files_to_consolidate.append(coverage_output_prefix+'*')
	
	#submit cluster job array for all bam file metrics
	#command file for job array
	if not args.NoPicardMetrics:
		command_file, command_filename = unique_file(os.path.join(output_dir,'BamMetrics.commands'))
		command_file.writelines("%s\n" % item for item in job_array_commands)
		command_file.close()
		try:
			met_s = drmaa.Session()
			met_s.initialize()
			met_jt = met_s.createJobTemplate()
			met_jt.WORKING_DIRECTORY = output_dir
			met_jt.remoteCommand = RUNNER_FILE
			met_jt.email = [args.email]
			met_jt.blockEmail = args.NoEmail
			met_jt.outputPath = ':'+qout_dir
			met_jt.errorPath = met_jt.outputPath
			met_jt.nativeSpecification = allq
			met_jt.args = [command_filename]
			met_jt.jobName = 'BamMetrics'
			joblist = met_s.runBulkJobs(met_jt,1,len(job_array_commands),1)
			print 'Bam Metrics job array commands: {0}\nBam Metrics job ids: {1}'.format(command_filename,joblist)
		except Exception as e:
			print EXCEPTION_MSG.format(e, 'bam metrics')
		finally:
			met_s.deleteJobTemplate(met_jt)
			met_s.exit()
			print 'bam metrics session exited'

	#submit cluster job array for all bam file depth of coverage
	#command file for job array
	if not args.NoDepthOfCoverage:
		doc_command_file, doc_command_filename = unique_file(os.path.join(output_dir,'DepthOfCoverage.commands'))
		doc_command_file.writelines("%s\n" % item for item in doc_job_array_commands)
		doc_command_file.close()
		try:
			doc_s = drmaa.Session()
			doc_s.initialize()
			doc_jt = doc_s.createJobTemplate()
			doc_jt.WORKING_DIRECTORY = output_dir
			doc_jt.remoteCommand = RUNNER_FILE
			doc_jt.email = [args.email]
			doc_jt.blockEmail = args.NoEmail
			doc_jt.outputPath = ':'+qout_dir
			doc_jt.errorPath = doc_jt.outputPath
			doc_jt.nativeSpecification = himemq
			doc_jt.args = [doc_command_filename]
			doc_jt.jobName = 'DepthOfCoverage'
			doc_joblist = doc_s.runBulkJobs(doc_jt,1,len(doc_job_array_commands),1)
			print 'Depth of Coverage job array commands: {0}\nDepth of Coverage job ids: {1}'.format(doc_command_filename,doc_joblist)
		except Exception as e:
			print EXCEPTION_MSG.format(e,'depth of coverage')
		finally:
			doc_s.deleteJobTemplate(doc_jt)
			doc_s.exit()
			print 'depth of coverage session exited'


if __name__ == "__main__": sys.exit(main())



