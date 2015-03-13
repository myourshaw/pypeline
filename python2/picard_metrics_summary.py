#!/usr/bin/env python

import os
import re
import traceback
from optparse import OptionParser, OptionGroup
from glob import iglob

# Init cmd-line args
description = """
foo bar
"""

parser = OptionParser(description=description)
parser.add_option("-d", "--dir", metavar="DIR", dest="dir", help="Directory for output files")
(options, args) = parser.parse_args()

def makedir(dir):
	try:
	  os.makedirs(dir)
	except OSError:
		if os.path.isdir(dir):
			pass
		else:
			raise

def get_files(names):
	return (file for file in names if os.path.isfile(file))

def is_number(s):
	try:
		n = float(s)
		if (math.isinf(n) or math.isnan(n)): return False
		return True
	except ValueError:
		return False

outputDir = options.dir
makedir(outputDir)

#-d /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/pre-rmdup_metrics/summary /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/pre-rmdup_metrics/* 
#-d /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/post-rmdup_metrics/summary /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/post-rmdup_metrics/* 

AlignmentSummaryMetrics_file = open(os.path.join(outputDir,'AlignmentSummaryMetrics.txt'),"w")
GcBiasSummaryMetrics_file = open(os.path.join(outputDir,'GcBiasSummaryMetrics.txt'),"w")
HsMetrics_file = open(os.path.join(outputDir,'HybridSelectionMetrics.txt'),"w")
InsertSizeMetrics_file = open(os.path.join(outputDir,'InsertSizeMetrics.txt'),"w")
DuplicationMetrics_file = open(os.path.join(outputDir,'DuplicationMetrics.txt'),"w")
AlignmentSummaryMetrics_header_written = False
GcBiasSummaryMetrics_header_written = False
HsMetrics_header_written = False
InsertSizeMetrics_header_written = False
DuplicationMetrics_header_written = False
id_header='metrics_class\texperiment\tdate\tmachine\trun\tlane\tbarcode'
for arg in args:
	for file in get_files(iglob(arg)):
		if file.endswith(".pdf"): continue
		file_base = os.path.basename(file)
		file_fields = file_base.split('.')
		id_string = '\t'.join(file_fields[0:6])
		linecount = 0
		metrics_class = ''
		for line in open(file):
			linecount += 1
			if (line.strip() == ''): continue
			if (metrics_class == ''):
				if (line.startswith("# net.sf.picard.analysis.")):
					analysis = re.sub(r'^# net.sf.picard.analysis\.','',line).strip()
					continue
				elif (line.startswith("## METRICS CLASS")):
					if (line.strip().endswith('.AlignmentSummaryMetrics')):
						metrics_class = 'AlignmentSummaryMetrics'
						#CATEGORY	TOTAL_READS	PF_READS	PCT_PF_READS	PF_NOISE_READS	PF_READS_ALIGNED	PCT_PF_READS_ALIGNED	PF_ALIGNED_BASES	PF_HQ_ALIGNED_READS	PF_HQ_ALIGNED_BASES	PF_HQ_ALIGNED_Q20_BASES	PF_HQ_MEDIAN_MISMATCHES	PF_MISMATCH_RATE	PF_HQ_ERROR_RATE	MEAN_READ_LENGTH	READS_ALIGNED_IN_PAIRS	PCT_READS_ALIGNED_IN_PAIRS	BAD_CYCLES	STRAND_BALANCE	PCT_CHIMERAS	PCT_ADAPTER
					elif (line.strip().endswith('.GcBiasSummaryMetrics')):
						metrics_class = 'GcBiasSummaryMetrics'
						#WINDOW_SIZE	TOTAL_CLUSTERS	ALIGNED_READS	TOTAL_BIAS	LOW_GC_BIAS	MID_GC_BIAS	HIGH_GC_BIAS	JAFFE_BIAS_METRIC
					elif (line.strip().endswith('.HsMetrics')):
						metrics_class = 'HsMetrics'
						#BAIT_SET	GENOME_SIZE	BAIT_TERRITORY	TARGET_TERRITORY	BAIT_DESIGN_EFFICIENCY	TOTAL_READS	PF_READS	PF_UNIQUE_READS	PCT_PF_READS	PCT_PF_UQ_READS	PF_UQ_READS_ALIGNED	PCT_PF_UQ_READS_ALIGNED	PF_UQ_BASES_ALIGNED	ON_BAIT_BASES	NEAR_BAIT_BASES	OFF_BAIT_BASES	ON_TARGET_BASES	PCT_SELECTED_BASES	PCT_OFF_BAIT	ON_BAIT_VS_SELECTED	MEAN_BAIT_COVERAGE	MEAN_TARGET_COVERAGE	PCT_USABLE_BASES_ON_BAIT	PCT_USABLE_BASES_ON_TARGET	FOLD_ENRICHMENT	ZERO_CVG_TARGETS_PCT	FOLD_80_BASE_PENALTY	PCT_TARGET_BASES_2X	PCT_TARGET_BASES_10X	PCT_TARGET_BASES_20X	PCT_TARGET_BASES_30X	HS_LIBRARY_SIZE	HS_PENALTY_10X	HS_PENALTY_20X	HS_PENALTY_30X
					elif (line.strip().endswith('.InsertSizeMetrics')):
						metrics_class = 'InsertSizeMetrics'
						InsertSizeMetrics_histogram = False
						#MEDIAN_INSERT_SIZE	MIN_INSERT_SIZE	MAX_INSERT_SIZE	MEAN_INSERT_SIZE	STANDARD_DEVIATION	READ_PAIRS	PAIR_ORIENTATION	WIDTH_OF_10_PERCENT	WIDTH_OF_20_PERCENT	WIDTH_OF_30_PERCENT	WIDTH_OF_40_PERCENT	WIDTH_OF_50_PERCENT	WIDTH_OF_60_PERCENT	WIDTH_OF_70_PERCENT	WIDTH_OF_80_PERCENT	WIDTH_OF_90_PERCENT	WIDTH_OF_99_PERCENT
					elif (line.strip().endswith('.DuplicationMetrics')):
						metrics_class = 'DuplicationMetrics'
						DuplicationMetrics_histogram = False
						#LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
				continue
			else:
				#metrics_class is known
				if(metrics_class == 'AlignmentSummaryMetrics'):
					if(not AlignmentSummaryMetrics_header_written):
						if(line.startswith('CATEGORY')):
							AlignmentSummaryMetrics_file.write("%s\t%s" % (id_header,line))
							AlignmentSummaryMetrics_header_written = True
					else:
						fields = line.rstrip("\n").split("\t")
						if (len(fields) > 0 and fields[0] in ('FIRST_OF_PAIR','SECOND_OF_PAIR','PAIR')):
							AlignmentSummaryMetrics_file.write("Alignment\t%s\t%s" % (id_string,line))
							
				elif(metrics_class == 'GcBiasSummaryMetrics'):
					if(not GcBiasSummaryMetrics_header_written):
						if(line.startswith('WINDOW_SIZE')):
							GcBiasSummaryMetrics_file.write("%s\t%s" % (id_header,line))
							GcBiasSummaryMetrics_header_written = True
					else:
						fields = line.rstrip("\n").split("\t")
						if (len(fields) > 0 and fields[0] != 'WINDOW_SIZE'):
							GcBiasSummaryMetrics_file.write("GcBias\t%s\t%s" % (id_string,line))
							
				elif(metrics_class == 'HsMetrics'):
					if(not HsMetrics_header_written):
						if(line.startswith('BAIT_SET')):
							HsMetrics_file.write("%s\t%s" % (id_header,line))
							HsMetrics_header_written = True
					else:
						fields = line.rstrip("\n").split("\t")
						if (len(fields) > 0 and fields[0] != 'BAIT_SET'):
							HsMetrics_file.write("HybridSelection\t%s\t%s" % (id_string,line))
							
				elif(metrics_class == 'InsertSizeMetrics'):
					if(line.startswith('## HISTOGRAM')): InsertSizeMetrics_histogram = True
					if (not InsertSizeMetrics_histogram):
						if(not InsertSizeMetrics_header_written):
							if(line.startswith('MEDIAN_INSERT_SIZE')):
								InsertSizeMetrics_file.write("%s\t%s" % (id_header,line))
								InsertSizeMetrics_header_written = True
						else:
							fields = line.rstrip("\n").split("\t")
							if (len(fields) > 0 and fields[0] != 'MEDIAN_INSERT_SIZE'):
								InsertSizeMetrics_file.write("InsertSize\t%s\t%s" % (id_string,line))
					
				elif(metrics_class == 'DuplicationMetrics'):
					if(line.startswith('## HISTOGRAM')): DuplicationMetrics_histogram = True
					if (not DuplicationMetrics_histogram):
						if(not DuplicationMetrics_header_written):
							if(line.startswith('LIBRARY')):
								DuplicationMetrics_file.write("%s\t%s" % (id_header,line))
								DuplicationMetrics_header_written = True
						else:
							fields = line.rstrip("\n").split("\t")
							if (len(fields) > 0 and fields[0] != 'LIBRARY'):
								DuplicationMetrics_file.write("Duplication\t%s\t%s" % (id_string,line))
				
