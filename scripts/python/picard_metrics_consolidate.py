#!/usr/bin/env python

# Copyright 2011 Michael Yourshaw, all rights reserved

#-d /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/pre-rmdup_metrics/summary /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/pre-rmdup_metrics/* 
#-d /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/post-rmdup_metrics/summary /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/post-rmdup_metrics/* 
#-d /Users/myourshaw/lab/GMD/2011-04_sequencing/metrics_test /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/post-rmdup_metrics/GMD* /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/rmdup_metrics/GMD*

import sys,os
import re
import traceback
from optparse import OptionParser, OptionGroup
import glob

# Init cmd-line args
description = '''
Consolidates multiple picard metrics files into a nice format for a spreadsheet.
'''
parser = OptionParser(description=description,usage="%prog -d <output dir> -i <metrics file list>", version="%prog 1.0.beta.1")
parser.add_option('-d', '--dir', action="store", type="string", dest='dir', help='Directory for output files')
parser.add_option("-v", action="store_true", dest="verbose")
parser.add_option("-q", action="store_false", dest="verbose", default=False)
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
	
def out(file,string):
	if (not string.endswith('\n')): string += '\n'
	if (options.verbose):
		print(string)
	file.write(string)
	
def parse_metrics(output_dir,input_files,metrics_class,header_startswith):
	if (len(input_files) > 0):
		output_file_name = os.path.join(output_dir,metrics_class+'.txt')
		output_file = open(output_file_name,'w')
		for file in input_files:
			out (output_file,'## %s' % (file))
		header_written = False
		for file in input_files:
			file_id = os.path.join(os.path.basename(os.path.dirname(file)),os.path.basename(file))
			linecount = 0
			metrics_class_found = False
			for line in open(file):
				linecount += 1
				if (line.strip() == ''): continue
				if (line.startswith('## HISTOGRAM')):
					break
				if (not metrics_class_found and line.startswith('## METRICS CLASS') and line.rstrip().endswith('.'+metrics_class)):
					metrics_class_found = True
					continue
				if(metrics_class_found):
					if(not header_written):
						if(line.startswith(header_startswith)):
							out(output_file,'%s\t%s\t%s' % ('METRICS_CLASS','FILE',line))
							header_written = True
							continue
					else:
						fields = line.rstrip('\n').split('\t')
						if (len(fields) > 0 and fields[0] != header_startswith):
							out(output_file,'%s\t%s\t%s' % (metrics_class,file_id,line))

def parse_histograms(output_dir,input_files,metrics_class,header_startswith):
	if (len(input_files) > 0):
		output_file_name = os.path.join(output_dir,metrics_class+'.txt')
		output_file = open(output_file_name,'w')
		parallel_output_file_name = os.path.join(output_dir,metrics_class+'.parallel_histograms.txt')
		parallel_output_file = open(parallel_output_file_name,'w')
		for file in input_files:
			out (output_file,'## %s' % (file))
			out (parallel_output_file,'## %s' % (file))
		header_written = False
		out_list=[]
		header = ()
		min_x = None
		max_x = None
		x_set = set()
		for file in input_files:
			file_id = os.path.join(os.path.basename(os.path.dirname(file)),os.path.basename(file))
			linecount = 0
			histogram_found = False
			for line in open(file):
				linecount += 1
				if (line.strip() == ''): continue
				if (not histogram_found and line.startswith('## HISTOGRAM')):
					histogram_found = True
					continue
				if(histogram_found):
					if(not header_written):
						if(line.startswith(header_startswith)):
							out(output_file,'%s\t%s\t%s' % ('METRICS_CLASS','FILE',line))
							header_fields = line.rstrip('\n').split('\t')
							header = (header_fields[0],header_fields[1])
							header_written = True
							out(parallel_output_file,"## %s for each file" % header_fields[1])
							out(parallel_output_file,"%s\t%s" % (header_fields[0],'\t'.join(input_files)))
							continue
					else:
						fields = line.rstrip('\n').split('\t')
						if (len(fields) > 0 and fields[0] != header_startswith):
							out(output_file,'%s\t%s\t%s' % (metrics_class,file_id,line))
							x = float(fields[0])
							y = float(fields[1])
							out_list.append((file,x,y))
							x_set.add(x)
							if (min_x == None or x < min_x): min_x = x
							if (max_x == None or x > max_x): max_x = x
		x_list = list(x_set)
		x_list.sort()
		for i in x_list:
			line_list = ['0']*len(input_files)
			these_x = [o for o in out_list if o[1] == i]
			for fxy in these_x:
				pos = input_files.index(fxy[0])
				line_list[pos] = str(fxy[2])
			out(parallel_output_file,"%s\t%s" % (i,'\t'.join(line_list)) )

def main(argv=None):
	if argv is None: argv = sys.argv
	output_dir = options.dir
	makedir(output_dir)
	
	files_alignment_summary_metrics = []
	files_gcbias_summary = []
	files_gcbias_table = []
	files_hsmetrics = []
	files_idxstats = []
	files_indexstats = []
	files_insert_size_metrics = []
	files_insert_size_histograms = []
	files_librarycomplexity = []
	files_librarycomplexityhistograms = []
	files_quality_by_cycle_metrics = []
	files_quality_distribution_metrics = []
	files_duplication_metrics = []
	files_duplication_histograms = []
	for arg in args:
		for file in get_files(glob.iglob(arg)):
			if file.endswith('.pdf'): continue
			elif file.endswith('alignment_summary_metrics'): files_alignment_summary_metrics.append(file)
			elif file.endswith('.gcbias.summary'): files_gcbias_summary.append(file)
			elif file.endswith('.gcbias.table'): files_gcbias_table.append(file)
			elif file.endswith('.hsmetrics'): files_hsmetrics.append(file)
			elif file.endswith('.idxstats'): files_idxstats.append(file)
			elif file.endswith('.indexstats'): files_indexstats.append(file)
			elif file.endswith('.insert_size_metrics'):
				files_insert_size_metrics.append(file)
				files_insert_size_histograms.append(file)
			elif file.endswith('.librarycomplexity'):
				files_librarycomplexity.append(file)
				files_librarycomplexityhistograms.append(file)
			elif file.endswith('.quality_by_cycle_metrics'): files_quality_by_cycle_metrics.append(file)
			elif file.endswith('.quality_distribution_metrics'): files_quality_distribution_metrics.append(file)
			elif file.endswith('.metrics'):
				files_duplication_metrics.append(file)
				files_duplication_histograms.append(file)
	
	parse_metrics(output_dir,files_alignment_summary_metrics,'AlignmentSummaryMetrics','CATEGORY')
	parse_metrics(output_dir,files_gcbias_summary,'GcBiasSummaryMetrics','WINDOW_SIZE')
	parse_metrics(output_dir,files_gcbias_table,'GcBiasDetailMetrics','GC')
	parse_metrics(output_dir,files_hsmetrics,'HsMetrics','BAIT_SET')
	#parse_metrics_class_files(output_dir,files_idxstats,'','')
	#parse_metrics_class_files(output_dir,files_indexstats,'','')
	parse_metrics(output_dir,files_insert_size_metrics,'InsertSizeMetrics','MEDIAN_INSERT_SIZE')
	parse_histograms(output_dir,files_insert_size_histograms,'InsertSizeHistograms','insert_size')
	parse_metrics(output_dir,files_librarycomplexity,'DuplicationMetrics','LIBRARY')
	parse_histograms(output_dir,files_librarycomplexityhistograms,'DuplicationHistograms','duplication_group_count')
	parse_histograms(output_dir,files_quality_by_cycle_metrics,'QualityByCycleHistograms','CYCLE')
	parse_histograms(output_dir,files_quality_distribution_metrics,'QualityDistributionsHistograms','QUALITY')
	parse_metrics(output_dir,files_duplication_metrics,'DuplicationMetrics','LIBRARY')
	parse_histograms(output_dir,files_duplication_histograms,'DuplicationHistograms','BIN')
	
if __name__ == "__main__": sys.exit(main())
