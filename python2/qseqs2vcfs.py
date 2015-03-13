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
import collections
import my
import job

class Bams2VcfsError(Exception): pass

#-d /scratch0/tmp/myourshaw/tmp-gmd -m /home/myourshaw/lab/git-myourshaw/python/metadata.txt -q /data/storage-1-02/archive/myourshaw/illumina/110623_SN860_0067_2011-100R_A81MVKABXX/Data/Intensities/BaseCalls/QSEQ/s_[1-8]_[123]_*_qseq.txt /data/storage-1-02/solexa/110511_SN430_0243_B817FLABXX/Data/Intensities/BaseCalls/s_[1-7]_[123]_*_qseq.txt /data/storage-1-04/archive/illumina/110616_SN860_0065_2011-101_B817EAABXX/Data/Intensities/BaseCalls/QSEQ/s_[1-7]_[123]_*_qseq.txt


def main():

	#command line arguments
	parser = argparse.ArgumentParser(parents=[my.get_directory_parser()],
		description = 'Starting with markdup bam files (that may be specified in metadata), realign indels, recalibrate base quality scores, genotype, and filter/recalibrate variants; optionally merge cohort(s)',
		epilog = 'pypeline.bams2vcfs version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--qseqs', '-q', nargs='+',
											help='raw reads in qseq[.gz] format',)
	parser.add_argument('--metadata', '-m', required=True,
											help='tab-delimited metadata file',)
	parser.add_argument('--barcode_file',
											help='path to novobarcode-style tag file (default: TruSeq/XT 12 barcodes)',)
	args = parser.parse_args()
	#my.set_directory_defaults(args)
	config = my.get_config(args)

	#get qseqs
	qseqs = set()
	for qs in args.qseqs:
		[qseqs.add(q) for q in glob.glob(qs)]
	qseqs = sorted(list(qseqs))
	qseq_info = []
	[qseq_info.append(my.qseq_peek(q)) for q in qseqs]
	filegroups = set()
	readgroups = set()
	for q in qseq_info:
		q['metadata'] = [m for m in my.get_metadata(args.metadata) if m['machine'] == q['machine'] and m['run'] == q['run'] and m['lane'] == str(q['lane'])]
		if q.get('metadata'):
			for m in q['metadata']:
				filegroups.add(my.FilegroupInfo(q['machine'], q['run'], q['lane']))
				readgroups.add(my.ReadgroupInfo(q['machine'], q['run'], q['lane'], m['barcode'].strip() if m['barcode'].strip() else None))
				if not m['sample'] or not m['library'] or not m['platform'] or not m['predicted_median_insert_size'] or not m['adapter_1'] or not m['adapter_2'] or not m['sequencing_center'] or not m['run_date'] or not m['run_folder']:
					missing_columns = [c for c in ['sample' if not m['sample'] else '', 'library' if not m['library'] else '', 'platform' if not m['platform'] else '', 'predicted_median_insert_size' if not m['predicted_median_insert_size'] else '', 'adapter_1' if not m['adapter_1'] else '', 'adapter_2' if not m['adapter_2'] else '', 'sequencing_center' if not m['sequencing_center'] else '', 'run_date' if not m['run_date'] else '', 'run_folder' if not m['run_folder'] else ''] if c]
					raise Bams2VcfsError('missing metadata information for qseq file {}  with machine={} run={} lane={}: {}'
						.format(q['path'], q['machine'], q['run'], q['lane'], ';'.join(missing_columns)))
		else:
			raise Bams2VcfsError('no metadata for qseq file {} with machine={} run={} lane={}'.format(q['path'], q['machine'], q['run'], q['lane']))
	filegroups = sorted(list(filegroups))
	readgroups = sorted(list(readgroups))
	for fg in filegroups:
		filegroup_qseq_info = [q for q in qseq_info if fg.machine == q['machine'] and fg.run == q['run'] and fg.lane == q['lane']]
		filegroup_read1_files = sorted([r['path'] for r in filegroup_qseq_info if r['read'] == 1])
		filegroup_read2_files = sorted([r['path'] for r in filegroup_qseq_info if r['read'] == 2])
		filegroup_read3_files = sorted([r['path'] for r in filegroup_qseq_info if r['read'] == 3])
		if (filegroup_read2_files and len(filegroup_read2_files) != len(filegroup_read1_files)) or (filegroup_read3_files and len(filegroup_read3_files) != len(filegroup_read1_files)):
			raise Bams2VcfsError('unequal number of qseq files for reads; 1 has {} files, 2 has {} files, 3 has {} files'.format(len(filegroup_read1_files), len(filegroup_read2_files), len(filegroup_read3_files)))
		barcodes = set()
		[barcodes.add(r.barcode) for r in readgroups if r.barcode]
		if barcodes:
			barcodes = sorted(list(barcodes))
			cmds.append('python {} --file1s {} --file2s {} --qseqtagfiles {} --barcodes {} --barcode_file {} --top_dir {} --email {}'
				.format(qseq2rg_novobarcode_and_cat.__file__, ' '.join(filegroup_read1_files), ' '.join(filegroup_read3_files), ' '.join(filegroup_read2_files), ' '.join(barcodes), args.barcode_file, args.top_dir, ' ' .join(args.email)))
		else:
			compressed = 0
			for fq in filegroup_qseq_info:
				compressed += int(fq['compressed'])
				if compressed > 0 and compressed < len(filegroup_qseq_info):
					raise Bams2VcfsError('impermissible mix of compressed and uncompressed qseq files in a filegroup:\n\t{}'.format('\n\t'.join(sorted([r['path'] for r in filegroup_qseq_info]))))
				cmd_str = 'cat {} | gzip > {}' if compressed == 0 else 'gzip -cd {} | cat - | gzip > {}'
				if filegroup_read1_files:
					cmds.append(cmd_str.format(' '.join(filegroup_read1_files), os.path.join(args.readgroup_qseqs_dir,'{}.{}.{}.{}_qseq.txt.gz'.format(fg.machine, fg.run, fg.lane, 1))))
				if filegroup_read2_files:
					cmds.append(cmd_str.format(' '.join(filegroup_read2_files), os.path.join(args.readgroup_qseqs_dir,'{}.{}.{}.{}_qseq.txt.gz'.format(fg.machine, fg.run, fg.lane, 2))))
				if filegroup_read3_files:
					cmds.append(cmd_str.format(' '.join(filegroup_read3_files), os.path.join(args.readgroup_qseqs_dir,'{}.{}.{}.{}_qseq.txt.gz'.format(fg.machine, fg.run, fg.lane, 3))))
	
	jobName = 'qseq2rg_{}'.format(my.localtime_squish())
	my.print_log('job {} started; commands:\n{}'.format(jobName, '\n'.join(cmds)))
	j = job.Job(jobName = jobName,
	outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = readgroup_tmpdir,
	email = args.email, blockEmail = False if args.email else True)

	try:
		j.executeCommands(cmds, synchronous=True)
	except Exception as e:
		my.print_err('There was an error running the job {}:\n{}\nintermediate files are in {}'.format(jobName, e, readgroup_tmpdir))
	else:
		my.print_log('job {} finished; status:\n{}'.format(jobName, '\n'.join(j.completionMsg)))
		with open(os.path.join(readgroup_tmpdir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
			s.write(format('\n'.join(j.completionMsg)))


if __name__ == "__main__": sys.exit(main())
