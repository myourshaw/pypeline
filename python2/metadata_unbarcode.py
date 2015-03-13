#!/usr/bin/env python

import sys
import os
import argparse
import collections
import glob
import re
import warnings
import pickle
import ConfigParser #configparser in python 3
import glob
import tempfile
from time import localtime, strftime
import my
import metadata
import job

#-n /share/apps/myourshaw/novocraft-current/novobarcode -m /home/myourshaw/lab/git-myourshaw/python/metadata.txt -d /scratch1/tmp/myourshaw/gmd/unbarcode-test -e myourshaw@ucla.edu

class QseqError(Exception): pass
class QseqIOError(QseqError): pass
class QseqFileNotFoundError(QseqIOError): pass
class BarcodeError(Exception): pass
class SequencingMetadataError(Exception): pass


def main():
		#command line arguments
	parser = argparse.ArgumentParser(
		description = 'Use novobarcode to segregate qseq files by barcode and then merge results to one file per readgroup')
	parser.add_argument('--novobarcode','-n', nargs='?',
											help='novobarcode executable (default: b37.ini->novobarcode)')
	parser.add_argument('--metadata', '-m',
											help='metadata file')
	parser.add_argument('--email', '-e', nargs='?',
											help='email for job completion notification')
	parser.add_argument('--output_dir', '-d', nargs='?', default=os.getcwd()+'/unbarcode',
											help='output directory (default: '+os.getcwd()+'/unbarcode)')
	args = parser.parse_args()

	owners, studies, samples, libraries, readgroups, barcodefilegroup = metadata.table2metadata(args.metadata)
	
	tmp_dir = my.unique_dir('unbarcode', args.output_dir)
	qout_dir = os.path.join(tmp_dir, 'qout')
	my.makedir(qout_dir)
	
	cmds = []
	output_dir_list = ['owner\texperiment_run_id\tlane\tnovobarcode_output_dir']
	output_dir_msg = ['owner\texperiment_run_id\tlane\tnovobarcode_output_dir']
	
	for bfg in sorted(barcodefilegroup):

		owner = bfg.owner
		exprun = bfg.experiment_run_id
		lane = bfg.lane
		
		#create barcode file with barcodes common to the set of read files
		barcode_set = barcodefilegroup[bfg].barcode_set
		barcode_file, barcode_filename = my.unique_file(os.path.join(tmp_dir,'unbarcode.barcodes'))
		barcode_file.writelines('Distance 2\nFormat N 5\n{}\n'.format(str(barcode_set)))
		barcode_file.close()
		qseq_trios = barcodefilegroup[bfg].qseq_trios
		
		#create job array
		for q in sorted(qseq_trios):
			novobarcode_output_dir = os.path.join(args.output_dir, bfg.owner, 'runs', bfg.experiment_run_id, bfg.lane)
			my.makedir(novobarcode_output_dir)
			output_dir_list.append(novobarcode_output_dir)
			output_dir_msg.append('{}\t{}\t{}\t{}'.format(bfg.owner, bfg.experiment_run_id, bfg.lane, novobarcode_output_dir))
			cmd='{0} -b {1} -d {2} -F QSEQ -f {3} {4} -i {5} --ILQ_SKIP --QSEQ_OUT --NC_OFF;'\
			 .format(args.novobarcode, barcode_filename, novobarcode_output_dir, q.read1, q.read3, q.read2)
			cmds.append(cmd)

	j = job.Job(jobName = 'novobarcode_{}'.format(args.output_dir),
		outputPath = qout_dir, errorPath = qout_dir, workingDirectory = tmp_dir,
		blockEmail = False if args.email else True, email = args.email)
	
	my.print_log('starting novobarcode\nthis run may take a while ...\nunmerged qseq files will be in:\n{}\nSTDOUT and STDERR will be in [{}]'\
	.format(time,'\n'.join(output_dir_msg),qout_dir))

	try:
		j.executeCommands(cmds, synchronous=True)
	except Exception as __e:
		my.print_err('There was an error running the novobarcode jobs on the cluster:\n{}\n'.format(__e))
	else:
		my.print_log('novobarcode completed\njob status report:\n{}\n'\
		.format('\n'.join(j.completionMsg)))


if __name__ == "__main__": sys.exit(main())
