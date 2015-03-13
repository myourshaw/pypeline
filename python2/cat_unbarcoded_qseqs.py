#!/usr/bin/env python

import sys
import os
import argparse
from glob import glob
import shutil
import cat
import job
from myourshaw import makedir, unique_dir

class MergeQseqError(Exception): pass

def cat_unbarcoded_qseqs(dirs, debug=False):
	cmds = []
	for dir in dirs:
		items = dir.split(os.path.sep)
		output_dir = os.path.sep.join(items[:len(items)-3])
		makedir(output_dir)
		tmp_dir = unique_dir('merge_qseqs_tmp', output_dir)
		qout_dir = os.path.join(tmp_dir,'qout')
		makedir(qout_dir)
		exptrun, lane, barcode = items[-3:]
		prefix = os.path.join(output_dir, exptrun, '{}.{}.{}.'.format(exptrun, lane, barcode))
		files1 = sorted(glob(os.path.join(dir,'s_{}_1_*_qseq.txt'.format(lane))))
		cmd = 'python {} -i {} -o {}'.format(cat.__file__, ' '.join(files1), '{}pe1_qseq.txt'.format(prefix))
		cmds.append(cmd)
		files2 = sorted(glob(os.path.join(dir,'s_{}_3_*_qseq.txt'.format(lane))))
		if files2:
			if len(files1) != len(files2):
				raise MergeQseqError('{} files for paired end 1 and {} files for paired end 2 in [{}]'.format(len(files1), len(files2), dir))
			for i in range(0,len(files1)-1):
				if os.path.basename(files1[i])[:4] != os.path.basename(files2[i])[:4] or os.path.basename(files1[i])[-14:] != os.path.basename(files2[i])[-14:]:
					raise MergeQseqError('paired end 1 file [{}] does not match paired end 2 file [{}]'.format(files1[i], files2[i]))
			cmd = 'python {} -i {} -o {}'.format(cat.__file__, ' '.join(files2), '{}pe2_qseq.txt'.format(prefix))
			cmds.append(cmd)
	j = job.Job(jobName = 'mergeqseqs',
	outputPath = qout_dir, errorPath = qout_dir, workingDirectory = tmp_dir,
	blockEmail = True)
	try:
		j.executeCommands(cmds, synchronous=True)
	except Exception as e:
		print_err('There was an error running the merge qseqs jobs on the cluster:\n{}\nintermediate files are in {}'.format(e, tmp_dir))
	else:
		if debug:
			with open(os.path.join(tmp_dir,'job_status_report.txt')) as stat:
				stat.write(format('\n'.join(j.completionMsg)))
		else:
			shutil.rmtree(tmp_dir)

#-i /scratch1/tmp/myourshaw/gmd/unbarcode-test/myourshaw/runs/110511_SN430_0243_B817FLABXX/1/TTAGGCA

def main():
		#command line arguments
	parser = argparse.ArgumentParser(
		description = 'cat /*/exptrun/lane/barcode/s_[1-8]_[13]_*_qseq.txt > /*/exptrun/exptrun.lane.barcode.pe[12]_qseq.txt',
		epilog = '\a9 2011 Michael Yourshaw all rights reserved')
	parser.add_argument('--input', '-i', nargs='+',
											help='list of barcode directories to merge')
	parser.add_argument('--debug', action='store_true', default=False,
											help='do not delete temp files')
	args = parser.parse_args()
	
	print 'begin merge qseq files. please wait ...'
	cat_unbarcoded_qseqs(args.input)
	print 'merge qseq files finished'
	
if __name__ == "__main__": sys.exit(main())