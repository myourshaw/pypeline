#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import shutil
import glob
import my
import job

class Vcf2TpaigeError(Exception): pass

#--keeptempfiles -i /scratch0/tmp/myourshaw/gmd/analysis/vcf/vep_test/GMD151.GMD151_XT_16.gatk.snpFiltered.indelFiltered.recalibrated.pos.vcf

def main():

	#command line arguments
	parser = argparse.ArgumentParser(
		description = 'run Variant Effect Predictor for multiple vcfs',
		epilog = 'pypeline.vcf2vep version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--input', '-i', nargs='+',
											help='input vcf file(s); output will be <vcf>.tpaige.vep')
	parser.add_argument('--vep',
											help='path to tpaige perl executable')
	parser.add_argument('--records_per_job', default=1000,
											help='number of vcf records to process in each job (default: 1000)')
	parser.add_argument('--keeptempfiles', action='store_true', default=False,
											help='keep temporary files')
	parser.add_argument('--email', nargs='*',
											help='email address for job completion notifications')
	parser.add_argument('--config',
											help='additional (overriding?) configuration file (default=[pypeline.cfg, ~/pypeline.cfg])')
	args = parser.parse_args()
	config = my.get_config(args)
	if not args.vep:
		args.vep = config.get('ensembl','tpaige')
	vep_cache = config.get('ensembl','vep_cache')
	if not os.path.isdir(vep_cache):
		vep_cache = os.path.join(os.path.expanduser('~'),'.vep')

	vcfs = sorted(my.flatten([glob.glob(v) for v in args.input]))
	all_lines_count = 0
	for vcf in vcfs:
		output_dir = os.path.dirname(vcf)
		if not output_dir:
			output_dir = os.getcwd()
		my.makedir(output_dir)
		vcf_base = os.path.basename(vcf)
		vep_out = os.path.join(output_dir, vcf_base+'.tpaige.vep')
		tempdir = my.unique_dir(os.path.join(output_dir,'tmp'))
		qout = os.path.join(tempdir, 'qout')
		my.makedir(qout)
		part_base = '{}.part.'.format(os.path.join(tempdir, vcf_base))
		cmd = 'split -a 8 -d -l {} {} {}'.format(args.records_per_job, vcf, part_base)
		j = my.run_job(cmd, jobName='split_'+vcf_base, job_dir=tempdir, qout_dir=qout, email=None, synchronous=True, processors=None, memory=None, hold_jid=None)
		#split_id = j.jobId
		cmds = []
		for part in glob.iglob('{}*'.format(part_base)):
			cmds.append('perl {} \
--no_progress \
--species human \
--input_file {} \
--format vcf \
--output_file {} \
--use_ncbi \
--check_existing=1 \
--failed=0 \
--use_uniprot \
--force_overwrite \
--host cortex.local \
--user ensembl \
--password ensembl \
--port 3306 \
--terms so \
--sift=b \
--polyphen=b \
--condel=b \
--regulatory \
--hgvs \
--gene \
--protein \
--hgnc \
--buffer_size 5000 \
'.format(args.vep, part, part+'.vep'))
		j = my.run_job(cmds, jobName='vep_'+vcf_base, job_dir=tempdir, qout_dir=qout, email=None, synchronous=True, processors=None, memory=None, hold_jid=None)
		#vep_id = j.jobId
		vep_parts = sorted(glob.glob('{}*.vep'.format(part_base)))
		header_fd, header_file = my.unique_file(os.path.join(tempdir,'header'))
		with open(vep_parts[0]) as header_vep:
			for line in header_vep:
				if line.startswith('#'):
					header_fd.write(line)
				else:
					break
		header_fd.close()
		cmd = """find {}*.vep -type f -print0 | xargs -0 cat | grep -v '#' | cat {} - > {}""".format(part_base, header_file, vep_out)
		#cmd = """cat {} | grep -v '#' | cat {} - > {}""".format(' '.join(vep_parts), header_file, vep_out)
		j = my.run_job(cmd, jobName='join_'+vcf_base+'.vep', job_dir=tempdir, qout_dir=qout, email=args.email, synchronous=True, processors=None, memory=None, hold_jid=None)
		#join_id = j[0]

		if not args.keeptempfiles:
			cmd = 'rm -rf {}'.format(tempdir)
			shutil.rmtree(tempdir, ignore_errors=True)

	
if __name__ == "__main__": sys.exit(main())
