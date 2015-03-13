#!/usr/bin/env python

import sys
import os
import argparse
from ConfigParser import SafeConfigParser
import shutil
import subprocess
from my import list_missing_files, makedir
from getpass import getuser

class CatFilesMissingInputError(Exception): pass
class CatFilePicardOptionsError(Exception): pass

#-i /scratch1/tmp/myourshaw/gmd/unbarcode-test/myourshaw/runs/110511_SN430_0243_B817FLABXX/1/TTAGGCA/s_1_3_0061_qseq.txt /scratch1/tmp/myourshaw/gmd/unbarcode-test/myourshaw/runs/110511_SN430_0243_B817FLABXX/1/TTAGGCA/s_1_3_0062_qseq.txt /scratch1/tmp/myourshaw/gmd/unbarcode-test/myourshaw/runs/110511_SN430_0243_B817FLABXX/1/TTAGGCA/s_1_3_0063_qseq.txt /scratch1/tmp/myourshaw/gmd/unbarcode-test/myourshaw/runs/110511_SN430_0243_B817FLABXX/1/TTAGGCA/s_1_3_0064_qseq.txt /scratch1/tmp/myourshaw/gmd/unbarcode-test/myourshaw/runs/110511_SN430_0243_B817FLABXX/1/TTAGGCA/s_1_3_0065_qseq.txt /scratch1/tmp/myourshaw/gmd/unbarcode-test/myourshaw/runs/110511_SN430_0243_B817FLABXX/1/TTAGGCA/s_1_3_0066_qseq.txt /scratch1/tmp/myourshaw/gmd/unbarcode-test/myourshaw/runs/110511_SN430_0243_B817FLABXX/1/TTAGGCA/s_1_3_0067_qseq.txt /scratch1/tmp/myourshaw/gmd/unbarcode-test/myourshaw/runs/110511_SN430_0243_B817FLABXX/1/TTAGGCA/s_1_3_0068_qseq.txt -o /scratch1/tmp/myourshaw/gmd/unbarcode-test/cat_test_qseq.txt
#--MergeSamFiles -i /scratch1/tmp/myourshaw/gmd/bams-2011-05-11.HWI-ST430.243/rg/GMD154A_XT_Illumina.2011-05-11.HWI-ST430.243.1.TTAGGCA/GMD154A_XT_Illumina.2011-05-11.HWI-ST430.243.1.TTAGGCA.novoalign.markdup.bam /scratch1/tmp/myourshaw/gmd/bams-2011-05-11.HWI-ST430.243/rg/GMD154A_XT_Illumina.2011-05-11.HWI-ST430.243.2.TTAGGCA/GMD154A_XT_Illumina.2011-05-11.HWI-ST430.243.2.TTAGGCA.novoalign.markdup.bam -o /scratch1/tmp/myourshaw/gmd/unbarcode-test/cat_test.bam

def main():
	
	#command line arguments
	parser = argparse.ArgumentParser(
		description = 'cat input > output',
		epilog = '\a9 2011 Michael Yourshaw all rights reserved')
	parser.add_argument('--input', '-i', nargs='+',
											help='list of files to merge')
	parser.add_argument('--output', '-o',
											help='output file')
	parser.add_argument('--config',
											help='configuration file (default=[pypeline.cfg,~/pypeline.cfg,pypeline.cfg])')
	parser.add_argument('--MergeSamFiles', action='store_true', default=False,
											help='use picard MergeSamFiles (requires )')
	parser.add_argument('--picard_jar',
											help='picard MergeSamFiles jar file (default from config->[picard]/MergeSamFiles)')
	parser.add_argument('--PicardOptions', nargs='?',
											help='picard options in the form OPTION=value [...] (default = )')
	parser.add_argument('--tmpdir',
											help='temporary directory (default from config->tmpscratch)')
	args = parser.parse_args()

	#configuration
	config = SafeConfigParser()
	config.readfp(open(os.path.join(os.path.dirname(__file__),'pypeline.cfg')))
	cfg_files = [os.path.expanduser('~/.pypeline.cfg')]
	if args.config:
		cfg_files.append(args.config)
	config.read(cfg_files)
	tmpdir = os.path.join(config.get('DEFAULT','tmpscratch'),getuser(),'tmp') if args.tmpdir == None else args.tmpdir
	makedir(tmpdir)

	missing = list_missing_files(args.input)
	if missing:
		raise CatFilesMissingInputError('\n'.join(missing))
	makedir(os.path.dirname(args.output))
	
	if args.MergeSamFiles:
		picard_jar = config.get('picard','MergeSamFiles') if args.picard_jar == None else args.picard_jar
		picard_options = {
			'TMP_DIR': tmpdir,
			'VERBOSITY': 'INFO',
			'QUIET': 'false',
			'VALIDATION_STRINGENCY': 'STRICT',
			'COMPRESSION_LEVEL': '5',
			'MAX_RECORDS_IN_RAM': '1000000',
			'CREATE_INDEX': 'true',
			'CREATE_MD5_FILE': 'true',
			'SORT_ORDER': 'coordinate',
			'ASSUME_SORTED': 'false',
			'MERGE_SEQUENCE_DICTIONARIES': 'false',
			'USE_THREADING': 'true'
		}
		if args.PicardOptions:
			for opt in args.PicardOptions:
				kv = opt.strip().split('=')
				if len(kv) == 2:
					picard_options[kv[0]] = kv[1]
				else:
					raise CatFilePicardOptionsError(opt)
#['java', '-Xmx5g', '-Djava.io.tmpdir=/scratch1/tmp/myourshaw/tmp', '-jar', '/share/apps/myourshaw/picard-current/MergeSamFiles.jar', 'CREATE_INDEX=true', 'TMP_DIR=/scratch1/tmp/myourshaw/tmp', 'VERBOSITY=INFO', 'QUIET=false', 'MERGE_SEQUENCE_DICTIONARIES=false', 'VALIDATION_STRINGENCY=STRICT', 'USE_THREADING=true', 'MAX_RECORDS_IN_RAM=1000000', 'CREATE_MD5_FILE=true', 'ASSUME_SORTED=false', 'SORT_ORDER=coordinate', 'COMPRESSION_LEVEL=5', 'INPUT=/scratch1/tmp/myourshaw/gmd/bams-2011-05-11.HWI-ST430.243/rg/GMD154A_XT_Illumina.2011-05-11.HWI-ST430.243.1.TTAGGCA/GMD154A_XT_Illumina.2011-05-11.HWI-ST430.243.1.TTAGGCA.novoalign.markdup.bam', 'INPUT=/scratch1/tmp/myourshaw/gmd/bams-2011-05-11.HWI-ST430.243/rg/GMD154A_XT_Illumina.2011-05-11.HWI-ST430.243.2.TTAGGCA/GMD154A_XT_Illumina.2011-05-11.HWI-ST430.243.2.TTAGGCA.novoalign.markdup.bam', 'OUTPUT=/scratch1/tmp/myourshaw/gmd/unbarcode-test/cat_test.bam']
		java_args = ['java', '-Xmx5g', '-Djava.io.tmpdir={}'.format(tmpdir), '-jar']
		java_args += [picard_jar]
		java_args += ['{}={}'.format(kv[0], kv[1]) for kv in picard_options.items()]
		java_args += ['INPUT={}'.format(i) for i in args.input]
		java_args += ['OUTPUT={}'.format(args.output)]
		return subprocess.check_call(java_args)

	else:
		with open(args.output, 'wb') as o:
			for i in args.input:
				shutil.copyfileobj(open(i, 'rb'), o)

if __name__ == "__main__": sys.exit(main())
