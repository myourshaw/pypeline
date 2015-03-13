#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from ConfigParser import SafeConfigParser #configparser in python 3
import shutil
import my
import job

class GatkError(Exception): pass


def run(tool, gatk_current=None, email=None, debug=False, synchronous=True, hold_jid=None, **kwargs):
	config = SafeConfigParser()
	config.readfp(open(os.path.join(os.path.dirname(__file__),'pypeline.cfg')))
	cfg_files = [os.path.expanduser('~/.pypeline.cfg')]
	if kwargs.get('config'):
		cfg_files.append(args.kwargs['config'])
	config.read(cfg_files)
	ref = config.get('reference','default_reference')
	gatk_info = {
		'ValidateSamFile': {'defaults': {'MODE': 'VERBOSE', 'REFERENCE_SEQUENCE': ref, 'VALIDATE_INDEX': 'true', 'INPUT': None, 'OUTPUT': None},
			'run_parameters': {'processors': 1, 'memory': '3G'}},
		'FixMateInformation': {'defaults': {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate', 'INPUT': None, 'OUTPUT': None},
			'run_parameters': {'processors': 1, 'memory': '3G'}},
		'MergeSamFiles': {'defaults': {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate', 'ASSUME_SORTED': 'false', 'MERGE_SEQUENCE_DICTIONARIES': 'false', 'USE_THREADING': 'true', 'INPUT': None, 'OUTPUT': None},
			'run_parameters': {'processors': 1, 'memory': '3G'}},
		'MarkDuplicates': {'defaults': {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'REMOVE_DUPLICATES': 'false', 'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP': 8000, 'READ_NAME_REGEX': "\"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*\"", 'OPTICAL_DUPLICATE_PIXEL_DISTANCE': 100, 'METRICS_FILE': None, 'INPUT': None, 'OUTPUT': None},
			'run_parameters': {'processors': 1, 'memory': '3G'}},
		'BuildBamIndex': {'defaults': {'INPUT': None},
			'run_parameters': {'processors': 1, 'memory': '3G'}},
	}
	if not gatk_info.get(tool):
		raise GatkError('unable to run {}'.format(tool))
	out_dir = os.path.dirname(kwargs['OUTPUT']) if kwargs.get('OUTPUT') else  None
	if out_dir:
		my.makedir(out_dir)
	tmp_dir = kwargs['TMP_DIR'] if kwargs.get('TMP_DIR') else my.unique_dir('tmp_picard',out_dir) if out_dir else my.unique_dir('tmp_picard')
	picard_current = picard_current if picard_current else config.get('picard','picard_current')
	jar = os.path.join(picard_current, tool+'.jar')
	picard_options = picard_info[tool]['defaults']
	if kwargs:
		for k in kwargs.keys():
			if isinstance(kwargs[k], (list, tuple)):
				picard_options[k] = ' {}='.format(k).join(kwargs[k])
			else:
				picard_options[k] = kwargs[k]
	for o in picard_options.keys():
		if picard_options[o] == None:
			raise PicardError('missing value for picard {} option {}'.format(tool, o))
	option_str = picard_options = ' '.join(['{}={}'.format(k,picard_options[k]) for k in picard_options.keys()])
	cmd = 'java -Xmx5g -Djava.io.tmpdir={tmp_dir} -jar {jar} TMP_DIR={tmp_dir} VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=1000000 {option_str}'.format(tmp_dir=tmp_dir, jar=jar, option_str=option_str)
	jobName = 'picard_{}{}'.format(tool, '_'+os.path.basename(kwargs['OUTPUT']) if kwargs.get('OUTPUT') else '')
	job_dir = tmp_dir
	qout_dir = tmp_dir
	processors = picard_info[tool]['run_parameters']['processors']
	memory = picard_info[tool]['run_parameters']['memory']
	completionMsg = my.run_job(cmd, jobName, job_dir, qout_dir, email, synchronous=synchronous, processors=processors, memory=memory, hold_jid=hold_jid)
	if synchronous and not debug:
		shutil.rmtree(tmp_dir, ignore_errors=True)
	return completionMsg
